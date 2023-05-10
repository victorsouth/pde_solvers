#pragma once


/// @brief Описание типов данных для метода конечных объемов на основе upstream differencing 
template <size_t Dimension>
struct upstream_fv_solver_traits
{
    typedef templated_layer<0, Dimension/*переменные - ячейки*/, 0, 0, 0, 0> var_layer_data;
    typedef templated_layer<Dimension /*потоки F*/, 0,
        0, 0,
        0, 0> specific_layer;
};

/// @brief Солвер на основе upstream differencing, только для размерности 1!
/// [Leonard 1979]
class upstream_fv_solver {
public:
    typedef typename upstream_fv_solver_traits<1>::var_layer_data var_layer_data;
    typedef typename upstream_fv_solver_traits<1>::specific_layer specific_layer;
    typedef typename fixed_system_types<1>::matrix_type matrix_type;
    typedef typename fixed_system_types<1>::var_type vector_type;
protected:
    /// @brief ДУЧП
    pde_t<1>& pde;
    /// @brief Сетка, полученная от ДУЧП
    const vector<double>& grid;
    /// @brief Количество точек сетки
    const size_t n;
    /// @brief Предыдущий слой переменных
    const var_layer_data& prev_vars;
    /// @brief Новый (рассчитываемый) слой переменных
    var_layer_data& curr_vars;
    /// @brief Предыдущий специфический слой (сейчас не нужен! нужен ли в будущем?)
    const specific_layer& prev_spec;
    /// @brief Текущий специфический слой
    specific_layer& curr_spec;
public:
    /// @brief Конструктор для буфера в котором простой слой:
    /// когда в слое только один блок целевых переменных и один блок служебных данных
    /// Из буфера берется current() и previous()
    /// @param pde ДУЧП
    /// @param buffer Буфер слоев
    upstream_fv_solver(pde_t<1>& pde,
        custom_buffer_t<composite_layer_t<var_layer_data, specific_layer>>& buffer)
        : upstream_fv_solver(pde, buffer.current(), buffer.previous())
    {}
    
    /// @brief Конструктор для простых слоев - 
    /// когда в слое только один блок целевых переменных и один блок служебных данных
    /// @param pde ДУЧП
    /// @param prev Предыдущий слой (уже рассчитанный)
    /// @param curr Следующий (новый), для которого требуется сделать расчет
    upstream_fv_solver(pde_t<1>& pde,
        const composite_layer_t<var_layer_data, specific_layer>& prev,
        composite_layer_t<var_layer_data, specific_layer>& curr)
        : pde(pde)
        , grid(pde.get_grid())
        , n(pde.get_grid().size())
        , prev_vars(prev.vars)
        , curr_vars(curr.vars)
        , prev_spec(std::get<0>(prev.specific))
        , curr_spec(std::get<0>(curr.specific))
    {

    }
    /// @brief Расчет шага
    /// @param dt Заданный период времени
    /// @param u_in Левое граничное условие
    /// @param u_out Правое граничное условие
    void step(double dt, double u_in, double u_out) {
        auto& F = curr_spec.point_double[0]; // потоки на границах ячеек
        const auto& U = prev_vars.cell_double[0];
        auto& U_new = curr_vars.cell_double[0]; 

        // Расчет потоков на границе на основе граничных условий
        double v_in = pde.getEquationsCoeffs(0, U[0]);
        double v_out = pde.getEquationsCoeffs(F.size() - 1, U[U.size() - 1]);
        if (v_in >= 0) {
            F.front() = v_in * u_in;
        }
        if (v_out <= 0) {
            F.back() = v_out * u_out;
        }


        // Расчет потоков на границе по правилу донорской ячейки
        for (size_t cell = 0; cell < U.size(); ++cell) {
            double u = U[cell];
            double v_cell = pde.getEquationsCoeffs(cell, u);//не совсем корректно, скорость в ячейке берется из скорости на ее левой границе
            if (v_cell > 0) {
                size_t right_point = cell + 1;
                double v_right = pde.getEquationsCoeffs(right_point, u);
                F[right_point] = u * v_right;
            }
            else {
                size_t left_point = cell;
                double v_left = pde.getEquationsCoeffs(left_point, u);
                F[left_point] = u * v_left;
            }
        }

        for (size_t cell = 0; cell < U.size(); ++cell) {
            double dx = grid[cell + 1] - grid[cell]; // ячейки обычно одинаковой длины, но мало ли..
            U_new[cell] = U[cell] + dt / dx * ((F[cell] - F[cell + 1]));
        }
    }
};

/// @brief Тесты для солвера upstream_fv_solver
class UpstreamDifferencing : public ::testing::Test {
protected:
    // Профиль переменных
    typedef upstream_fv_solver_traits<1>::var_layer_data target_var_t;
    typedef upstream_fv_solver_traits<1>::specific_layer specific_data_t;

    // Слой: переменных Vars + сколько угодно служебных Specific
    typedef composite_layer_t<target_var_t, specific_data_t> layer_t;
protected:
    /// @brief Параметры трубы
    PipeProperties pipe;
    /// @brief Профиль расхода
    vector<double> Q;
    std::unique_ptr<PipeQAdvection> advection_model;
    std::unique_ptr<custom_buffer_t<layer_t>> buffer;
protected:
    
    /// @brief Подготовка к расчету для семейства тестов
    virtual void SetUp() override {
        // Упрощенное задание трубы - 50км, с шагом разбиения для расчтной сетки 1км, диаметром 700мм
        simple_pipe_properties simple_pipe;
        simple_pipe.length = 50e3;
        simple_pipe.diameter = 0.7;
        simple_pipe.dx = 1000;
        pipe = PipeProperties::build_simple_pipe(simple_pipe);

        Q = vector<double> (pipe.profile.getPointCount(), 0.5);
        advection_model = std::make_unique<PipeQAdvection>(pipe, Q);
        buffer = std::make_unique<custom_buffer_t<layer_t>>(2, pipe.profile.getPointCount());

        layer_t& prev = buffer->previous();
        prev.vars.cell_double[0] = vector<double>(prev.vars.cell_double[0].size(), 850);
    }
};

/// @brief Проверка, правильно ли учитывается инверсия потока
TEST_F(UpstreamDifferencing, CanConsiderFlowSwap) {
    Q = vector<double>(pipe.profile.getPointCount(), -0.5);

    //Получение текущего/предыдущего слоя
    layer_t& prev = buffer->previous();
    layer_t& next = buffer->current();

    double rho_in = 860;
    double rho_out = 870;
    double dt = 60; // 1 минута

    upstream_fv_solver solver(*advection_model, prev, next);
    solver.step(dt, rho_in, rho_out);

    const auto& rho_prev = prev.vars.cell_double[0];
    const auto& rho_curr = next.vars.cell_double[0];
    ASSERT_GT(rho_curr.back(), rho_prev.back()); // плотность в конце выросла
    ASSERT_NEAR(rho_curr.front(), rho_prev.front(), 1e-8); // плотность в начале не изменилась
}

/// @brief Разработка метода прямых разностей по [Leonard 1979]
TEST_F(UpstreamDifferencing, Develop)
{
    //Получение текущего/предыдущего слоя
    layer_t& prev = buffer->previous();
    layer_t& next = buffer->current();

    double rho_in = 860;
    double rho_out = 870;
    double dt = 60; // 1 минута

    upstream_fv_solver solver(*advection_model, prev, next);
    solver.step(dt, rho_in, rho_out);

}


