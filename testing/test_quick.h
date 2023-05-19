#pragma once

using namespace std;

/// @brief Описание типов данных для метода конечных объемов на основе upstream differencing 
template <size_t Dimension>
struct upstream_fv_solver_traits
{
    typedef templated_layer<0, Dimension/*переменные - ячейки*/, 0, 0, 0, 0> var_layer_data;
    typedef templated_layer<Dimension /*потоки F*/, 0,
        0, 0,
        0, 0> specific_layer;
};

/// @brief Описание типов данных для метода конечных объемов на основе QUICK 
template <size_t Dimension>
struct quick_fv_solver_traits
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
        : upstream_fv_solver(pde, buffer.previous(), buffer.current())
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
            double v_cell = pde.getEquationsCoeffs(cell, u); // не совсем корректно, скорость в ячейке берется из скорости на ее левой границе
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

double quick_border_approximation(double U_L, double U_C, double U_R)
{
    double Ub_linear = (U_C + U_R) / 2;
    double Ub_cubic_correction = -(U_L + U_R - 2 * U_C) / 8;
    return Ub_linear + Ub_cubic_correction;
}

/// @brief Солвер на основе QUICK, только для размерности 1!
/// [Leonard 1979]
class quick_fv_solver {
public:
    typedef typename quick_fv_solver_traits<1>::var_layer_data var_layer_data;
    typedef typename quick_fv_solver_traits<1>::specific_layer specific_layer;
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
    quick_fv_solver(pde_t<1>& pde,
        custom_buffer_t<composite_layer_t<var_layer_data, specific_layer>>& buffer)
        : quick_fv_solver(pde, buffer.previous(), buffer.current())
    {}

    /// @brief Конструктор для простых слоев - 
    /// когда в слое только один блок целевых переменных и один блок служебных данных
    /// @param pde ДУЧП
    /// @param prev Предыдущий слой (уже рассчитанный)
    /// @param curr Следующий (новый), для которого требуется сделать расчет
    quick_fv_solver(pde_t<1>& pde,
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

        double v_pipe = pde.getEquationsCoeffs(0, U[0]);//не совсем корректно, скорость в ячейке берется из скорости на ее левой границе
        // Расчет потоков на границе по правилу QUICK
        if (v_pipe >= 0) {
            for (size_t cell = 0; cell < U.size(); ++cell) {
                size_t right_border = cell + 1;
                double Vb = v_pipe; // предположили, что скорость на границе во всех точках трубы одна и та же
                double Ub;
                if (cell == 0) {
                    Ub = quick_border_approximation(U[cell], U[cell], U[cell + 1]); // костыль U_L = U_C
                }
                else if (cell == U.size() - 1) {
                    Ub = quick_border_approximation(U[cell - 1], U[cell], U[cell]); // костыль U_R = U_C
                }
                else {
                    Ub = quick_border_approximation(U[cell - 1], U[cell], U[cell + 1]); // честный расчет
                }
                F[right_border] = Ub * Vb;
            }
        }
        else {
            for (size_t cell = 0; cell < U.size(); ++cell) {
                size_t left_border = cell;
                double Vb = v_pipe; // предположили, что скорость на границе во всех точках трубы одна и та же
                double Ub;
                if (cell == 0) {
                    Ub = quick_border_approximation(U[cell + 1], U[cell], U[cell]); // костыль U_L = U_C
                }
                else if (cell == U.size() - 1) {
                    Ub = quick_border_approximation(U[cell], U[cell], U[cell - 1]); // костыль U_R = U_C
                }
                else {
                    Ub = quick_border_approximation(U[cell + 1], U[cell], U[cell - 1]); // честный расчет
                }
                F[left_border] = Ub * Vb;
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

/// @brief Тесты для солвера quick_fv_solver
class QUICK : public ::testing::Test {
protected:
    // Профиль переменных
    typedef quick_fv_solver_traits<1>::var_layer_data target_var_t;
    typedef quick_fv_solver_traits<1>::specific_layer specific_data_t;

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

        Q = vector<double>(pipe.profile.getPointCount(), 0.5);
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

/// @brief Пример вывода в файл через
TEST_F(UpstreamDifferencing, UseCaseStepDensity)
{
    string path = prepare_test_folder();
    std::ofstream output(path + "output.csv");

    vector<double> rho_in{ 860, 860, 860, 860, 860, 860, 860, 860, 860, 860 };
    double rho_out = 870;
    double dt = 60; // 1 минута
    double t = 0;

    for (size_t index = 0; index < rho_in.size(); ++index) {
        if (index == 0) {
            layer_t& prev = buffer->previous();
            prev.vars.print(t, output);
        }

        t += dt;
       
        upstream_fv_solver solver(*advection_model, *buffer);
        solver.step(dt, rho_in[index], rho_out);

        layer_t& next = buffer->current();
        next.vars.print(t, output);

        buffer->advance(+1);

    }
    output.flush();
    output.close();
}

/// @brief Пример вывода в файл через
TEST_F(QUICK, UseCaseStepDensity)
{
    string path = prepare_test_folder();


    double rho_in = 860;
    double rho_out = 870;
    double T = 50000; // период моделирования

    const auto& x = advection_model->get_grid();
    double dx = x[1] - x[0];
    double v = advection_model->getEquationsCoeffs(0, 0);
    double dt_ideal = dx / v;


    for (double Cr = 0.05; Cr < 1.01; Cr += 0.05) {
        advection_model = std::make_unique<PipeQAdvection>(pipe, Q);
        buffer = std::make_unique<custom_buffer_t<layer_t>>(2, pipe.profile.getPointCount());

        layer_t& prev = buffer->previous();
        prev.vars.cell_double[0] = vector<double>(prev.vars.cell_double[0].size(), 850);

        double t = 0; // текущее время
        //double dt = 60; // 1 минута
        double dt = Cr * dt_ideal; // время в долях от Куранта

        std::stringstream filename;
        filename << path << "output Cr=" << Cr << ".csv";
        std::ofstream output(filename.str());


        size_t N = static_cast<int>(T / dt);
        for (size_t index = 0; index < N; ++index) {
            if (index == 0) {
                layer_t& prev = buffer->previous();
                prev.vars.print(t, output);
            }

            t += dt;

            quick_fv_solver solver(*advection_model, *buffer);
            solver.step(dt, rho_in, rho_out);

            layer_t& next = buffer->current();
            next.vars.print(t, output);

            buffer->advance(+1);

        }
        output.flush();
        output.close();
    }

}

/// @brief Разработка метода QUICK по [Leonard 1979]
TEST(QUICK, Develop)
{
    // Профиль переменных
    typedef templated_layer<0, 1> target_var_t;
    typedef templated_layer<1, 0> specific_data_t;

    // Слой: переменных Vars + сколько угодно служебных Specific
    typedef composite_layer_t<target_var_t, specific_data_t> layer_t;

    // Упрощенное задание трубы - 50км, с шагом разбиения для расчтной сетки 1км, диаметром 700мм
    simple_pipe_properties simple_pipe;
    simple_pipe.length = 50e3;
    simple_pipe.diameter = 0.7;
    simple_pipe.dx = 1000;

    PipeProperties pipe = PipeProperties::build_simple_pipe(simple_pipe);

    vector<double> Q(pipe.profile.getPointCount(), -0.5); // задаем по трубе расход 0.5 м3/с
    PipeQAdvection advection_model(pipe, Q);

    //Текущий и предыдущий слой, каждый из которых представляет собой composite_layer (Var+Specific)
    custom_buffer_t<layer_t> buffer(2, 10);

    //Получение текущего/предыдущего слоя
    layer_t& prev = buffer.previous();
    layer_t& next = buffer.current();

    double rho_in = 860;
    double rho_out = 870;

    prev.vars.cell_double[0] = vector<double>(prev.vars.cell_double[0].size(), 850);

    auto& next_spec = std::get<0>(next.specific);
    auto& F = next_spec.point_double[0]; // потоки на границах ячеек
    auto& U = prev.vars.cell_double[0];
    auto& U_new = next.vars.cell_double[0]; //

    double dt = 60; // 1 минута

    // Расчет потоков на границе на основе граничных условий
    double v_in = advection_model.getEquationsCoeffs(0, U[0]);
    double v_out = advection_model.getEquationsCoeffs(F.size() - 1, U[U.size() - 1]);

    if (v_in >= 0) {
        F.front() = v_in * rho_in;
    }
    if (v_out <= 0) {
        F.back() = v_out * rho_out;
    }

    double v_pipe = advection_model.getEquationsCoeffs(0, U[0]);//не совсем корректно, скорость в ячейке берется из скорости на ее левой границе

    if (v_pipe >= 0) {
        // Расчет потоков на границе по правилу QUICK
        for (size_t cell = 0; cell < U.size(); ++cell) {
            size_t right_border = cell + 1;
            double Vb = v_pipe; // предположили, что скорость на границе во всех точках трубы одна и та же
            double Ub;
            if (cell == 0) {
                Ub = quick_border_approximation(U[cell], U[cell], U[cell + 1]); // костыль U_L = U_C
            }
            else if (cell == U.size() - 1) {
                Ub = quick_border_approximation(U[cell - 1], U[cell], U[cell]); // костыль U_R = U_C
            }
            else {
                Ub = quick_border_approximation(U[cell - 1], U[cell], U[cell + 1]); // честный расчет
            }
            F[right_border] = Ub * Vb;
        }
    }
    else {
        for (size_t cell = 0; cell < U.size(); ++cell) {
            size_t left_border = cell;
            double Vb = v_pipe; // предположили, что скорость на границе во всех точках трубы одна и та же
            double Ub;
            if (cell == 0) {
                Ub = quick_border_approximation(U[cell + 1], U[cell], U[cell]); // костыль U_L = U_C
            }
            else if (cell == U.size() - 1) {
                Ub = quick_border_approximation(U[cell], U[cell], U[cell - 1]); // костыль U_R = U_C
            }
            else {
                Ub = quick_border_approximation(U[cell + 1], U[cell], U[cell - 1]); // честный расчет
            }
            F[left_border] = Ub * Vb;
        }

    }

    const auto& grid = advection_model.get_grid();
    for (size_t cell = 0; cell < U.size(); ++cell) {
        double dx = grid[cell + 1] - grid[cell]; // ячейки обычно одинаковой длины, но мало ли..
        U_new[cell] = U[cell] + dt / dx * ((F[cell] - F[cell + 1]));
    }

    //Движение на слой вперед 
    //buffer.advance(+1);
}