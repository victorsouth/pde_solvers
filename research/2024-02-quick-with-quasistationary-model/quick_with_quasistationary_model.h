#include <random>
#include <algorithm>
using namespace std;

// Определение алиасов типов
using target_var_t = quickest_ultimate_fv_solver_traits<1>::var_layer_data;
using specific_data_t = quickest_ultimate_fv_solver_traits<1>::specific_layer;
using layer_t = composite_layer_t<target_var_t, specific_data_t>;



/// @brief Тесты для солвера quickest_ultimate_fv_solver
/// @brief Тесты для солвера quickest_ultimate_fv_solver
class QuickWithQuasiStationaryModel : public ::testing::Test {
protected:
    // Профиль переменных
    typedef quickest_ultimate_fv_solver_traits<1>::var_layer_data target_var_t;
    typedef quickest_ultimate_fv_solver_traits<1>::specific_layer specific_data_t;

    // Слой: переменных Vars + сколько угодно служебных Specific
    typedef composite_layer_t<target_var_t, specific_data_t> layer_t;
protected:
    /// @brief Параметры трубы
    pipe_properties_t pipe;
    /// @brief Профиль расхода
    vector<double> Q_profile;
    std::unique_ptr<PipeQAdvection> advection_model;
    std::unique_ptr<ring_buffer_t<layer_t>> buffer;
protected:

    /// @brief Подготовка к расчету для семейства тестов
    virtual void SetUp() override {

        simple_pipe_properties simple_pipe;
        simple_pipe.length = 200e3; // тест трубы 700км
        simple_pipe.diameter = 0.514; // тест трубы 700км
        simple_pipe.dx = 100; // тест трубы 700км
        pipe = pipe_properties_t::build_simple_pipe(simple_pipe);
        Q_profile = vector<double>(pipe.profile.getPointCount(), 0.2);
        advection_model = std::make_unique<PipeQAdvection>(pipe, Q_profile);
        buffer = std::make_unique<ring_buffer_t<layer_t>>(2, pipe.profile.getPointCount());
    }
};

TEST_F(QuickWithQuasiStationaryModel, UseCase_Advection)
{
    const auto& x = advection_model->get_grid();
    double dx = x[1] - x[0];
    double v = advection_model->getEquationsCoeffs(0, 0);

    layer_t& prev = buffer->previous();
    layer_t& next = buffer->current();

    prev.vars.cell_double[0] = vector<double>(prev.vars.cell_double[0].size(), 850);
    
    quickest_ultimate_fv_solver solver(*advection_model, *buffer);
    double dt = abs(dx / advection_model->getEquationsCoeffs(0, 0));
    double rho_in = 840; // плотность нефти, закачиваемой на входе трубы при положительном расходе
    double rho_out = 860; // плотность нефти, закачиваемой с выхода трубы при отрицательном расходе
    solver.step(dt, rho_in, rho_out);


    auto& c_new = next.vars.cell_double[0];
}



// Проблемно-ориентированный слой
struct density_viscosity_layer_for_quick {
    quickest_ultimate_fv_solver_traits<1>::var_layer_data density;
    quickest_ultimate_fv_solver_traits<1>::var_layer_data viscosity;
    quickest_ultimate_fv_solver_traits<1>::specific_layer specific;

    density_viscosity_layer_for_quick(size_t cell_count)
        : density(cell_count)
        , viscosity(cell_count)
        , specific(cell_count)
    {
    }

    // Обертка для плотности
    static quickest_ultimate_fv_wrapper<1> get_density_quick_wrapper(density_viscosity_layer_for_quick& layer)
    {
        return quickest_ultimate_fv_wrapper<1>(layer.density, layer.specific);
    }

    static quickest_ultimate_fv_wrapper<1> get_viscosity_quick_wrapper(density_viscosity_layer_for_quick& layer)
    {
        return quickest_ultimate_fv_wrapper<1>(layer.viscosity, layer.specific);
    }
};


TEST_F(QuickWithQuasiStationaryModel, UseCase_Advection_Density_Viscosity)
{
    const auto& x = advection_model->get_grid();
    double dx = x[1] - x[0];
    double v = advection_model->getEquationsCoeffs(0, 0);

    ring_buffer_t<density_viscosity_layer_for_quick> buffer(2, pipe.profile.getPointCount());

    auto& rho_initial = buffer[0].density.cell_double[0];
    auto& viscosity_initial = buffer[0].viscosity.cell_double[0];
    rho_initial = vector<double>(rho_initial.size(), 850); // инициализация начальной плотности
    viscosity_initial = vector<double>(viscosity_initial.size(), 1e-5); // инициализация начальной плотности

    {
        auto density_buffer = buffer.get_custom_buffer(&density_viscosity_layer_for_quick::get_density_quick_wrapper);
        auto& prev = density_buffer.previous();
        auto& next = density_buffer.current();

        prev.vars.cell_double[0] = vector<double>(prev.vars.cell_double[0].size(), 850);
 
        quickest_ultimate_fv_solver solver(*advection_model, density_buffer);
        double dt = abs(dx / v);
        double rho_in = 840; // плотность нефти, закачиваемой на входе трубы при положительном расходе
        double rho_out = 860; // плотность нефти, закачиваемой с выхода трубы при отрицательном расходе
        solver.step(dt, rho_in, rho_out);

        auto& c_new = next.vars.cell_double[0];
    }

    {
        auto viscosity_buffer = buffer.get_custom_buffer(&density_viscosity_layer_for_quick::get_viscosity_quick_wrapper);
        auto& prev = viscosity_buffer.previous();
        auto& next = viscosity_buffer.current();
        
        prev.vars.cell_double[0] = vector<double>(prev.vars.cell_double[0].size(), 1e-5);

        quickest_ultimate_fv_solver solver(*advection_model, viscosity_buffer);
        double dt = abs(dx / v);
        double visc_in = 2e-5; // вязкость нефти, закачиваемой на входе трубы при положительном расходе
        double visc_out = 0.5e-5;; // вязкость нефти, закачиваемой с выхода трубы при отрицательном расходе
        solver.step(dt, visc_in, visc_out);

        auto& c_new = next.vars.cell_double[0];
    }
    buffer.advance(+1);
    auto& curr = buffer[0];
}

/// @brief Уравнение трубы для задачи PQ
class pipe_model_PQ_cell_parties_t : public ode_t<1>
{
public:
    using ode_t<1>::equation_coeffs_type;
    using ode_t<1>::right_party_type;
    using ode_t<1>::var_type;
protected:
    const vector<double>& rho_profile;
    const vector<double>& nu_profile;
    const pipe_properties_t& pipe;
    const double flow;
    const int solver_direction;
public:
    /// @brief Констуктор уравнения трубы
    /// @param pipe Ссылка на сущность трубы
    /// @param oil Ссылка на сущность нефти
    /// @param flow Объемный расход
    /// @param solver_direction Направление расчета по Эйлеру, должно обязательно совпадать с параметром солвера Эйлера
    pipe_model_PQ_cell_parties_t(pipe_properties_t& pipe, vector<double>& rho_profile, vector<double>& nu_profile, double flow,
        int solver_direction)
        : pipe(pipe)
        , rho_profile(rho_profile)
        , nu_profile(nu_profile)
        , flow(flow)
        , solver_direction(solver_direction)
    {

    }

    /// @brief Возвращает известную уравнению сетку
    virtual const vector<double>& get_grid() const override {
        return pipe.profile.coordinates;
    }

    /// @brief Возвращает значение правой части ДУ
    /// @param grid_index Обсчитываемый индекс расчетной сетки
    /// @param point_vector Начальные условия
    /// @return Значение правой части ДУ в точке point_vector
    virtual right_party_type ode_right_party(
        size_t grid_index, const var_type& point_vector) const override
    {
       
        /// Обработка индекса в случае расчетов на границах трубы
        /// Чтобы не выйти за массив высот, будем считать dz/dx в соседней точке
        //grid_index = grid_index == 0 ? grid_index + 1 : grid_index;
        //grid_index = grid_index == pipe.profile.heights.size() - 1 ? grid_index - 1 : grid_index;

        size_t reo_index = solver_direction == +1
            ? grid_index
            : grid_index - 1;

        double rho = rho_profile[reo_index];
        double S_0 = pipe.wall.getArea();
        double v = flow / (S_0);
        double Re = v * pipe.wall.diameter / nu_profile[reo_index];
        double lambda = pipe.resistance_function(Re, pipe.wall.relativeRoughness());
        double tau_w = lambda / 8 * rho * v * abs(v);

        double height_derivative = pipe.profile.get_height_derivative(grid_index, solver_direction);
        double result = -4*tau_w / pipe.wall.diameter - rho * M_G * height_derivative;
        return result;
    }
};

