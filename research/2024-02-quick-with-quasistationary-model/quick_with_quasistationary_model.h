#include <random>
#include <algorithm>
using namespace std;

// Определение алиасов типов
using target_var_t = quickest_ultimate_fv_solver_traits<1>::var_layer_data;
using specific_data_t = quickest_ultimate_fv_solver_traits<1>::specific_layer;
using layer_t = composite_layer_t<target_var_t, specific_data_t>;

struct my_task_parameters
{
    /// @brief Плотность жидкости, (кг/м^3)
    double rho;
    /// @brief Кинематическая вязкость, (сСТ)
    double nu;
    /// @brief Давление в начале участка, (Па)
    double p_0;
    /// @brief Давление в конце участка, (Па)
    double p_L;
    /// @brief Расход жидкости, (м^3/с)
    double Q;
    /// @brief Профиль плотностей по всей трубе, (кг/м^3)
    vector<double> rho_profile;
    /// @brief Профиль вязкостей по всей трубе, (сСТ)
    vector<double> nu_profile;
    my_task_parameters(double rho, double nu, double p_0, double p_L, double Q, const vector<double> rho_profile, const vector<double> nu_profile) :
        rho{ rho }, nu{ nu }, p_0{ p_0 }, p_L{ p_L }, Q{ Q }, rho_profile{ rho_profile }, nu_profile{ nu_profile }
    {
    }
};



/// @brief Функция расчета скорости из расхода
/// @param Q Расход, (м^3/с)
/// @param internal_diameter Внутренний диаметр, (м)
/// @return Скорость, (м/с)
double calc_speed(double Q, double internal_diameter)
{
    double speed = (4 * Q) / (pow(internal_diameter, 2) * M_PI);
    return speed;
}

/// @brief Функция расчета расхода из скорости
/// @param speed Скорость, (м/с)
/// @param internal_diameter Внутренний диаметр, (м) 
/// @return Расход, (м^3/с)
double calc_flow(double speed, double internal_diameter)
{
    double flow = (pow(internal_diameter, 2) * M_PI * speed) / 4;
    return flow;
}

double linear_interpolator(vector<double> original_time, vector<double> original_value, double new_time_step)
{
    size_t index1 = 0;
    size_t index2 = 1;
    while (original_time[index2] < new_time_step)
    {
        ++index1;
        ++index2;
    }
    double t1 = original_time[index1];
    double t2 = original_time[index2];
    double value1 = original_value[index1];
    double value2 = original_value[index2];
    return value1 + (value2 - value1) * (new_time_step - t1) / (t2 - t1);
}

void print_data_to_csv(const vector<double>& time_rho,
    const vector<double>& rho_and_nu_in_0,
    const vector<double>& time_nu,
    const vector<double>& rho_and_nu_in_1,
    const vector<double>& time_p_in,
    const vector<double>& p_in,
    const vector<double>& time_p_out,
    const vector<double>& p_out,
    const vector<double>& time_Q,
    const vector<double>& Q,
    const string& filename)
{

    // Определяем максимальную длину вектора
    size_t maxLength = max({ time_rho.size(), rho_and_nu_in_0.size(), time_nu.size(), rho_and_nu_in_1.size(),
                                      time_p_in.size(), p_in.size(), time_p_out.size(), p_out.size(), Q.size(), time_Q.size() });

    // Открываем файл для записи
    ofstream file(filename);

    // Проверяем, открыт ли файл успешно
    if (file.is_open()) {
        // Записываем заголовки столбцов
        file << "Time Density; Density; Time Viscosity; Viscosity; Time Pressure In; Pressure In; Time Pressure Out; Pressure Out; Time Flow Rate; Flow Rate\n";

        // Записываем данные из векторов
        for (size_t i = 0; i < maxLength; ++i) {
            // Если индекс находится в пределах длины вектора, записываем значение, иначе записываем пустую ячейку
            file << (i < time_rho.size() ? to_string(time_rho[i]) : "") << ";"
                << (i < rho_and_nu_in_0.size() ? to_string(rho_and_nu_in_0[i]) : "") << ";"
                << (i < time_nu.size() ? to_string(time_nu[i]) : "") << ";"
                << (i < rho_and_nu_in_1.size() ? to_string(rho_and_nu_in_1[i]) : "") << ";"
                << (i < time_p_in.size() ? to_string(time_p_in[i]) : "") << ";"
                << (i < p_in.size() ? to_string(p_in[i]) : "") << ";"
                << (i < time_p_out.size() ? to_string(time_p_out[i]) : "") << ";"
                << (i < p_out.size() ? to_string(p_out[i]) : "") << ";"
                << (i < time_Q.size() ? to_string(time_Q[i]) : "") << ";"
                << (i < Q.size() ? to_string(Q[i]) : "") << "\n";
        }
        // Закрываем файл
        file.close();
    }
}


void print_layers(const double dt,
    const vector<double>& layer,
    const string& filename)
{
    ofstream  file(filename, ios::app);
    if (file.is_open()) {
        file << to_string(dt) << ";";
        for (int j = 0; j < layer.size(); j++)
        {
            file << to_string(layer[j]) << ";";
        }
        file << "\n";
        file.close();
    }
}

void clear_directory(const filesystem::path& dir_path) {
    for (const auto& entry : filesystem::directory_iterator(dir_path)) {
        filesystem::remove_all(entry.path());
    }
}


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
template<typename VarLayerType, typename SpecificLayerType>
struct density_viscosity_layer_for_quick {
    VarLayerType density;
    VarLayerType viscosity;
    SpecificLayerType specific;

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

    ring_buffer_t<density_viscosity_layer_for_quick<target_var_t, specific_data_t>> buffer(2, pipe.profile.getPointCount());

    auto& prev_case = buffer.get_custom_buffer(&density_viscosity_layer_for_quick<target_var_t, specific_data_t>::get_density_quick_wrapper).previous();

    auto& rho_initial = buffer[0].density.cell_double[0];
    auto& viscosity_initial = buffer[0].viscosity.cell_double[0];
    rho_initial = vector<double>(rho_initial.size(), 850); // инициализация начальной плотности
    viscosity_initial = vector<double>(viscosity_initial.size(), 1e-5); // инициализация начальной плотности

    {
        auto density_buffer_wrapped = buffer.get_custom_buffer(&density_viscosity_layer_for_quick<target_var_t, specific_data_t>::get_density_quick_wrapper);
        auto& prev = density_buffer_wrapped.previous();
        auto& next = density_buffer_wrapped.current();

        prev.vars.cell_double[0] = vector<double>(prev.vars.cell_double[0].size(), 850);
        //ring_buffer_t<layer_t> density_buffer_wrapped = buffer.get_custom_buffer(&density_viscosity_layer_for_quick<target_var_t, specific_data_t>::get_density_quick_wrapper);
        quickest_ultimate_fv_solver solver(*advection_model, density_buffer_wrapped);
        double dt = abs(dx / v);
        double rho_in = 840; // плотность нефти, закачиваемой на входе трубы при положительном расходе
        double rho_out = 860; // плотность нефти, закачиваемой с выхода трубы при отрицательном расходе
        solver.step(dt, rho_in, rho_out);

        auto& c_new = next.vars.cell_double[0];
    }
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

class newton_solver_PP_with_euler : public fixed_system_t<1>
{
    using fixed_system_t<1>::var_type;
protected:
    pipe_properties_t& pipe;
    my_task_parameters& task;
    double Q_approx;
    const int solver_direction;
public:
    newton_solver_PP_with_euler(pipe_properties_t& pipe, my_task_parameters& task, double Q_approx, const int solver_direction)
        : pipe(pipe)
        , task(task)
        , Q_approx(Q_approx)
        , solver_direction(solver_direction)
    {

    }
    /// @brief Задание функции невязок
    /// @param x Искомый расход
    /// @return Функция невязок
    var_type residuals(const var_type& x)
    {
        Q_approx = x; // во временной структуре используем Q для нашего уравнения невязки, эта Q будет идти в солвер
        // Объявляем переменную класса солвера Эйлером
        pipe_model_PQ_cell_parties_t pipeModel(pipe, task.rho_profile, task.nu_profile, Q_approx, solver_direction);

        solve_euler_corrector<1>(pipeModel, solver_direction, task.p_0, &p_profile);
        p_profile.pop_back();
        return (p_profile.back() - task.p_L);
    }
private:
    vector<double> p_profile;
};