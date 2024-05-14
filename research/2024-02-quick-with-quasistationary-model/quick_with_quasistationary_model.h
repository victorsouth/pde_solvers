#pragma once

#include <pde_solvers/timeseries.h>

/// @brief Проблемно-ориентированный слой для расчета методом конечных объемов 
struct density_viscosity_cell_layer {
    /// @brief Профиль плотности
    std::vector<double> density;
    /// @brief Профиль вязкости
    std::vector<double> viscosity;
    /// @brief Профиль давления
    std::vector<double> pressure;
    /// @brief Дифференциальный профиль давления
    std::vector<double> pressure_delta;
    /// @brief Профиль вспомогательных расчетов для метода конечных объемов (и для вязкости, и для плотности)
    quickest_ultimate_fv_solver_traits<1>::specific_layer specific;
    /// @brief Инициализация профилей
    /// @param point_count Количество точек
    density_viscosity_cell_layer(size_t point_count)
        : density(point_count - 1)
        , viscosity(point_count - 1)
        , specific(point_count)
        , pressure(point_count)
        , pressure_delta(point_count)
    {}

    // @brief Подготовка плотности для расчета методом конечных объемов 
    static quickest_ultimate_fv_wrapper<1> get_density_wrapper(density_viscosity_cell_layer& layer)
    {
        return quickest_ultimate_fv_wrapper<1>(layer.density, layer.specific);
    }
    /// @brief Подготовка вязкости для расчета методом конечных объемов 
    static quickest_ultimate_fv_wrapper<1> get_viscosity_wrapper(density_viscosity_cell_layer& layer)
    {
        return quickest_ultimate_fv_wrapper<1>(layer.viscosity, layer.specific);
    }
};
/// @brief Проблемно-ориентированный слой для расчета методом характеристик
struct density_viscosity_layer_moc{
    /// @brief Профиль плотности
    std::vector<double> density;
    /// @brief Профиль вязкости
    std::vector<double> viscosity;
    /// @brief Профиль давления
    std::vector<double> pressure;
    /// @brief Дифференциальный профиль давления
    std::vector<double> pressure_delta;
    /// @brief Профиль вспомогательных расчетов для МХ (и для вязкости, и для плотности)
    moc_solver<1>::specific_layer specific;
    /// @brief Конструктор на заданное количество точек
    density_viscosity_layer_moc(size_t point_count)
        : density(point_count)
        , viscosity(point_count)
        , specific(point_count)
        , pressure(point_count)
        , pressure_delta(point_count)
    {

    }
    /// @brief Подготовка плотности для расчета методом характеристик
    /// Оборачивает профиль плотности и вспомогательный расчет МХ в обертку для МХ
    static moc_layer_wrapper<1> get_density_wrapper(density_viscosity_layer_moc& layer)
    {
        return moc_layer_wrapper<1>(layer.density, layer.specific);
    }
    /// @brief Подготовка вязкости для расчета методом характеристик
    static moc_layer_wrapper<1> get_viscosity_wrapper(density_viscosity_layer_moc& layer)
    {
        return moc_layer_wrapper<1>(layer.viscosity, layer.specific);
    }
};

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
    size_t flag_for_points;
public:
    /// @brief Констуктор уравнения трубы
    /// @param pipe Ссылка на сущность трубы
    /// @param oil Ссылка на сущность нефти
    /// @param flow Объемный расход
    /// @param solver_direction Направление расчета по Эйлеру, должно обязательно совпадать с параметром солвера Эйлера
    pipe_model_PQ_cell_parties_t(const pipe_properties_t& pipe, const vector<double>& rho_profile, const vector<double>& nu_profile, double flow,
        int solver_direction, size_t flag_for_points = 0)
        : pipe(pipe)
        , rho_profile(rho_profile)
        , nu_profile(nu_profile)
        , flow(flow)
        , solver_direction(solver_direction)
        , flag_for_points(flag_for_points)
    {}

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
        size_t reo_index = grid_index;
        if (pipe.profile.getPointCount() == rho_profile.size())
        { 
            if (solver_direction == +1)
                reo_index += 1;
            else
                reo_index -= 1;
        }
        else
        { 
            reo_index = solver_direction == +1
            ? grid_index
            : grid_index - 1;
        }
        double rho = rho_profile[reo_index];
        double S_0 = pipe.wall.getArea();
        double v = flow / (S_0);
        double Re = v * pipe.wall.diameter / nu_profile[reo_index];
        double lambda = pipe.resistance_function(Re, pipe.wall.relativeRoughness());
        double tau_w = lambda / 8 * rho * v * abs(v);
        if (reo_index == 1999)
            double stop = 0;
        double height_derivative = pipe.profile.get_height_derivative(grid_index, solver_direction);
        double result = -4 * tau_w / pipe.wall.diameter - rho * M_G * height_derivative;
        return result;
    }
};

/// @brief Формирует имя файл для результатов исследования разных численных метов
/// @tparam Solver Класс солвера
/// @param path Путь, в котором формируется файл
/// @return Путь к файлу в заивисимости от указанного класса солвера
static std::string get_courant_research_filename_for_qsm(const string& path, const string& layer_name)
{
    std::stringstream filename;
    filename << path << "output " << layer_name << ".csv";
    return filename.str();
}

/// @brief Структура, созданная для хранения в себе начального профиля давлений и буфера с расчетными данными
template <typename Layer>
struct isothermal_quasistatic_task_buffer_t {
    /// @brief Изначальный профиль давления
    vector<double> pressure_initial;
    /// @brief Буфер профилей давления, плотности, вязкости
    ring_buffer_t<Layer> buffer;
    isothermal_quasistatic_task_buffer_t(size_t point_count)
        : pressure_initial(point_count)
        , buffer(2, point_count)
    {}
};

/// @brief Структура, содержащая в себе начальные условия задачи PQ
struct isothermal_quasistatic_task_boundaries_t {
    /// @brief Изначальный объемный расход
    double volumetric_flow;
    /// @brief Изначальное давление на входе
    double pressure_in;
    /// @brief Изначальная плотность на входе
    double density;
    /// @brief Изначальная вязкость на входе
    double viscosity;

    isothermal_quasistatic_task_boundaries_t() = default;

    isothermal_quasistatic_task_boundaries_t(const vector<double>& values) {
        volumetric_flow = values[0];
        pressure_in = values[1];
        density = values[2];
        viscosity = values[3];

        //boundaries.volumetric_flow = values_in_time_model[0];
        //boundaries.pressure_in = values_in_time_model[1];
        //boundaries.density = values_in_time_model[2];
        //boundaries.viscosity = values_in_time_model[3];

    }

    static isothermal_quasistatic_task_boundaries_t default_values() {
        isothermal_quasistatic_task_boundaries_t result;
        result.volumetric_flow = 0.2;
        result.pressure_in = 6e6;
        result.density = 850;
        result.viscosity = 15e-6;
        return result;
    }
};
/// @brief Проблемно-ориентированный класс для задачи PQ


/// @brief Расчетная задача (task) для гидравлического изотермического 
/// квазистационарного расчета в условиях движения партий с разной плотностью и вязкостью
/// Расчет партий делается методом характеристик или Quickest-Ultimate
/// @tparam Layer Тип слоя, содержащего профили плотности, вязкости, давления
/// Для партий методом характеристик = density_viscosity_layer_moc, для партий методом Quickest-Ultimate =density_viscosity_cell_layer
/// @tparam Solver Тип солвера партий (moc_solver или quickest_ultimate_fv_solver)
template <typename Layer, typename Solver>
class isothermal_quasistatic_task_t {
    pipe_properties_t pipe;
    isothermal_quasistatic_task_buffer_t<Layer> buffer;

public:
    isothermal_quasistatic_task_t(const pipe_properties_t& pipe)
        : pipe(pipe)
        , buffer(pipe.profile.getPointCount())
    {

    }

    void solve(const isothermal_quasistatic_task_boundaries_t& initial_conditions)
    {
        size_t n = pipe.profile.getPointCount();

        // Инициализация реологии
        auto& current = buffer.buffer.current();

        // Инициализация начального профиля плотности (не важно, ячейки или точки)
        for (double& density : current.density) {
            density = initial_conditions.density;
        }
        // Инициализация начального профиля вязкости (не важно, ячейки или точки)
        for (double& viscosity : current.viscosity) {
            viscosity = initial_conditions.viscosity;
        }

        // Начальный гидравлический расчет
        int euler_direction = +1;
        pipe_model_PQ_cell_parties_t pipeModel(pipe, current.density, current.viscosity, initial_conditions.volumetric_flow, euler_direction);
        solve_euler<1>(pipeModel, euler_direction, initial_conditions.pressure_in, &current.pressure);

        buffer.pressure_initial = current.pressure; // Получаем изначальный профиль давлений
    }
public:
    double get_time_step_assuming_max_speed(double v_max) const {
        const auto& x = pipe.profile.coordinates;
        double dx = x[1] - x[0]; // Шаг сетки
        double dt = abs(dx / v_max); // Постоянный шаг по времени для Куранта = 1
        return dt;
    }
private:
    void make_rheology_step(double dt, const isothermal_quasistatic_task_boundaries_t& boundaries) {
        size_t n = pipe.profile.getPointCount();
        vector<double>Q_profile(n, boundaries.volumetric_flow); // задаем по трубе новый расход из временного ряда

        PipeQAdvection advection_model(pipe, Q_profile);

        auto density_wrapper = buffer.buffer.get_buffer_wrapper(
            &Layer::get_density_wrapper);

        auto viscosity_wrapper = buffer.buffer.get_buffer_wrapper(
            &Layer::get_viscosity_wrapper);

        if constexpr (std::is_same<Solver, moc_solver<1>>::value) {
            // Шаг по плотности
            moc_solver<1> solver_rho(advection_model, density_wrapper);
            solver_rho.step_optional_boundaries(dt, boundaries.density, boundaries.density);
            // Шаг по вязкости
            moc_solver<1> solver_nu(advection_model, viscosity_wrapper);
            solver_nu.step_optional_boundaries(dt, boundaries.viscosity, boundaries.viscosity);

        }
        else {
            // Шаг по плотности
            quickest_ultimate_fv_solver solver_rho(advection_model, density_wrapper);
            solver_rho.step(dt, boundaries.density, boundaries.density);
            // Шаг по вязкости
            quickest_ultimate_fv_solver solver_nu(advection_model, viscosity_wrapper);
            solver_nu.step(dt, boundaries.viscosity, boundaries.viscosity);

        }
    }
    void calc_pressure_layer(const isothermal_quasistatic_task_boundaries_t& boundaries) {
        // Получаем новый профиль давлений

        auto& current = buffer.buffer.current();

        vector<double>& p_profile = current.pressure;
        int euler_direction = +1; // Задаем направление для Эйлера

        pipe_model_PQ_cell_parties_t pipeModel(pipe, current.density, current.viscosity, boundaries.volumetric_flow, euler_direction);
        solve_euler<1>(pipeModel, euler_direction, boundaries.pressure_in, &p_profile);
        // Получаем дифференциальный профиль давлений
        std::transform(buffer.pressure_initial.begin(), buffer.pressure_initial.end(), p_profile.begin(),
            current.pressure_delta.begin(),
            [](double initial, double current) {return initial - current;  });

    }
public:
    void step(double dt, const isothermal_quasistatic_task_boundaries_t& boundaries) {
        make_rheology_step(dt, boundaries);
        calc_pressure_layer(boundaries);
    }
    void advance()
    {
        buffer.buffer.advance(+1);
    }

    auto& get_buffer()
    {
        return buffer.buffer;
    }
    void print(const vector<double>& layer, const time_t& dt, const string& path, const string& layer_name)
    {
        string filename = get_courant_research_filename_for_qsm(path, layer_name);

        ofstream  file(filename, ios::app);
        if (file.is_open()) {
            file << std::put_time(std::localtime(&dt), "%c") << ";";
            for (int j = 0; j < layer.size(); j++)
            {
                file << to_string(layer[j]) << ";";
            }
            file << "\n";
            file.close();
        }
    }
    void print_all(const double& dt, const string& path) {
        auto& current = buffer.buffer.current();
        print(current.density, dt, path, "density");
        print(current.viscosity, dt, path, "viscosity");
        print(current.pressure, dt, path, "pressure");
        print(current.pressure_delta, dt, path, "pressure_delta");
    }
};

inline std::string prepare_research_folder_for_qsm_model()
{
    auto test_info = ::testing::UnitTest::GetInstance()->current_test_info();
    std::string research_name = std::string(test_info->test_case_name());
    std::string case_name = std::string(test_info->name());

    std::string path = std::string("../research_out/QSM_models/") +
        research_name + "/" + case_name + "/";
    std::filesystem::create_directories(path);
    for (const auto& entry : filesystem::directory_iterator(path)) {
        filesystem::remove_all(entry.path());
    }
    return path;
}

/// @brief Тесты для солвера
class QuasiStationaryModel : public ::testing::Test {
protected:
    /// @brief Параметры трубы
    pipe_properties_t pipe;
protected:

    /// @brief Подготовка к расчету для семейства тестов
    virtual void SetUp() override {
        // Упрощенное задание трубы - 200км, с шагом разбиения для расчтной сетки 100 м, диаметром 514мм
        simple_pipe_properties simple_pipe;
        simple_pipe.length = 200e3;
        simple_pipe.diameter = 0.514;
        simple_pipe.dx = 100;

        pipe = pipe_properties_t::build_simple_pipe(simple_pipe);
    }

public:
    vector_timeseries_t generate_timeseries(const vector<pair<string, double>>& timeseries_initial_values, 
        timeseries_generator_settings settings = timeseries_generator_settings::default_settings(), 
        double jump_time = 0, double jump_value = 0, string name_parameter = "")
    {
        // Задаем время 04.08.2024  16:42:53
        settings.start_time = 1712583773;
        // Генерируем данные
        synthetic_time_series_generator data_time_series(timeseries_initial_values, settings);
        if (jump_time != 0 && name_parameter != "")
        {
            data_time_series.apply_jump(jump_time, jump_value, name_parameter);
        }
        // Получаем данные
        const auto& data = data_time_series.get_data();
        // Помещаем временные ряды в вектор
        vector_timeseries_t params(data);
        return params;
    }


    /// @brief Считаем квазистационарную модель
    /// @tparam Layer Слой для расчета плотности, вязкости для численного метода 
    /// @tparam Solver Численный метод расчета движения партий
    /// @param path Путь с результатом
    /// @param timeseries_initial_values Исходные условия для генерации временных рядов
    template <typename Layer, typename Solver>
    void calc_quasistationary_model(const string& path, 
        const isothermal_quasistatic_task_boundaries_t& initial_boundaries,
        const vector_timeseries_t& params, double dt = std::numeric_limits<double>::quiet_NaN())
    {
        isothermal_quasistatic_task_t<Layer, Solver> task(pipe);
        task.solve(initial_boundaries);

        task.advance();

        time_t t = params.get_start_date(); // Момент времени начала моделирования
        do
        {
            // Интерополируем значения параметров в заданный момент времени
            vector<double> values_in_time_model = params(t);
            isothermal_quasistatic_task_boundaries_t boundaries(values_in_time_model);

            double time_step = dt;
            if (std::isnan(time_step)) {
                double v = boundaries.volumetric_flow / pipe.wall.getArea(); 
                time_step = task.get_time_step_assuming_max_speed(v);
            }
            t += time_step;

            task.step(time_step, boundaries);
            task.advance();
            task.print_all(t - time_step, path);
        } while (t < params.get_end_date());
    }
};

/// @brief Пример испольования метода Quickest Ultimate с гидравлическим расчетом  
TEST_F(QuasiStationaryModel, QuickWithQuasiStationaryModel)
{
    // Создаём папку с результатами и получаем путь к ней
    string path = prepare_research_folder_for_qsm_model();

    // Исходные данные для начального стационарного расчета
    constexpr double density_initial = 850;
    isothermal_quasistatic_task_boundaries_t initial_boundaries({ 0.2, 6e6, density_initial, 15e-6 });

    // Временные ряды краевых условий для квазистационарного расчета
    // Даем скачок по плотности на +10 кг/м^3
    vector<pair<string, double>> timeseries_initial_values = {
        { "Q", initial_boundaries.volumetric_flow }, // "Q" Расход по всей трубе (опционально), (м^3/с)
        { "p_in", initial_boundaries.pressure_in }, // "p_in" Давление на входе (опционально), (Па)
        { "rho_in", 10 + density_initial }, // "rho_in" Плотность жидкости, (кг/м3)
        { "visc_in", initial_boundaries.viscosity }, // "visc_in" Вязкость жидкости, (м2/сек)
    };
    // Вызываем метод расчета квазистационарной модели с помощью Quickest Ultimate
    vector_timeseries_t time_series = generate_timeseries(timeseries_initial_values);
    double dt = 75;
    calc_quasistationary_model<density_viscosity_cell_layer, quickest_ultimate_fv_solver>(
        path, initial_boundaries, time_series, dt);
}
/// @brief Пример испольования метода Quickest Ultimate с гидравлическим расчетом  
TEST_F(QuasiStationaryModel, IdealQuickWithQuasiStationaryModel)
{
    // Создаём папку с результатами и получаем путь к ней
    string path = prepare_research_folder_for_qsm_model();

    // Исходные данные для начального стационарного расчета
    constexpr double density_initial = 850;
    isothermal_quasistatic_task_boundaries_t initial_boundaries({ 0.2, 6e6, density_initial, 15e-6 });

    // Временные ряды краевых условий для квазистационарного расчета
    // Даем скачок по плотности на +10 кг/м^3
    vector<pair<string, double>> timeseries_initial_values = {
        { "Q", initial_boundaries.volumetric_flow }, // "Q" Расход по всей трубе (опционально), (м^3/с)
        { "p_in", initial_boundaries.pressure_in }, // "p_in" Давление на входе (опционально), (Па)
        { "rho_in", 10 + density_initial }, // "rho_in" Плотность жидкости, (кг/м3)
        { "visc_in", initial_boundaries.viscosity }, // "visc_in" Вязкость жидкости, (м2/сек)
    };

        // Объявляем структуру с исходными данными и настроечными параметрами
    timeseries_generator_settings settings = timeseries_generator_settings::default_settings();
    settings.value_relative_decrement = 0;
    settings.value_relative_increment = 0;
    settings.sample_time_max = 200;
    settings.sample_time_min = 200;
    vector_timeseries_t time_series = generate_timeseries(timeseries_initial_values, settings);
    // Вызываем метод расчета квазистационарной модели с помощью Quickest Ultimate
    calc_quasistationary_model<density_viscosity_cell_layer, quickest_ultimate_fv_solver>(
        path, initial_boundaries, time_series);
}
/// @brief Пример испольования метода характеристик с гидравлическим расчетом  
TEST_F(QuasiStationaryModel, MocWithQuasiStationaryModel)
{
    // Создаём папку с результатами и получаем путь к ней
    string path = prepare_research_folder_for_qsm_model();
    // Исходные данные для начального стационарного расчета
    constexpr double density_initial = 850;
    isothermal_quasistatic_task_boundaries_t initial_boundaries({ 0.2, 6e6, density_initial, 15e-6 });

    // Временные ряды краевых условий для квазистационарного расчета
    // Даем скачок по плотности на +10 кг/м^3
    vector<pair<string, double>> timeseries_initial_values = {
        { "Q", initial_boundaries.volumetric_flow }, // "Q" Расход по всей трубе (опционально), (м^3/с)
        { "p_in", initial_boundaries.pressure_in }, // "p_in" Давление на входе (опционально), (Па)
        { "rho_in", 10 + density_initial }, // "rho_in" Плотность жидкости, (кг/м3)
        { "visc_in", initial_boundaries.viscosity }, // "visc_in" Вязкость жидкости, (м2/сек)
    };
    // Вызываем метод расчета квазистационарной модели с помощью МХ
    vector_timeseries_t time_series = generate_timeseries(timeseries_initial_values);
    double dt = 75;
    calc_quasistationary_model<density_viscosity_layer_moc, moc_solver<1>>(
        path, initial_boundaries, time_series, dt);
}
/// @brief Пример испольования метода характеристик с гидравлическим расчетом  
TEST_F(QuasiStationaryModel, OptionalStepMocWithQuasiStationaryModel)
{
    // Создаём папку с результатами и получаем путь к ней
    string path = prepare_research_folder_for_qsm_model();
    // Исходные данные для начального стационарного расчета
    constexpr double density_initial = 850;
    isothermal_quasistatic_task_boundaries_t initial_boundaries({ 0.2, 6e6, density_initial, 15e-6 });

    // Временные ряды краевых условий для квазистационарного расчета
    // Даем скачок по плотности на +10 кг/м^3
    vector<pair<string, double>> timeseries_initial_values = {
        { "Q", initial_boundaries.volumetric_flow }, // "Q" Расход по всей трубе (опционально), (м^3/с)
        { "p_in", initial_boundaries.pressure_in }, // "p_in" Давление на входе (опционально), (Па)
        { "rho_in", 10 + density_initial }, // "rho_in" Плотность жидкости, (кг/м3)
        { "visc_in", initial_boundaries.viscosity }, // "visc_in" Вязкость жидкости, (м2/сек)
    };
    // Вызываем метод расчета квазистационарной модели с помощью МХ
    vector_timeseries_t time_series = generate_timeseries(timeseries_initial_values);
    calc_quasistationary_model<density_viscosity_layer_moc, moc_solver<1>>(
        path, initial_boundaries, time_series);
}

/// @brief Пример испольования метода характеристик с гидравлическим расчетом  
TEST_F(QuasiStationaryModel, IdealMocWithQuasiStationaryModel)
{
    // Создаём папку с результатами и получаем путь к ней
    string path = prepare_research_folder_for_qsm_model();
    // Исходные данные для начального стационарного расчета
    constexpr double density_initial = 850;
    isothermal_quasistatic_task_boundaries_t initial_boundaries({ 0.2, 6e6, density_initial, 15e-6 });

    // Временные ряды краевых условий для квазистационарного расчета
    // Даем скачок по плотности на +10 кг/м^3
    vector<pair<string, double>> timeseries_initial_values = {
        { "Q", initial_boundaries.volumetric_flow }, // "Q" Расход по всей трубе (опционально), (м^3/с)
        { "p_in", initial_boundaries.pressure_in }, // "p_in" Давление на входе (опционально), (Па)
        { "rho_in", 10 + density_initial }, // "rho_in" Плотность жидкости, (кг/м3)
        { "visc_in", initial_boundaries.viscosity }, // "visc_in" Вязкость жидкости, (м2/сек)
    };
    // Объявляем структуру с исходными данными и настроечными параметрами
    timeseries_generator_settings settings = timeseries_generator_settings::default_settings();
    settings.value_relative_decrement = 0;
    settings.value_relative_increment = 0;
    settings.sample_time_max = 200;
    settings.sample_time_min = 200;
    // Вызываем метод расчета квазистационарной модели с помощью МХ
    vector_timeseries_t time_series = generate_timeseries(timeseries_initial_values, settings);
    calc_quasistationary_model<density_viscosity_layer_moc, moc_solver<1>>(
        path, initial_boundaries, time_series);
}

/// @brief Пример испольования метода характеристик с гидравлическим расчетом  
TEST_F(QuasiStationaryModel, IdealImpulsMocWithQuasiStationaryModel)
{
    // Создаём папку с результатами и получаем путь к ней
    string path = prepare_research_folder_for_qsm_model();
    // Исходные данные для начального стационарного расчета
    constexpr double density_initial = 850;
    isothermal_quasistatic_task_boundaries_t initial_boundaries({ 0.2, 6e6, density_initial, 15e-6 });

    // Временные ряды краевых условий для квазистационарного расчета
    // Даем скачок по плотности на +10 кг/м^3
    vector<pair<string, double>> timeseries_initial_values = {
        { "Q", initial_boundaries.volumetric_flow }, // "Q" Расход по всей трубе (опционально), (м^3/с)
        { "p_in", initial_boundaries.pressure_in }, // "p_in" Давление на входе (опционально), (Па)
        { "rho_in", 10 + density_initial }, // "rho_in" Плотность жидкости, (кг/м3)
        { "visc_in", initial_boundaries.viscosity }, // "visc_in" Вязкость жидкости, (м2/сек)
    };
    // Объявляем структуру с исходными данными и настроечными параметрами
    timeseries_generator_settings settings = timeseries_generator_settings::default_settings();
    settings.value_relative_decrement = 0;
    settings.value_relative_increment = 0;
    settings.sample_time_max = 200;
    settings.sample_time_min = 200;
    double jump_time = 5000;
    double jump_value = -10;
    // Вызываем метод расчета квазистационарной модели с помощью МХ
    vector_timeseries_t time_series = generate_timeseries(timeseries_initial_values, settings, jump_time, jump_value, "rho_in");
    calc_quasistationary_model<density_viscosity_layer_moc, moc_solver<1>>(
        path, initial_boundaries, time_series);
}