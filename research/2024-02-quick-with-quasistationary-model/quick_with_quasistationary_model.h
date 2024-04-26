#pragma once

#include <pde_solvers/timeseries.h>
#include <iostream>
#include <vector>
#include <string>

using std::vector;
using std::pair;
using std::string;
using std::time_t;

/// @brief Слой для расчета плотности, вязкости методом конечных объемов 
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

/// @brief Проблемно-ориентированный слой
/// Задача про плотность и вязкость
struct density_viscosity_layer_moc
{
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
    /// @brief Подготовка плотности для расчета по методу характеристик
    /// Оборачивает профиль плотности и вспомогательный расчет МХ в обертку для МХ
    static moc_layer_wrapper<1> get_density_wrapper(density_viscosity_layer_moc& layer)
    {
        return moc_layer_wrapper<1>(layer.density, layer.specific);
    }
    /// @brief Подготовка вязкости для расчета методом характеристик moc_solver
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
public:
    /// @brief Констуктор уравнения трубы
    /// @param pipe Ссылка на сущность трубы
    /// @param oil Ссылка на сущность нефти
    /// @param flow Объемный расход
    /// @param solver_direction Направление расчета по Эйлеру, должно обязательно совпадать с параметром солвера Эйлера
    pipe_model_PQ_cell_parties_t(const pipe_properties_t& pipe, const vector<double>& rho_profile, const vector<double>& nu_profile, double flow,
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
struct task_buffer_t {
    /// @brief Изначальный профиль давления
    vector<double> pressure_initial;
    /// @brief Буфер профилей давления, плотности, вязкости
    ring_buffer_t<Layer> buffer;
    task_buffer_t(size_t point_count)
        : pressure_initial(point_count)
        , buffer(2, point_count)
    {}
};

/// @brief Структура, содержащая в себе начальные условия задачи PQ
struct quasistatic_task_boundaries_t {
    /// @brief Изначальный объемный расход
    double volumetric_flow;
    /// @brief Изначальное давление на входе
    double pressure_in;
    /// @brief Изначальная плотность на входе
    double density;
    /// @brief Изначальная вязкость на входе
    double viscosity;

    static quasistatic_task_boundaries_t default_values() {
        quasistatic_task_boundaries_t result;
        result.volumetric_flow = 0.2;
        result.pressure_in = 6e6;
        result.density = 850;
        result.viscosity = 15e-6;
        return result;
    }
};
/// @brief Проблемно-ориентированный класс для задачи PQ
template <typename Layer, typename Solver>
class quasistatic_task_t {
    pipe_properties_t pipe;
    task_buffer_t<Layer> buffer;

public:
    quasistatic_task_t(const pipe_properties_t& pipe,
        const quasistatic_task_boundaries_t& initial_conditions)
        : pipe(pipe)
        , buffer(pipe.profile.getPointCount())
    {
        size_t n = pipe.profile.getPointCount();

        // Инициализация реологии
        auto& current = buffer.buffer.current();
        if constexpr (std::is_same<Solver, moc_solver<1>>::value) {
            current.density = vector<double>(n, initial_conditions.density); // Инициализация начального профиля плотности
            current.viscosity = vector<double>(n, initial_conditions.viscosity); // Инициализация начального профиля вязкости
        }
        else {
            current.density = vector<double>(n - 1, initial_conditions.density); // Инициализация начального профиля плотности
            current.viscosity = vector<double>(n - 1, initial_conditions.viscosity); // Инициализация начального профиля вязкости
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

    void step(double dt, const quasistatic_task_boundaries_t& boundaries) {
        size_t n = pipe.profile.getPointCount();
        auto& current = buffer.buffer.current();

        vector<double>& p_profile = current.pressure;
        vector<double> density_buffer;
        vector<double> viscosity_buffer;
        auto density_wrapper = buffer.buffer.get_buffer_wrapper(
            &Layer::get_density_wrapper);

        auto viscosity_wrapper = buffer.buffer.get_buffer_wrapper(
            &Layer::get_viscosity_wrapper);



        int euler_direction = +1; // Задаем направление для Эйлера


        vector<double>Q_profile(n, boundaries.volumetric_flow); // задаем по трубе новый расход из временного ряда

        PipeQAdvection advection_model(pipe, Q_profile);

        if constexpr (std::is_same<Solver, moc_solver<1>>::value) {
            // Шаг по плотности
            Solver solver_rho(advection_model, density_wrapper);
            solver_rho.step_optional_boundaries(dt, boundaries.density, boundaries.density);
            // Шаг по вязкости
            Solver solver_nu(advection_model, viscosity_wrapper);
            solver_nu.step_optional_boundaries(dt, boundaries.viscosity, boundaries.viscosity);
            density_buffer = density_wrapper.current().values;
            viscosity_buffer = viscosity_wrapper.current().values;
        }
        else {
            // Шаг по плотности
            Solver solver_rho(advection_model, density_wrapper);
            solver_rho.step(dt, boundaries.density, boundaries.density);
            // Шаг по вязкости
            Solver solver_nu(advection_model, viscosity_wrapper);
            solver_nu.step(dt, boundaries.viscosity, boundaries.viscosity);
            density_buffer = density_wrapper.current().vars;
            viscosity_buffer = viscosity_wrapper.current().vars;
        }


        pipe_model_PQ_cell_parties_t pipeModel(pipe, density_buffer, viscosity_buffer, Q_profile[0], euler_direction);
        // Получаем новый профиль давлений
        solve_euler<1>(pipeModel, euler_direction, boundaries.pressure_in, &p_profile);
        // Получаем дифференциальный профиль давлений
        std::transform(buffer.pressure_initial.begin(), buffer.pressure_initial.end(), p_profile.begin(),
            current.pressure_delta.begin(),
            [](double initial, double current) {return initial - current;  });

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
    /// @brief Считаем квазистационарную модель
    /// @tparam Layer Слой для расчета плотности, вязкости для численного метода 
    /// @tparam Solver Численный метод расчета движения партий
    /// @param path Путь с результатом
    /// @param timeseries_initial_values Исходные условия для генерации временных рядов
    template <typename Layer, typename Solver>
    void calc_quasistationary_model(string path, vector<pair<string, double>> timeseries_initial_values)
    {
        // Объявляем структуру с исходными данными и настроечными параметрами
        timeseries_generator_settings settings = timeseries_generator_settings::default_values();
        // Задаем время 04.08.2024  16:42:53
        settings.start_time = 1712583773;
        // Генерируем данные
        synthetic_time_series_generator data_time_series(timeseries_initial_values, settings);

        // Получаем данные
        const auto& data = data_time_series.get_data();
        // Помещаем временные ряды в вектор
        vector_timeseries_t params(data);

        task_buffer_t<Layer> buffer(pipe.profile.getPointCount());

        quasistatic_task_boundaries_t initial_boundaries = { 0.2, 6e6, 850, 15e-6 };
        quasistatic_task_t<Layer, Solver> task(pipe, initial_boundaries);

        task.advance();

        double v_max = 1; // Предполагаем скорость для Куранта = 1, скорость, больше чем во временных рядах и в профиле
        time_t dt = task.get_time_step_assuming_max_speed(v_max);
        time_t t = settings.start_time; // Момент времени начала моделирования
        do
        {
            t += dt;
            // Интерополируем значения параметров в заданный момент времени
            vector<double> values_in_time_model = params(t);

            quasistatic_task_boundaries_t boundaries;
            boundaries.volumetric_flow = values_in_time_model[0];
            boundaries.pressure_in = values_in_time_model[1];
            boundaries.density = values_in_time_model[2];
            boundaries.viscosity = values_in_time_model[3];

            task.step(dt, boundaries);
            task.advance();
            task.print_all(t - dt, path);
        } while (t < settings.start_time + settings.duration - dt);
    }
};

/// @brief Пример испольования метода Quickest Ultimate с гидравлическим расчетом  
TEST_F(QuasiStationaryModel, QuickWithQuasiStationaryModel)
{
    // Создаём папку с результатами и получаем путь к ней
    string path = prepare_research_folder_for_qsm_model();
    // Создаём исходные данные
    vector<pair<string, double>> timeseries_initial_values = {
        { "Q", 0.2 }, // "Q" Расход по всей трубе (опционально), (м^3/с)
        { "p_in", 6e6}, // "p_in" Давление на входе (опционально), (Па)
        { "rho_in", 860 }, // "rho_in" Плотность жидкости, (кг/м3)
        { "visc_in", 15e-6},
    };
    // Вызываем метод расчета квазистационарной модели с помощью Quickest Ultimate
    calc_quasistationary_model<density_viscosity_cell_layer, quickest_ultimate_fv_solver>(path, timeseries_initial_values);
}
/// @brief Пример испольования метода характеристик с гидравлическим расчетом  
TEST_F(QuasiStationaryModel, MocWithQuasiStationaryModel)
{
    // Создаём папку с результатами и получаем путь к ней
    string path = prepare_research_folder_for_qsm_model();
    // Создаём исходные данные
    vector<pair<string, double>> timeseries_initial_values = {
        { "Q", 0.2 }, // "Q" Расход по всей трубе (опционально), (м^3/с)
        { "p_in", 6e6}, // "p_in" Давление на входе (опционально), (Па)
        { "rho_in", 860 }, // "rho_in" Плотность жидкости, (кг/м3)
        { "visc_in", 15e-6},
    };
    // Вызываем метод расчета квазистационарной модели с помощью МХ
    calc_quasistationary_model<density_viscosity_layer_moc, moc_solver<1>>(path, timeseries_initial_values);
}