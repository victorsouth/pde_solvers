#pragma once

#include <iostream>
#include <vector>
#include <string>
#include <random>
#include <ctime>
#include <algorithm>
#include <iomanip>

using std::vector;
using std::pair;
using std::string;
using std::time_t;

/// @brief Функция расчета скорости из расхода
/// @param Q Расход, (м^3/с)
/// @param internal_diameter Внутренний диаметр, (м)
/// @return Скорость, (м/с)
double calc_speed(double Q, double internal_diameter)
{
    double speed = (4 * Q) / (pow(internal_diameter, 2) * M_PI);
    return speed;
}

// Типы синонимов для улучшения читаемости кода
using TimeVector = vector<time_t>;
using ParamVector = vector<double>;
using ParamPair = pair<TimeVector, ParamVector>;

/// @brief Исходны данные и настроечные параметры
struct timeseries_generator_settings {
    /// @brief Время моделирования, с
    std::time_t duration;
    /// @brief Минимальное значение размаха шага, с
    std::time_t sample_time_min;
    /// @brief Максимальное значение размаха шага, с
    std::time_t sample_time_max;
    /// @brief Относительное минимальное отклонение значения параметров, доли
    double value_relative_decrement;
    /// @brief Относительное максимальное отклонение значения параметров, доли
    double value_relative_increment;
    /// @brief Исходные данные, обязательно должны присутствовать два опциональных параметра
    /// @brief "rho_in" Плотность жидкости, (кг/м3)
    /// @brief "visc_in" Кинематическая вязкость, (м2/с)
    /// @brief "p_in" Давление на входе (опционально), (Па)
    /// @brief "p_out" Давление на выходе (опционально), (Па)
    /// @brief "Q" Расход по всей трубе (опционально), (м^3/с)
    vector<pair<string, double>> timeseries_initial_values = {
        { "Q", 0.2 },
        { "p_in", 6e6},
        { "rho_in", 860 },
        { "visc_in", 15e-6},
    };
    /// @brief Настроечные параметры по умолчанию 
    static timeseries_generator_settings default_values() {
        timeseries_generator_settings result;
        result.duration = 300000;
        result.sample_time_min = 200;
        result.sample_time_max = 400;
        result.value_relative_decrement = 0.0002;
        result.value_relative_increment = 0.0002;
        return result;
    }

};

/// @brief Класс для генерации синтетических временных рядов
class synthetic_time_series_generator {
public:
    /// @brief Конструктор класса
    /// @param settings Настройки генератора временных рядов
    synthetic_time_series_generator(const timeseries_generator_settings& settings)
        : settings_(settings),
        gen(rd()) {
        std::uniform_real_distribution<double> timeDis(settings_.sample_time_min, settings_.sample_time_max);

        for (const auto& param : settings_.timeseries_initial_values) {
            TimeVector timeValues;
            ParamVector paramValues;

            std::uniform_real_distribution<double> normalDis(param.second * (1 - settings_.value_relative_increment), param.second * (1 + settings_.value_relative_decrement));

            time_t timeStep = timeDis(gen);
            for (time_t time = std::time(nullptr); time <= std::time(nullptr) + settings_.duration; time += timeStep) {
                timeValues.push_back(time);
                timeStep = timeDis(gen);
            }

            std::transform(timeValues.begin(), timeValues.end(), std::back_inserter(paramValues),
                [&](time_t) { return normalDis(gen); });

            data.push_back({ timeValues, paramValues });
        }
    }
    /// @brief Применение скачка к временному ряду
    /// @param jump_time Время, когда происходит скачок
    /// @param jump_value Значение скачка
    /// @param paramName Имя параметра, к которому применяется скачок
    void apply_jump(time_t jump_time, double jump_value, const string& paramName) {
        for (size_t i = 0; i < settings_.timeseries_initial_values.size(); ++i) {
            if (settings_.timeseries_initial_values[i].first == paramName) {
                auto it = std::lower_bound(data[i].first.begin(), data[i].first.end(), time(nullptr) + jump_time);
                size_t position = std::distance(data[i].first.begin(), it);
                std::uniform_real_distribution<double> normalDis(jump_value * (1 - settings_.value_relative_increment), jump_value * (1 + settings_.value_relative_increment));
                for (size_t j = position; j < data[i].first.size(); ++j) {
                    double value = normalDis(gen);
                    data[i].second[j] = value;
                }
            }
        }
    }
    /// @brief Получение сгенерированных данных
    /// @return Вектор временных рядов
    const vector<ParamPair>& get_data() const {
        return data;
    }

private:
    /// @brief Настройки генератора
    timeseries_generator_settings settings_;
    /// @brief Данные временных рядов
    vector<ParamPair> data;
    /// @brief Генератор случайных чисел
    std::random_device rd;
    /// @brief Генератор псевдослучайных чисел
    std::mt19937 gen;
};

/// @brief Тесты для солвера quickest_ultimate_fv_solver
class QuickWithQuasiStationaryModel : public ::testing::Test {
protected:
    /// @brief Параметры трубы
    pipe_properties_t pipe;
    /// @brief Профиль расхода
    vector<double> Q_profile;
    /// @brief Модель адвекции
    std::unique_ptr<PipeQAdvection> advection_model;
protected:

    /// @brief Подготовка к расчету для семейства тестов
    virtual void SetUp() override {
        // Упрощенное задание трубы - 200км, с шагом разбиения для расчтной сетки 100 м, диаметром 514мм
        simple_pipe_properties simple_pipe;
        simple_pipe.length = 200e3;
        simple_pipe.diameter = 0.514;
        simple_pipe.dx = 100;

        pipe = pipe_properties_t::build_simple_pipe(simple_pipe);
        Q_profile = vector<double>(pipe.profile.getPointCount(), 0.2); // задаем по трубе расход 0.2 м3/с
        advection_model = std::make_unique<PipeQAdvection>(pipe, Q_profile);
    }
};

/// @brief Тесты для солвера moc_solver
class MocWithQuasiStationaryModel : public ::testing::Test {
protected:
    /// @brief Параметры трубы
    pipe_properties_t pipe;
    /// @brief Профиль расхода
    vector<double> Q_profile;
    /// @brief Модель адвекции
    std::unique_ptr<PipeQAdvection> advection_model;
protected:

    /// @brief Подготовка к расчету для семейства тестов
    virtual void SetUp() override {
        // Упрощенное задание трубы - 200км, с шагом разбиения для расчтной сетки 100 м, диаметром 514мм
        simple_pipe_properties simple_pipe;
        simple_pipe.length = 200e3;
        simple_pipe.diameter = 0.514;
        simple_pipe.dx = 100;

        pipe = pipe_properties_t::build_simple_pipe(simple_pipe);
        Q_profile = vector<double>(pipe.profile.getPointCount(), 0.2); // задаем по трубе расход 0.2 м3/с
        advection_model = std::make_unique<PipeQAdvection>(pipe, Q_profile);
    }
};

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
/// @param Cr Число куранта
/// @return Путь к файлу в заивисимости от указанного класса солвера
template <typename Solver>
static std::string get_courant_research_filename_for_qsm(const string& path, const string& layer_name)
{
    std::stringstream filename;
    if constexpr (std::is_same<Solver, quickest_ultimate_fv_solver>::value) {
        filename << path << "output ultimate " << layer_name << ".csv";
    }
    else if constexpr (std::is_same<Solver, moc_solver<1>>::value)
    {
        filename << path << "output MOC " << layer_name << ".csv";
    }
    else {
        throw std::runtime_error("Please specify filename for solver type");
    }
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
    template <typename Solver>
    void print(const vector<double>& layer, const time_t& dt, const string& path, const string& layer_name)
    {
        string filename = get_courant_research_filename_for_qsm<Solver>(path, layer_name);
        
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
    template <typename Solver>
    void print_all(const double& dt, const string& path) {
        auto& current = buffer.buffer.current();
        print<Solver>(current.density, dt, path, "density");
        print<Solver>(current.viscosity, dt, path, "viscosity");
        print<Solver>(current.pressure, dt, path, "pressure");
        print<Solver>(current.pressure_delta, dt, path, "pressure_delta");
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

/// @brief Базовый пример использования метода Quickest Ultimate для уравнения адвекции
TEST_F(QuickWithQuasiStationaryModel, UseCase_Advection_Density_Viscosity)
{
    const auto& x = advection_model->get_grid();
    double dx = x[1] - x[0];
    double v = advection_model->getEquationsCoeffs(0, 0);

    ring_buffer_t<density_viscosity_cell_layer> buffer(2, pipe.profile.getPointCount());

    auto& rho_initial = buffer[0].density;
    auto& viscosity_initial = buffer[0].viscosity;
    rho_initial = vector<double>(rho_initial.size(), 850); // инициализация начальной плотности
    viscosity_initial = vector<double>(viscosity_initial.size(), 1e-5); // инициализация начальной плотности

    buffer.advance(+1);

    {
        auto density_wrapper = buffer.get_buffer_wrapper(
            &density_viscosity_cell_layer::get_density_wrapper);
        quickest_ultimate_fv_solver solver(*advection_model, density_wrapper);

        double dt = abs(dx / v);
        double rho_in = 840; // плотность нефти, закачиваемой на входе трубы при положительном расходе
        double rho_out = 860; // плотность нефти, закачиваемой с выхода трубы при отрицательном расходе
        solver.step(dt, rho_in, rho_out);
    }

    {
        auto viscosity_wrapper = buffer.get_buffer_wrapper(
            &density_viscosity_cell_layer::get_viscosity_wrapper);
        quickest_ultimate_fv_solver solver(*advection_model, viscosity_wrapper);

        double dt = abs(dx / v);
        double visc_in = 2e-5; // вязкость нефти, закачиваемой на входе трубы при положительном расходе
        double visc_out = 0.5e-5;; // вязкость нефти, закачиваемой с выхода трубы при отрицательном расходе
        solver.step(dt, visc_in, visc_out);
    }

    auto& curr = buffer[0];
}
/// @brief Базовый пример генерации синтетических временных рядов
TEST(Random, PrepareTimeSeries)
{
    timeseries_generator_settings settings;
    settings.duration = 350000; // Задаю время моделирования (опционально, по умолчанию 300000)
    settings.sample_time_min = 300; // Задаю минимальное значение размаха шага (опционально, по умолчанию 200)
    settings.sample_time_max = 450; // Задаю максимальное значение размаха шага (опционально, по умолчанию 400)
    settings.value_relative_decrement = 0.005; // Задаю относительное минимальное отклонение значения параметров (опционально, по умолчанию 0.0002)
    settings.value_relative_increment = 0.005; // Задаю относительное максимальное отклонение значения параметров (опционально, по умолчанию 0.0002)
    settings.timeseries_initial_values = {
        { "Q", 0.3 },
        { "p_in", 5e6},
        { "rho_in", 850 },
        { "visc_in", 17e-6},


    };
    synthetic_time_series_generator data_generator(settings);

    const time_t jump_time_rho = 100000;
    const double jump_value_rho = 870;
    data_generator.apply_jump(jump_time_rho, jump_value_rho, "rho_in");

    const time_t jump_time_Q = 150000;
    const double jump_value_Q = 0.1;
    data_generator.apply_jump(jump_time_Q, jump_value_Q, "Q");

    const auto& data = data_generator.get_data();

    vector_timeseries_t params(data);

    // Задаём интересующий нас момент времени
    time_t test_time = static_cast<time_t>(std::time(nullptr) + 200000);

    // Интерополируем значения параметров в заданный момент времени
    vector<double> values_in_test_time = params(test_time);
}
/// @brief Пример испольования метода Quickest Ultimate с гидравлическим расчетом  
TEST_F(QuickWithQuasiStationaryModel, WorkingWithTimeSeries)
{
    // Создаём папку с результатами и получаем путь к ней
    string path = prepare_research_folder_for_qsm_model();
    // Объявляем структуру с исходными данными и настроечными параметрами
    timeseries_generator_settings settings = timeseries_generator_settings::default_values();
    // Генерируем данные
    synthetic_time_series_generator data_time_series(settings);

    const time_t jump_time_Q = 100000;
    const double jump_value_Q = 0.1;
    data_time_series.apply_jump(jump_time_Q, jump_value_Q, "Q");

    // Получаем данные
    const auto& data = data_time_series.get_data();
    // Помещаем временные ряды в вектор
    vector_timeseries_t params(data);

    task_buffer_t<density_viscosity_cell_layer> buffer(pipe.profile.getPointCount());

    quasistatic_task_boundaries_t initial_boundaries = { 0.2, 6e6, 850, 15e-6 };
    quasistatic_task_t<density_viscosity_cell_layer, quickest_ultimate_fv_solver> task(pipe, initial_boundaries);

    task.advance();

    double v_max = 1; // Предполагаем скорость для Куранта = 1, скорость, больше чем во временных рядах и в профиле
    time_t dt = task.get_time_step_assuming_max_speed(v_max);
    time_t t = std::time(nullptr); // Момент времени начала моделирования

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
        task.print_all<quickest_ultimate_fv_solver>(t - dt, path);
    } while (t < std::time(nullptr) + settings.duration - dt);
}
/// @brief Пример испольования метода характеристик с гидравлическим расчетом  
TEST_F(MocWithQuasiStationaryModel, WorkingWithTimeSeries)
{
    // Создаём папку с результатами и получаем путь к ней
    string path = prepare_research_folder_for_qsm_model();
    // Объявляем структуру с исходными данными и настроечными параметрами
    timeseries_generator_settings settings = timeseries_generator_settings::default_values();
    // Генерируем данные
    synthetic_time_series_generator data_time_series(settings);

    const time_t jump_time_Q = 100000;
    const double jump_value_Q = 0.1;
    data_time_series.apply_jump(jump_time_Q, jump_value_Q, "Q");

    // Получаем данные
    const auto& data = data_time_series.get_data();
    // Помещаем временные ряды в вектор
    vector_timeseries_t params(data);

    task_buffer_t<density_viscosity_layer_moc> buffer(pipe.profile.getPointCount());

    quasistatic_task_boundaries_t initial_boundaries = { 0.2, 6e6, 850, 15e-6 };
    quasistatic_task_t<density_viscosity_layer_moc, moc_solver<1>> task(pipe, initial_boundaries);

    task.advance();

    double v_max = 1; // Предполагаем скорость для Куранта = 1, скорость, больше чем во временных рядах и в профиле
    time_t dt = task.get_time_step_assuming_max_speed(v_max);
    time_t t = std::time(nullptr); // Момент времени начала моделирования

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
        task.print_all<moc_solver<1>>(t - dt, path);
    } while (t < std::time(nullptr) + settings.duration - dt);
}