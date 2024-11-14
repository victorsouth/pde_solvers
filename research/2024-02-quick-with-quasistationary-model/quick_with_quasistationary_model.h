#pragma once

#include <pde_solvers/timeseries.h>

inline std::string prepare_research_folder_for_qsm_model(std::string dop_path = "")
{
    auto test_info = ::testing::UnitTest::GetInstance()->current_test_info();
    std::string research_name = std::string(test_info->test_case_name());
    std::string case_name = std::string(test_info->name());

    std::string path = std::string("../research_out/QSM_models/") +
        research_name + "/" + case_name + "/" + dop_path + "/";
    std::filesystem::create_directories(path);
    for (const auto& entry : filesystem::directory_iterator(path)) {
        filesystem::remove_all(entry.path());
    }
    return path;
}


template <typename Solver = advection_moc_solver>
struct matlab_printer {
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

    /// @brief Запись в файл 
    /// @param layer Слой
    /// @param dt Временной шаг моделирования
    /// @param path Путь к файлу
    /// @param layer_name Тип профиля
    void print(const std::vector<double>& layer, const std::time_t dt, const std::string& path, const std::string& layer_name)
    {
        std::string filename = get_courant_research_filename_for_qsm(path, layer_name);

        std::ofstream  file(filename, std::ios::app);
        if (file.is_open()) {
            //std::tm tm_buf;
            //localtime_s(&tm_buf, &dt);
            //file << std::put_time(&tm_buf, "%c") << ";";
            file << UnixToString(dt, "%c") << ";";
            for (int j = 0; j < layer.size(); j++)
            {
                file << std::to_string(layer[j]) << ";";
            }
            file << "\n";
            file.close();
        }
    }

    /// @brief Вывод профилей плотности, вязкости, давления и отклонения давления от начального 
    /// @param path Путь к файлам с профилями
    /// @param t Момент моделирвоания
    /// @param pipe МОдель трубы
    /// @param layer Проблемно-ориентированный слой
    /// @param etalon_values Эталонные значения
    void print_all(
        std::string path,
        const time_t& t,
        const pipe_properties_t& pipe,
        const density_viscosity_quasi_layer<std::is_same<Solver, advection_moc_solver>::value ? false : true >& layer,
        const vector<double>& etalon_values = {}
    ) {
        print(layer.density, t, path, "density");
        print(layer.viscosity, t, path, "viscosity");
        print(layer.pressure, t, path, "pressure");
        print(layer.pressure_delta, t, path, "pressure_delta");

    }

};

/// @brief Тесты для солвера
class IsothermalQuasistaticModel : public ::testing::Test {
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
    /// @brief Создание профиля отличающегося от прямой линии
    /// Возрастает в первойи четвёртой четвертях, убывает во второй и в третьей
    vector<double> create_honest_profile() const
    {
        double z_min = 140;
        double z_max = 180;
        double z_start = (z_max + z_min) / 2;
        double dz = (z_max - z_start) / pipe.profile.get_point_count();

        vector<double> heights(pipe.profile.get_point_count(), z_start);

        const vector<double>& coordinates = pipe.profile.coordinates;
        for (size_t index = 1; index < heights.size(); index++)
        {
            size_t proc = 100 * index / heights.size();
            if (proc <= 25)
            {
                heights[index] = heights[index - 1] + dz;
            } 
            else if (proc <= 75) 
            {
                heights[index] = heights[index - 1] - 1.2 * dz;
            }
            else
            {
                heights[index] = heights[index - 1] + 1.8 * dz;
            }
        }
        
        return heights;
    }

    vector_timeseries_t generate_timeseries(const vector<pair<string, double>>& timeseries_initial_values, 
        timeseries_generator_settings settings = timeseries_generator_settings::default_settings(), 
        time_t jump_time = 0, double jump_value = 0, string name_parameter = "")
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
};

/// @brief Пример использования метода Quickest Ultimate с гидравлическим расчетом  
TEST_F(IsothermalQuasistaticModel, QuickWithQuasiStationaryModel)
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
    perform_quasistatic_simulation<quickest_ultimate_fv_solver, matlab_printer<quickest_ultimate_fv_solver>>(
        path, pipe, initial_boundaries, time_series, ModelType::FullQuasi, dt);
}
/// @brief Пример использования метода Quickest Ultimate с гидравлическим расчетом (идеальные настройки)
TEST_F(IsothermalQuasistaticModel, IdealQuickWithQuasiStationaryModel)
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
    perform_quasistatic_simulation<quickest_ultimate_fv_solver, matlab_printer<quickest_ultimate_fv_solver>>(
        path, pipe, initial_boundaries, time_series, ModelType::FullQuasi);
}
/// @brief Пример использования метода характеристик с гидравлическим расчетом  
TEST_F(IsothermalQuasistaticModel, MocWithQuasiStationaryModel)
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
    perform_quasistatic_simulation<advection_moc_solver, matlab_printer<advection_moc_solver>>(
        path, pipe, initial_boundaries, time_series, ModelType::FullQuasi, dt);
}
/// @brief Пример использования метода характеристик (переменный шаг) с гидравлическим расчетом  
TEST_F(IsothermalQuasistaticModel, OptionalStepMocWithQuasiStationaryModel)
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
    perform_quasistatic_simulation<advection_moc_solver, matlab_printer<advection_moc_solver>>(
        path, pipe, initial_boundaries, time_series, ModelType::FullQuasi);
}

/// @brief Пример использования метода характеристик с гидравлическим расчетом (идеальные настройки)  
TEST_F(IsothermalQuasistaticModel, IdealMocWithQuasiStationaryModel)
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
    perform_quasistatic_simulation<advection_moc_solver, matlab_printer<advection_moc_solver>>(
        path, pipe, initial_boundaries, time_series, ModelType::FullQuasi);
}

/// @brief Пример использования метода характеристик с гидравлическим расчетом (идеальные настройки)
/// Рассматривается пример с импульсной партией нефти
TEST_F(IsothermalQuasistaticModel, IdealImpulsMocWithQuasiStationaryModel)
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
    time_t jump_time = 5000;
    double jump_value = -10;
    // Вызываем метод расчета квазистационарной модели с помощью МХ
    vector_timeseries_t time_series = generate_timeseries(timeseries_initial_values, settings, jump_time, jump_value, "rho_in");
    perform_quasistatic_simulation<advection_moc_solver, matlab_printer<advection_moc_solver>>(
        path, pipe, initial_boundaries, time_series, ModelType::FullQuasi);
}


/// @brief Наглядное влияние выбора профиля на квазистационарный гидравлический расчёт  
TEST_F(IsothermalQuasistaticModel, ShowProfileImpactInQuasiStationaryModel)
{
    // Создаём папку с результатами для расчёта c профилем по первой и последней точкам трубопровода
    string path_start_end_profile = prepare_research_folder_for_qsm_model("start_end_profile");
    // Создаём папку с результатами для расчёта c полным профилем
    string path_full_profile = prepare_research_folder_for_qsm_model("path_full_profile");
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

    vector<double> profile = create_honest_profile();

    pipe.profile.heights = profile;

    // Вызываем метод расчета квазистационарной модели с помощью МХ для полного профиля 
    vector_timeseries_t time_series = generate_timeseries(timeseries_initial_values, settings);
    perform_quasistatic_simulation<advection_moc_solver, matlab_printer<advection_moc_solver>>(
        path_full_profile, pipe, initial_boundaries, time_series, ModelType::FullQuasi);

    pipe.profile = pipe_profile_t::create(
        pipe.profile.get_point_count(), 
        pipe.profile.coordinates.front(),
        pipe.profile.coordinates.back(),
        profile.front(),
        profile.back(),
        10e6
    );
    // Вызываем метод расчета квазистационарной модели с помощью МХ для профиля по первой и последней точкам
    vector_timeseries_t time_series_2 = generate_timeseries(timeseries_initial_values, settings);
    perform_quasistatic_simulation<advection_moc_solver, matlab_printer<advection_moc_solver>>(
        path_start_end_profile, pipe, initial_boundaries, time_series_2, ModelType::FullQuasi);
}



