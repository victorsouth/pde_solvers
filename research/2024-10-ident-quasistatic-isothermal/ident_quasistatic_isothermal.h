#pragma once

#include <Eigen/Sparse>

/// @brief Тесты для расчёта на реальных данных
class IdentIsothermalQSM : public ::testing::Test {
protected:
    /// @brief Путь к реальным данным с Линейного участка трубопровода
    const std::string data_path = "../research/2024-08-quasistationary-with-real-data/data/";

    /// @brief Создание модели трубопровода по реальным данным 
    /// @param path Путь к файлу с профилем реального ЛУ
    /// @return Модель трубы с профилем на основе профиля реального участка трубы
    static pipe_properties_t prepare_pipe(const std::string& path)
    {
        // Указываем имя файла
        std::string folder = path + "coord_heights.csv";
        //Желаемый шаг
        double desired_dx = 200;

        pipe_properties_t pipe;

        // Создаём новый профиль с постоянным шагом
        pipe.profile = pipe_profile_uniform::get_uniform_profile_from_csv(desired_dx, folder);
        // Номинальный диаметр трубы
        pipe.wall.diameter = 1;

        return pipe;
    };

    /// @brief Для моделирования работы ЛУ считываем значения параметров с реального трубопровода
    /// и на основе этих данных предпосчитываем интерполированные значения 
    /// краевых условий (плотность, вязкость, расход и давление в начале ЛУ) 
    /// и эталонного параметра (давление на выходе трубопровода)
    /// для временной сетки с шагом равным одной минуте
    /// @param path_to_real_data Путь к данным с реального трубопровода
    /// @return Кортёж, состоящий из 3-х элементов
    /// 1. Временная сетка
    /// 2. Предпосчитанные краевые условия
    /// 3. Предпосчитанные эталонные значения
    static std::tuple<std::vector<double>, std::vector<std::vector<double>>, std::vector<double>> prepare_real_data(const std::string& path_to_real_data)
    {
        // Временные ряды краевых условий
        std::vector<std::pair<std::vector<time_t>, std::vector<double>>> control_tag_data;
        // Временные ряды эталонных данных
        std::vector<std::pair<std::vector<time_t>, std::vector<double>>> etalon_tag_data;
        // Задаём период
        std::string start_period = "01.08.2021 00:00:00";
        std::string end_period = "01.09.2021 00:00:00"; 
        using namespace std::string_literals;

        // Прописываем названия файлов и единицы измерения параметров
        std::vector<std::pair<std::string, std::string>>parameters =
        {
            { path_to_real_data + "Q_in", "m3/h-m3/s"s },
            { path_to_real_data + "p_in", "MPa"s },
            { path_to_real_data + "rho_in", "kg/m3"s },
            { path_to_real_data + "visc_in", "mm^2/s-m^2/s"s },
            { path_to_real_data + "p_out", "MPa"s}

        };

        // Считываем временные ряды параметров
        csv_multiple_tag_reader tags(parameters);
        control_tag_data = tags.read_csvs(start_period, end_period);
        etalon_tag_data = { control_tag_data.back() };
        control_tag_data.pop_back();

        // Помещаем временные ряды в вектор
        vector_timeseries_t control_parameters_time_series(control_tag_data);
        vector_timeseries_t etalon_parameters_time_series(etalon_tag_data);

        // Задаём шаг для временной сетки - 1 минута
        double step = 60;

        // Определяем начало периода
        time_t start_period_time = std::max(control_parameters_time_series.get_start_date(), etalon_parameters_time_series.get_start_date());
        // Определяем конец периода
        time_t end_period_time = std::min(control_parameters_time_series.get_end_date(), etalon_parameters_time_series.get_end_date());
        // Определяем продолжительность периода
        time_t duration = (end_period_time - start_period_time);

        // Считаем количество точек в сетке
        size_t dots_count = static_cast<size_t>(ceil(duration / step) + 0.00001);

        std::vector<double>  times(dots_count);
        std::vector<std::vector<double>> control_data(dots_count);
        std::vector<double> etalon_pressure(dots_count);

        for (size_t i = 0; i < dots_count; i++)
        {
            // Заполняем временную сетку
            times[i] = step * i;
            // Определяем момент времени
            time_t t = start_period_time + static_cast<time_t>(times[i] + 0.5);

            // Получаем интерполированные значения краевых условий и эталонных значений
            control_data[i] = control_parameters_time_series(t);
            etalon_pressure[i] = etalon_parameters_time_series(t).front();

        };

        return std::make_tuple(std::move(times), std::move(control_data), std::move(etalon_pressure));
    }
};

/// @brief Проверяет идентификацию эффективного диаметра модели трубопровода по реальным данным СДКУ
TEST_F(IdentIsothermalQSM, IdentifiesDiameterFromTelemetry)
{
    // Arrange: подготовка модели трубопровода, краевых условий и эталона СДКУ
    pipe_properties_t pipe = prepare_pipe(data_path);
    auto [times, control_data, etalon_pressure] = prepare_real_data(data_path);

    ident_isothermal_qsm_pipe_settings ident_settings;
    ident_settings.ident_diameter = true;

    ident_isothermal_qsm_pipe_parameters_t test_ident(ident_settings, pipe, times, control_data, etalon_pressure);
    
    fixed_optimizer_result_t result;
    fixed_optimizer_result_analysis_t analysis;

    Eigen::VectorXd initial_arg = Eigen::VectorXd::Ones(1);
    double initial_residuals_norm = test_ident.residuals(initial_arg).norm();

    // Act: запуск идентификации по диаметру
    double result_d = test_ident.ident(&result, &analysis);

    double final_residuals_norm = test_ident.residuals(result.argument).norm();

    // Assert: проверка успешной сходимости, физического диапазона относительного диаметра и уменьшения невязок
    ASSERT_EQ(result.result_code, numerical_result_code_t::Converged);
    ASSERT_FALSE(std::isnan(result_d));
    ASSERT_NEAR(result_d, 0.99996, 0.99998);
    ASSERT_LT(final_residuals_norm, initial_residuals_norm);
}

/// @brief Проверяет идентификацию коэффициента гидравлического сопротивления модели трубопровода по реальным данным СДКУ
TEST_F(IdentIsothermalQSM, IdentifiesFrictionFromTelemetry)
{
    // Arrange: подготовка модели трубопровода, краевых условий и эталона СДКУ
    pipe_properties_t pipe = prepare_pipe(data_path);
    auto [times, control_data, etalon_pressure] = prepare_real_data(data_path);

    ident_isothermal_qsm_pipe_settings ident_settings;
    ident_settings.ident_friction = true;

    ident_isothermal_qsm_pipe_parameters_t test_ident(ident_settings, pipe, times, control_data, etalon_pressure);

    fixed_optimizer_result_t result;
    fixed_optimizer_result_analysis_t analysis;

    Eigen::VectorXd initial_arg = Eigen::VectorXd::Ones(1);
    double initial_residuals_norm = test_ident.residuals(initial_arg).norm();

    // Act: запуск идентификации по гидравлическому сопротивлению
    double result_friction = test_ident.ident(&result, &analysis);

    double final_residuals_norm = test_ident.residuals(result.argument).norm();

    // Assert: проверка успешной сходимости, физического диапазона поправки сопротивления и уменьшения невязок
    ASSERT_EQ(result.result_code, numerical_result_code_t::Converged);
    ASSERT_FALSE(std::isnan(result_friction));
    ASSERT_NEAR(result_friction, 1.00019, 1e-4);
    EXPECT_LT(final_residuals_norm, initial_residuals_norm);
}

/// @brief Функция вывода в csv историчности алгоритма идентификации,
/// а также временного ряда невязок до и после идентификации
/// @param times Временная сетка исследуемого периода
/// @param test_ident Класс идентификации
/// @param analysis Сущность с результатами идентификации:
/// значения аргумента, целевой функции и величину шага изменения 
/// @param path_to_ident_results 
void print_identification_result(
    const std::vector<double>& times,
    ident_isothermal_qsm_pipe_parameters_t& test_ident,
    const fixed_optimizer_result_analysis_t& analysis,
    const std::string& path_to_ident_results
) 
{
    auto& target_function = analysis.target_function;
    auto& argument_history = analysis.argument_history;
    auto& steps = analysis.steps;

    auto& init_arg = argument_history.front();
    auto& res_arg = argument_history.back();

    Eigen::VectorXd residuals_initial = test_ident.residuals(init_arg);
    Eigen::VectorXd residuals_result = test_ident.residuals(res_arg);

    std::vector<double>  initial(residuals_initial.data(), residuals_initial.data() + residuals_initial.size());
    std::vector<double>  results(residuals_result.data(), residuals_result.data() + residuals_result.size());

    python_printer<quickest_ultimate_fv_solver> printer;

    // Вывод невязок до и после идентификации
    printer.print_profiles<double>(static_cast<time_t>(0),
        times,
        std::vector<std::vector<double>>{ initial, results },
        "time,time,diff_press_before_ident,diff_press_after_ident",
        path_to_ident_results + "ident_diff_press.csv");


    // Вывод истории изменения параметров
    for (size_t index = 0; index < target_function.size(); index++)
    {
        double step = (index == target_function.size() - 1) ? 0.0 : steps[index];
        printer.print_profiles<size_t>(static_cast<time_t>(0),
            std::vector<size_t>{ index },
            std::vector<std::vector<double>>{ { target_function[index].front()  }, { argument_history[index][0]}, {step} },
            "time,step,target_function,argument_history,steps",
            path_to_ident_results + "ident_history.csv");
    }
}

/// @brief Проверяет идентификацию диаметра c сохранением невязок и истории оптимизации в файлы
TEST_F(IdentIsothermalQSM, IdentifiesDiameterAndPrintsResults)
{
    // Arrange: подготовка модели трубопровода, краевых условий и эталона СДКУ
    pipe_properties_t pipe = prepare_pipe(data_path);
    auto [times, control_data, etalon_pressure] = prepare_real_data(data_path);

    ident_isothermal_qsm_pipe_settings ident_settings;
    ident_settings.ident_diameter = true;

    ident_isothermal_qsm_pipe_parameters_t test_ident(ident_settings, pipe, times, control_data, etalon_pressure);

    fixed_optimizer_result_t result;
    fixed_optimizer_result_analysis_t analysis;

    Eigen::VectorXd initial_arg = Eigen::VectorXd::Ones(1);
    double initial_residuals_norm = test_ident.residuals(initial_arg).norm();

    // Act: проведение идентификации и вывод отчетов в файл
    double result_d = test_ident.ident(&result, &analysis);

    std::string path_to_ident_results = prepare_research_folder_for_qsm_model2();
    print_identification_result(times, test_ident, analysis, path_to_ident_results);

    double final_residuals_norm = test_ident.residuals(result.argument).norm();

    // Assert: проверка успешной оптимизации и сохранения отчетов
    ASSERT_EQ(result.result_code, numerical_result_code_t::Converged);
    ASSERT_FALSE(std::isnan(result_d));
    EXPECT_LT(final_residuals_norm, initial_residuals_norm);
    EXPECT_TRUE(std::filesystem::exists(path_to_ident_results + "ident_diff_press.csv"));
}

/// @brief Проверяет идентификацию коэффициента трения c сохранением невязок и истории оптимизации в файлы
TEST_F(IdentIsothermalQSM, IdentifiesFrictionAndPrintsResults)
{
    // Arrange: подготовка модели трубопровода, краевых условий и эталона СДКУ
    pipe_properties_t pipe = prepare_pipe(data_path);
    auto [times, control_data, etalon_pressure] = prepare_real_data(data_path);

    ident_isothermal_qsm_pipe_settings ident_settings;
    ident_settings.ident_friction = true;

    ident_isothermal_qsm_pipe_parameters_t test_ident(ident_settings, pipe, times, control_data, etalon_pressure);

    fixed_optimizer_result_t result;
    fixed_optimizer_result_analysis_t analysis;

    Eigen::VectorXd initial_arg = Eigen::VectorXd::Ones(1);
    double initial_residuals_norm = test_ident.residuals(initial_arg).norm();

    // Act: проведение идентификации по трению и вывод отчетов в файл
    double result_friction = test_ident.ident(&result, &analysis);

    std::string path_to_ident_results = prepare_research_folder_for_qsm_model2();
    print_identification_result(times, test_ident, analysis, path_to_ident_results);

    double final_residuals_norm = test_ident.residuals(result.argument).norm();

    // Assert: проверка успешной оптимизации и сохранения отчетов
    ASSERT_EQ(result.result_code, numerical_result_code_t::Converged);
    ASSERT_FALSE(std::isnan(result_friction));
    EXPECT_LT(final_residuals_norm, initial_residuals_norm);
    EXPECT_TRUE(std::filesystem::exists(path_to_ident_results + "ident_diff_press.csv"));
}


