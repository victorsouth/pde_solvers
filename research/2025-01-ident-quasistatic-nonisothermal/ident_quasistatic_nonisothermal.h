#pragma once


/// @brief Тесты для расчёта на реальных данных
class IdentNonisothermalQSM : public ::testing::Test {
protected:
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
    static std::tuple<std::vector<double>, std::vector<std::vector<double>>, std::vector<double>> 
        prepare_real_data(const std::string& path_to_real_data)
    {
        std::string start_period = "01.08.2021 00:00:00";
        std::string end_period = "01.09.2021 00:00:00";
        // Задаём шаг для временной сетки - 1 минута
        double dt = 60;

        using namespace std;
        // Прописываем названия файлов и единицы измерения параметров
        std::vector<std::pair<std::string, std::string>>parameters =
        {
            { path_to_real_data + "Q_in", "m3/h-m3/s"s },
            { path_to_real_data + "t_in_nps", "C"s },
            { path_to_real_data + "t_out_n", "C"s}
        };

        return prepare_timeseries_data(dt, start_period, end_period, parameters);
    }
};



/// @brief Пример идентификации модели трубопровода по коэффициенту теплообмена
TEST_F(IdentNonisothermalQSM, HTC) 
{
    /// @brief Путь к реальным данным с Линейного участка трубопровода
    const std::string data_path = "../research_out/data/";

    // Подготавливаем модель трубопровода и параметры для идентификации 
    pipe_noniso_properties_t pipe = get_research_pipe_heatmodel(data_path);
    auto [times, control_data, etalon_temperature] = prepare_real_data(data_path);
    oil_parameters_t oil = get_noniso_research_oil();

    // Выбираем в настройках параметр идентификации - коэффициент теплообмена
    ident_nonisothermal_qsm_pipe_settings ident_settings;
    ident_settings.ident_htc = true;

    // Создаём класс для идентификации
    ident_nonisothermal_qsm_pipe_parameters_t test_ident(
        ident_settings, pipe, oil, 
        times, control_data, etalon_temperature, noniso_qsm_model_type::Dynamic);

    // Создаём сущности для хранения результата и аналитики
    fixed_optimizer_result_t result;
    fixed_optimizer_result_analysis_t analysis;

    double result_d = test_ident.ident(&result, &analysis);
}
