#pragma once


/// @brief Тесты для расчёта на реальных данных
class IdentNonisothermalQSM : public ::testing::Test {
protected:
    /// @brief Путь к реальным данным с Линейного участка трубопровода
    const std::string data_path = "../research/2025-01-ident-quasistatic-nonisothermal/data/";
    
    
    // @brief Создание модели трубопровода по реальным данным 
    /// @param path Путь к файлу с профилем реального ЛУ
    /// @return Модель трубы с профилем на основе профиля реального участка трубы 
    static pipe_noniso_properties_t prepare_pipe(const std::string& path)
    {
        // Указываем имя файла
        std::string folder = path + "coord_heights.csv";
        //Желаемый шаг
        double desired_dx = 200;

        pipe_noniso_properties_t pipe;

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
    static std::tuple<vector<double>, vector<vector<double>>, vector<double>> prepare_real_data(const std::string& path_to_real_data)
    {
        // Временные ряды краевых условий
        vector<pair<vector<time_t>, vector<double>>> control_tag_data;
        // Временные ряды эталонных данных
        vector<pair<vector<time_t>, vector<double>>> etalon_tag_data;
        // Задаём период
        string start_period = "01.08.2021 00:00:00";
        string end_period = "01.09.2021 00:00:00";

        // Прописываем названия файлов и единицы измерения параметров
        vector<pair<string, string>>parameters =
        {
            { path_to_real_data + "Q_in", "m3/h-m3/s"s },
            { path_to_real_data + "t_in_n", "C"s },
            { path_to_real_data + "rho_in", "kg/m3"s },
            { path_to_real_data + "visc_in", "mm^2/s-m^2/s"s },
            { path_to_real_data + "t_out_n", "C"s}

        };
        // Считываем временные ряды параметров
        csv_multiple_tag_reader tags(parameters);
        
        control_tag_data = tags.read_csvs(start_period, end_period); //чувствителен к отсутствию явного формата секунд для некоторых csv файлов (по опыту, когда секунды :00)
        etalon_tag_data = { control_tag_data.back() };
        control_tag_data.pop_back();

        // Помещаем временные ряды в вектор
        vector_timeseries_t control_parameters_time_series(control_tag_data);
        vector_timeseries_t etalon_parameters_time_series(etalon_tag_data);

        // Задаём шаг для временной сетки - 1 минута
        double step = 60;

        // Определяем начало периода
        time_t start_period_time = max(control_parameters_time_series.get_start_date(), etalon_parameters_time_series.get_start_date());
        // Определяем конец периода
        time_t end_period_time = min(control_parameters_time_series.get_end_date(), etalon_parameters_time_series.get_end_date());
        // Определяем продолжительность периода
        time_t duration = (end_period_time - start_period_time);

        // Считаем количество точек в сетке
        size_t dots_count = static_cast<size_t>(ceil(duration / step) + 0.00001);

        vector<double>  times = vector<double>(dots_count);
        vector<vector<double>> control_data = vector<vector<double>>(dots_count);
        vector<double> etalon_pressure = vector<double>(dots_count);

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

/// @brief Пример идентификации модели трубопровода по коэффициенту теплообмена
TEST_F(IdentNonisothermalQSM, HTC) 
{
    using namespace pde_solvers;

    // Подготавливаем модель трубопровода и параметры для идентификации 
    //pipe_noniso_properties_t pipe = get_noniso_default_pipe();
    pipe_noniso_properties_t pipe = prepare_pipe(data_path);
    auto [times, control_data, etalon_pressure] = prepare_real_data(data_path);

    // Выбираем в настройках параметр идентификации - коэффициент теплообмена
    ident_nonisothermal_qsm_pipe_settings ident_settings;
    ident_settings.ident_htc = true;

    // Создаём класс для идентификации
    ident_nonisothermal_qsm_pipe_parameters_t test_ident(ident_settings, pipe, times, control_data, etalon_pressure);

    // Создаём сущности для хранения результата и аналитики
    fixed_optimizer_result_t result;
    fixed_optimizer_result_analysis_t analysis;

    double result_d = test_ident.ident(&result, &analysis);
}
