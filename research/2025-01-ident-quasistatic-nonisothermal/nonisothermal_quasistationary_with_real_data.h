#pragma once
#include "pde_solvers/tasks/nonisothermal_quasistatic_task_p.h"

/// @brief Класс, содержащий функции для вывода
/// профилей в файл в формате для плоттеров на Python
template <typename Solver = advection_moc_solver>
struct python_t_printer {
    /// @brief Вывод профилей в файл
    /// @tparam OX Тип данных по оси ОX
    /// @param time_moment Время моделирования
    /// @param x Ось ОХ
    /// @param params Вектор параметров для вывода
    /// @param params_name Названия параметров для вывода
    /// @param filename Название файла
    template <typename OX>
    static void print_t_profiles(const time_t& time_moment, const vector<OX>& x, const vector<vector<double>>& params, const std::string& params_name, const std::string& filename) {
        std::ofstream output_file;
        size_t profCount = params.size();
        std::string time_str = UnixToString(time_moment);
        if (!std::ifstream(filename))
        {
            output_file.open(filename);
            output_file << params_name << std::endl;
        }
        else
            output_file.open(filename, std::ios::app);

        for (int i = 0; i < params[0].size(); i++)
        {
            output_file << time_str << "," << x[i];
            for (int index = 0; index < profCount; index++)
            {
                output_file << ',' << params[index][i];
                if (index == (profCount - 1))
                    output_file << std::endl;
            }
        }
        output_file.close();
    };

    /// @brief Вывод профилей плотности, вязкости и темп
    /// @param t Время моделирования
    /// @param pipe Модель трубы
    /// @param layer Слой, хранящий информацию о 
    /// профилях плотности, вязкости и давления
    /// @param folder Название папки с результатом
    static void print_t_all(
        std::string folder,
        const time_t& t,
        const pipe_noniso_properties_t& pipe,
        const density_viscosity_temp_quasi_layer<std::is_same<Solver, advection_moc_solver>::value ? false : true >& layer,
        const vector<double>& etalon_values = {}
    ) {
        print_t_profiles<double>(
            t,
            pipe.profile.coordinates,
            { layer.density, layer.viscosity, layer.temp, layer.temp_shukhov},
            "time,coordinates,density,viscosity,temp,temp_shukhov",
            folder + "results.csv"
        );

        // Вывод временного ряда отклонения расчётного температуры от фактического
        if (!etalon_values.empty())
        {
            double temp_delta = etalon_values[0] - layer.temp.back();
            double temp_delta_shukhov = etalon_values[0] - layer.temp_shukhov.back();
            double etalon = etalon_values[0];
            double calculated = layer.temp.back();
            double calculated_shukhov = layer.temp_shukhov.back();
            print_t_profiles<std::string>(static_cast<time_t>(0),
                vector<string>{ UnixToString(t) },
                vector<vector<double>>{ { etalon }, { calculated }, { temp_delta }, { calculated_shukhov }, { temp_delta_shukhov } },
                "time,time,etalon,calculated,diff_temp,calculated_shukhov,temp_delta_shukhov",
                folder + "diff_temp.csv");
        }
    }
};

/// @brief Класс, содержащий функции для вывода
/// профилей в файл в формате для плоттеров на Python
template <typename Solver = advection_moc_solver>
struct python_t_p_printer {
    /// @brief Вывод профилей в файл
    /// @tparam OX Тип данных по оси ОX
    /// @param time_moment Время моделирования
    /// @param x Ось ОХ
    /// @param params Вектор параметров для вывода
    /// @param params_name Названия параметров для вывода
    /// @param filename Название файла
    template <typename OX>
    static void print_t_p_profiles(const time_t& time_moment, const vector<OX>& x, const vector<vector<double>>& params, const std::string& params_name, const std::string& filename) {
        std::ofstream output_file;
        size_t profCount = params.size();
        std::string time_str = UnixToString(time_moment);
        if (!std::ifstream(filename))
        {
            output_file.open(filename);
            output_file << params_name << std::endl;
        }
        else
            output_file.open(filename, std::ios::app);

        for (int i = 0; i < params[0].size(); i++)
        {
            output_file << time_str << "," << x[i];
            for (int index = 0; index < profCount; index++)
            {
                output_file << ',' << params[index][i];
                if (index == (profCount - 1))
                    output_file << std::endl;
            }
        }
        output_file.close();
    };

    /// @brief Вывод профилей плотности, вязкости и темп
    /// @param t Время моделирования
    /// @param pipe Модель трубы
    /// @param layer Слой, хранящий информацию о 
    /// профилях плотности, вязкости и давления
    /// @param folder Название папки с результатом
    static void print_t_p_all(
        std::string folder,
        const time_t& t,
        const pipe_noniso_properties_t& pipe,
        const density_viscosity_temp_p_quasi_layer<std::is_same<Solver, advection_moc_solver>::value ? false : true >& layer,
        const vector<double>& etalon_values = {}
    ) {
        print_t_p_profiles<double>(
            t,
            pipe.profile.coordinates,
            { layer.density, layer.viscosity, layer.temp, layer.temp_shukhov, layer.pressure },
            "time,coordinates,density,viscosity,temp,temp_shukhov,pressure",
            folder + "results.csv"
        );

        // Вывод временного ряда отклонения расчётного температуры от фактического
        if (!etalon_values.empty())
        {
            double etalon = etalon_values[0];
            double calculated = layer.temp.back();
            double calculated_shukhov = layer.temp_shukhov.back();
            double calculated_p = layer.pressure.back();
            double pressure_delta = etalon_values[0] - layer.pressure.back();
            print_t_p_profiles<std::string>(static_cast<time_t>(0),
                vector<string>{ UnixToString(t) },
                vector<vector<double>>{ { calculated }, { calculated_shukhov }, { etalon }, { calculated_p }, { pressure_delta } },
                "time,time,calculated,calculated_shukhov,etalon,calculated_p,diff_press",
                folder + "diff_temp.csv");
        }
    }
};

/// @brief Тесты для расчёта на реальных данных
class NonisothermalQuasistaticModelWithRealData : public ::testing::Test {
protected:
    // Параметры трубы
    pipe_noniso_properties_t pipe;
    // Путь к результатам ресёрча
    string path;
    // Временные ряды краевых условий
    vector<pair<vector<time_t>, vector<double>>> tag_data;
    // Временные ряды эталонных данных
    vector<pair<vector<time_t>, vector<double>>> etalon_tag_data;
    // Путь к реальным данным с трубопровода
    std::string folder = "../research/2025-01-ident-quasistatic-nonisothermal/data/";

    /// @brief Подготовка к расчету для семейства тестов
    virtual void SetUp() override {
        // Указываем имя файла и желаемый шаг новой сетки
        string file_name = folder + "coord_heights.csv";
        //Желаемый шаг
        double desired_dx = 200;

        // Создаём новый профиль с постоянным шагом
        pipe.profile = pipe_profile_uniform::get_uniform_profile_from_csv(desired_dx, file_name);
        pipe.wall.diameter = 1;

        // Создаём папку с результатами и получаем путь к ней
        path = prepare_research_folder_for_qsm_model();

        vector<pair<string, string>>parameters =
        {
            { folder + "Q_in", "m3/h-m3/s"s },
            { folder + "t_in_nps", "C"s },
            { folder + "rho_in", "kg/m3"s },
            { folder + "visc_in", "mm^2/s-m^2/s"s },
            { folder + "t_out_n", "C"s}

        };

        // Задаём период
        string start_period = "01.08.2021 00:00:00";
        string end_period = "01.09.2021 00:00:00";

        // Считываем временные ряды параметров
        csv_multiple_tag_reader tags(parameters);
        tag_data = tags.read_csvs(start_period, end_period);
        etalon_tag_data = { tag_data.back() };
        tag_data.pop_back();

        // Считываем временные ряды параметров
        /*csv_multiple_tag_reader etalon_tags(etalon_parameters);
        etalon_tag_data = etalon_tags.read_csvs(start_period, end_period);*/

        
    }
};

/// @brief Тесты для расчёта на реальных данных
class NonisothermalQuasistaticModelWithRealData_Pressure : public ::testing::Test {
protected:
    // Параметры трубы
    pipe_noniso_properties_t pipe;
    // Параметры нефти
    oil_parameters_t oil;
    // Путь к результатам ресёрча
    string path;
    // Временные ряды краевых условий
    vector<pair<vector<time_t>, vector<double>>> tag_data;
    // Временные ряды эталонных данных давления
    vector<pair<vector<time_t>, vector<double>>> etalon_tag_data;
    // Путь к реальным данным с трубопровода
    std::string folder = "../research/2025-01-ident-quasistatic-nonisothermal/data/";

    /// @brief Подготовка к расчету для семейства тестов
    virtual void SetUp() override {
        // Указываем имя файла и желаемый шаг новой сетки
        string file_name = folder + "coord_heights.csv";
        //Желаемый шаг
        double desired_dx = 200;

        // Создаём новый профиль с постоянным шагом
        pipe.profile = pipe_profile_uniform::get_uniform_profile_from_csv(desired_dx, file_name);
        pipe.wall.diameter = 1;

        // Создаём папку с результатами и получаем путь к ней
        path = prepare_research_folder_for_qsm_model();

        vector<pair<string, string>>parameters =
        {
            { folder + "Q_in", "m3/h-m3/s"s },
            { folder + "t_in_n", "C"s },
            { folder + "rho_in", "kg/m3"s },
            { folder + "visc_in", "mm^2/s-m^2/s"s },
            { folder + "p_in", "MPa"s },
            { folder + "p_out", "MPa"s}

        };

        // Задаём период
        string start_period = "01.08.2021 00:00:00";
        string end_period = "01.09.2021 00:00:00";

        // Считываем временные ряды параметров
        csv_multiple_tag_reader tags(parameters);
        tag_data = tags.read_csvs(start_period, end_period);
        etalon_tag_data = { tag_data.back() };
        tag_data.pop_back();

    }
};

/// @brief Пример использования метода Quickest Ultimate с гидравлическим расчетом  
TEST_F(NonisothermalQuasistaticModelWithRealData, QuasiStationaryFullReology)
{
    // Помещаем временные ряды в вектор
    vector_timeseries_t params(tag_data);

    // Помещаем временные ряды в вектор
    vector_timeseries_t etalon_params(etalon_tag_data);

    perform_noniso_quasistatic_simulation<quickest_ultimate_fv_solver, python_t_printer<quickest_ultimate_fv_solver>>(
        path, pipe, params, noniso_qsm_model_type::FullQuasi, etalon_params
    );
}

/// @brief Пример использования метода Quickest Ultimate с гидравлическим расчетом  
TEST_F(NonisothermalQuasistaticModelWithRealData_Pressure, QuasiStationaryFullReology_Pressure)
{
    // Помещаем временные ряды в вектор
    vector_timeseries_t params(tag_data);

    // Помещаем временные ряды в вектор
    vector_timeseries_t etalon_params(etalon_tag_data);

    perform_noniso_quasistatic_simulation_p<quickest_ultimate_fv_solver, python_t_p_printer<quickest_ultimate_fv_solver>>(
        path, pipe, oil, params, noniso_qsm_model_type::FullQuasi, etalon_params
    );
}


/*
/// @brief Стационарный расчёт с текущей реологией 
TEST_F(IsothermalQuasistaticModelWithRealData, StationaryCurrentReology)
{
    // Помещаем временные ряды в вектор
    vector_timeseries_t params(tag_data);

    // Помещаем временные ряды в вектор
    vector_timeseries_t etalon_params(etalon_tag_data);

    // Производим расчёт и записываем результаты в файлы
    perform_quasistatic_simulation<advection_moc_solver, python_printer<advection_moc_solver>>(
        path, pipe, params, QuasistaticModelType::Stationary, etalon_params
    );
}

/// @brief Стационарный расчёт с реологией из начала периода
TEST_F(IsothermalQuasistaticModelWithRealData, StationaryInitialReology)
{
    auto tag_data_initial = tag_data;

    auto& density_time_series = tag_data_initial[2];
    auto& density_values = density_time_series.second; // second - это значения (first - время, здесь не надо)

    // ставим везде начальное значение плотности
    double initial_density = density_values.front(); 
    for (double& value : density_values) {
        value = initial_density;
    }

    auto& viscosity_time_series = tag_data_initial[3];
    auto& viscosity_values = viscosity_time_series.second; // second - это значения (first - время, здесь не надо)

    // ставим везде начальное значение плотности
    double initial_viscosity = viscosity_values.front();
    for (double& value : viscosity_values) {
        value = initial_viscosity;
    }

    // Помещаем временные ряды в вектор
    vector_timeseries_t params(tag_data_initial);

    // Помещаем временные ряды в вектор
    vector_timeseries_t etalon_params(etalon_tag_data);

    perform_quasistatic_simulation<advection_moc_solver, python_printer<advection_moc_solver>>(
        path, pipe, params, QuasistaticModelType::Stationary, etalon_params
    );
}

/// @brief Стационарный расчёт со средней за период реологией 
TEST_F(IsothermalQuasistaticModelWithRealData, StationaryMeanReology)
{
    auto tag_data_initial = tag_data;

    auto& density_values = tag_data_initial[2].second; // second - это значения (first - время, здесь не надо)
    double sum = std::accumulate(density_values.begin(), density_values.end(), 0.0);
    double mean_density = sum / density_values.size();

    for (double& value : density_values) {
        value = mean_density;
    }

    auto& viscosity_values = tag_data_initial[3].second;
    sum = std::accumulate(viscosity_values.begin(), viscosity_values.end(), 0.0);
    double mean_viscosity = sum / viscosity_values.size();

    for (double& value : viscosity_values) {
        value = mean_viscosity;
    }

    // Помещаем временные ряды в вектор
    vector_timeseries_t params(tag_data_initial);

    // Помещаем временные ряды в вектор
    vector_timeseries_t etalon_params(etalon_tag_data);

    perform_quasistatic_simulation<advection_moc_solver, python_printer<advection_moc_solver>>(
        path, pipe, params, QuasistaticModelType::Stationary, etalon_params
    );
}

/// @brief Пример использования метода Quickest Ultimate с гидравлическим расчетом  
TEST_F(IsothermalQuasistaticModelWithRealData, QuasiStationaryFullReology)
{
    // Помещаем временные ряды в вектор
    vector_timeseries_t params(tag_data);

    // Помещаем временные ряды в вектор
    vector_timeseries_t etalon_params(etalon_tag_data);

    perform_quasistatic_simulation<advection_moc_solver, python_printer<advection_moc_solver>>(
        path, pipe, params, QuasistaticModelType::FullQuasi, etalon_params
    );
}

/// @brief Пример использования метода Quickest Ultimate с гидравлическим расчетом  
TEST_F(IsothermalQuasistaticModelWithRealData, QuasiStationaryDensityOnly)
{
    // Помещаем временные ряды в вектор
    vector_timeseries_t params(tag_data);

    // Помещаем временные ряды в вектор
    vector_timeseries_t etalon_params(etalon_tag_data);

    perform_quasistatic_simulation<advection_moc_solver, python_printer<advection_moc_solver>>(
        path, pipe, params, QuasistaticModelType::DensityQuasi, etalon_params
    );
}

/// @brief Пример использования метода Quickest Ultimate с гидравлическим расчетом  
TEST_F(IsothermalQuasistaticModelWithRealData, QuasiStationaryViscosityOnly)
{
    // Помещаем временные ряды в вектор
    vector_timeseries_t params(tag_data);

    // Помещаем временные ряды в вектор
    vector_timeseries_t etalon_params(etalon_tag_data);

    perform_quasistatic_simulation<advection_moc_solver, python_printer<advection_moc_solver>>(
        path, pipe, params, QuasistaticModelType::ViscosityQuasi, etalon_params
    );
}
*/

