#pragma once
#include "pde_solvers/tasks/nonisothermal_quasistatic_task_p.h"

/// @brief Класс, содержащий функции для вывода
/// профилей в файл в формате для плоттеров на Python
template <typename Solver = advection_moc_solver>
struct python_temperature_printer {
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

    /// @brief Вывод профилей температур
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
        const vector<double>& etalon_values = {},
        const nonisothermal_quasistatic_PQ_task_boundaries_t& boundaries = nonisothermal_quasistatic_PQ_task_boundaries_t::default_values()
    ) {
        print_t_profiles<double>(
            t,
            pipe.profile.coordinates,
            { layer.density, layer.viscosity, layer.temperature, layer.temp_shukhov_bias},
            "time,coordinates,density,viscosity,temperature,temp_shukhov_bias",
            folder + "results.csv"
        );

        // Вывод временного ряда отклонения расчётного температуры от фактического
        if (!etalon_values.empty())
        {
            double temp_delta = etalon_values[0] - layer.temperature.back();
            double temp_delta_shukhov_bias = etalon_values[0] - layer.temp_shukhov_bias.back();
            double etalon = etalon_values[0];
            double calculated = layer.temperature.back();
            double temp_in = boundaries.temperature;
            double volum = boundaries.volumetric_flow;
            double calculated_shukhov_bias = layer.temp_shukhov_bias.back();
            print_t_profiles<std::string>(static_cast<time_t>(0),
                vector<string>{ UnixToString(t) },
                vector<vector<double>>{ { etalon }, { calculated }, { temp_delta }, { calculated_shukhov_bias }, { temp_delta_shukhov_bias }, { temp_in }, { volum }},
                "time,time,etalon,calculated,diff_temp,calculated_shukhov_bias,temp_delta_shukhov_bias,temp_in,volum",
                folder + "diff_temp.csv");
        }
    }
};

/// @brief Класс, содержащий функции для вывода
/// профилей в файл в формате для плоттеров на Python
template <typename Solver = advection_moc_solver>
struct python_pressure_printer {
    /// @brief Вывод профилей в файл
    /// @tparam OX Тип данных по оси ОX
    /// @param time_moment Время моделирования
    /// @param x Ось ОХ
    /// @param params Вектор параметров для вывода
    /// @param params_name Названия параметров для вывода
    /// @param filename Название файла
    template <typename OX>
    static void print_p_profiles(const time_t& time_moment, const vector<OX>& x, const vector<vector<double>>& params, const std::string& params_name, const std::string& filename) {
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

    /// @brief Вывод профилей давления
    /// @param t Время моделирования
    /// @param pipe Модель трубы
    /// @param layer Слой, хранящий информацию о 
    /// профилях плотности, вязкости и давления
    /// @param folder Название папки с результатом
    static void print_p_all(
        std::string folder,
        const time_t& t,
        const pipe_noniso_properties_t& pipe,
        const density_viscosity_temp_p_quasi_layer<std::is_same<Solver, advection_moc_solver>::value ? false : true >& layer,
        const vector<double>& etalon_values = {}
    ) {
        print_p_profiles<double>(
            t,
            pipe.profile.coordinates,
            { layer.density, layer.viscosity, layer.temperature, layer.temp_shukhov, layer.pressure },
            "time,coordinates,density,viscosity,temperature,temp_shukhov,pressure",
            folder + "results.csv"
        );

        // Вывод временного ряда отклонения расчётного температуры от фактического
        if (!etalon_values.empty())
        {
            double etalon = etalon_values[0];
            double calculated = layer.temperature.back();
            double calculated_shukhov = layer.temp_shukhov.back();
            double calculated_p = layer.pressure.back();
            double pressure_delta = etalon_values[0] - layer.pressure.back();
            print_p_profiles<std::string>(static_cast<time_t>(0),
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
    // Параметры нефти
    oil_parameters_t oil;
    // Путь к результатам ресёрча
    string path;
    // Временные ряды краевых условий
    vector<pair<vector<time_t>, vector<double>>> tag_data;
    /// @brief Временные ряды эталонных данных
    vector<pair<vector<time_t>, vector<double>>> etalon_tag_data;
    /// @brief Путь к реальным данным с трубопровода
    std::string folder = "../research_out/data/";

    /// @brief Подготовка к расчету для семейства тестов
    virtual void SetUp() override {
        using namespace std;
        // Указываем имя файла и желаемый шаг новой сетки
        string file_name = folder + "coord_heights.csv";
        //Желаемый шаг
        double desired_dx = 200;

        // Создаём новый профиль с постоянным шагом
        pipe.profile = pipe_profile_uniform::get_uniform_profile_from_csv(desired_dx, file_name);
        pipe.wall.diameter = 1;
        //Температура грунта по таблице
        pipe.heat.ambientTemperature = 286.15;
        // Настройки взякости
        oil = get_noniso_default_oil();
        oil.density.nominal_density = 860;
        // Посчитано по Филонову для температур 0, 20, 50
        std::array<double, 3> visc{ 27.028120732908523e-6, 10.937124971104140e-6, 2.815361458795777e-6 };
        oil.viscosity = oil_viscosity_parameters_t(visc);

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
    std::string folder = "../research_out/data/";

    /// @brief Подготовка к расчету для семейства тестов
    virtual void SetUp() override {
        using namespace std;

        // Указываем имя файла и желаемый шаг новой сетки
        string file_name = folder + "coord_heights.csv";
        //Желаемый шаг
        double desired_dx = 200;

        // Создаём новый профиль с постоянным шагом
        pipe.profile = pipe_profile_uniform::get_uniform_profile_from_csv(desired_dx, file_name);
        pipe.wall.diameter = 1;
        //Температура грунта по таблице
        pipe.heat.ambientTemperature = 286.15;
        // Настройки взякости
        oil = get_noniso_default_oil();
        oil.density.nominal_density = 860;
        // Посчитано по Филонову для температур 0, 20, 50
        //std::array<double, 3> visc{ 35.2166964842424e-6, 15.1959389818927e-6, 4.30720885400170e-6 };
        std::array<double, 3> visc{ 27.028120732908523e-6, 10.937124971104140e-6, 2.815361458795777e-6 };
        oil.viscosity = oil_viscosity_parameters_t(visc);

        // Создаём папку с результатами и получаем путь к ней
        path = prepare_research_folder_for_qsm_model();

        vector<pair<string, string>>parameters =
        {
            { folder + "Q_in", "m3/h-m3/s"s },
            { folder + "t_in_nps", "C"s },
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

    perform_noniso_quasistatic_simulation<quickest_ultimate_fv_solver, python_temperature_printer<quickest_ultimate_fv_solver>>(
        path, pipe, oil, params, noniso_qsm_model_type::FullQuasi, etalon_params
    );
}

/// @brief Пример использования метода Quickest Ultimate с гидравлическим расчетом  
TEST_F(NonisothermalQuasistaticModelWithRealData_Pressure, QuasiStationaryFullReology_Pressure)
{
    // Помещаем временные ряды в вектор
    vector_timeseries_t params(tag_data);

    // Помещаем временные ряды в вектор
    vector_timeseries_t etalon_params(etalon_tag_data);

    perform_noniso_quasistatic_simulation_p<quickest_ultimate_fv_solver, python_pressure_printer<quickest_ultimate_fv_solver>>(
        path, pipe, oil, params, noniso_qsm_model_type::FullQuasi, etalon_params
    );
}
