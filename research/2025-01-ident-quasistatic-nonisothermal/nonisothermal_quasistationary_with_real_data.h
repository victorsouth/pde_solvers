#pragma once

/// @brief Класс, содержащий функции для вывода
/// профилей в файл в формате для плоттеров на Python
class python_temperature_printer 
    : public batch_processor_precalculated_times<qsm_noniso_T_layer>
{
private:
    const std::vector<double> times;
    /// @brief Управляющие параметры - расход и входная температура
    std::vector<std::vector<double>> control_data;
    /// @brief Эталонные (фактические) значения выхода
    std::vector<double> etalon_temperature_values;
    const pipe_noniso_properties_t pipe;

private:
    std::ofstream boundary_temperatures_file;
public:
    python_temperature_printer(
        const std::string& boundary_temperatures_file,
        const std::vector<double>& times,
        const std::vector<std::vector<double>>& control_data,
        const std::vector<double>& etalon_values,
        const pipe_noniso_properties_t& pipe)
        : times(times)
        , control_data(control_data)
        , etalon_temperature_values(etalon_values)
        , pipe(pipe)
        , boundary_temperatures_file(boundary_temperatures_file)
    {

    }
    /// @brief На каждом шаге моделирования пишет конечные температуры
    virtual void process_data(size_t step_index, const qsm_noniso_T_layer& layer) override
    {
        //print_t_profiles<double>(
        //    output_file,
        //    t,
        //    pipe.profile.coordinates,
        //    { layer.temperature, layer.temperature_shukhov},
        //    "time,coordinates,temperature,temperature_shukhov"
        //);
        
        time_t t = (time_t)times[step_index];
        double etalon_temperature = etalon_temperature_values[step_index];

        // Входные данные модели, которые выведутся в файл
        qsm_noniso_T_task_boundaries_t boundaries(control_data[step_index]);
        double temp_in = boundaries.temperature;
        double vol_flow = boundaries.volumetric_flow;

        // Выходные (расчетные) данные модели, которые выведутся в файл
        // температуры по двум способа расчета и их отклонения от факта
        double temp_out_calc = layer.temperature.back();
        double temp_out_delta = etalon_temperature - layer.temperature.back();
        double temp_out_calc_shukhov = layer.temperature_shukhov.back();
        double temp_out_delta_shukhov = etalon_temperature - layer.temperature_shukhov.back();

        std::vector<std::vector<double>>values_to_write{
            { etalon_temperature },
            { temp_out_calc },
            { temp_out_delta },
            { temp_out_calc_shukhov },
            { temp_out_delta_shukhov },
            { temp_in },
            { vol_flow }
        };

        print_t_profiles<std::string>(boundary_temperatures_file,
            static_cast<time_t>(0),
            std::vector<std::string>{ UnixToString(t) },
            values_to_write,
            "time,time,etalon,temp_out_calc_advection,temp_out_delta_advection,temp_out_calc_shukhov,temp_out_delta_shukhov,temp_in,volum");


    }

private:
    /// @brief Вывод профилей в файл
    /// @tparam OX Тип данных по оси ОX
    /// @param time_moment Время моделирования
    /// @param x Ось ОХ
    /// @param params Вектор параметров для вывода
    /// @param params_name Названия параметров для вывода
    /// @param filename Название файла
    template <typename OX>
    static void print_t_profiles(
        std::ofstream& output_file,
        const time_t& time_moment, 
        const std::vector<OX>& x, 
        const std::vector<std::vector<double>>& params, 
        const std::string& params_name)
    {
        static bool header_written = false;

        if (!header_written && !params_name.empty()) {
            output_file << params_name << std::endl;
            header_written = true;
        }

        size_t profCount = params.size();
        std::string time_str = UnixToString(time_moment);
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
        //output_file.close();
    };
};

/// @brief Класс, содержащий функции для вывода
/// профилей в файл в формате для плоттеров на Python
struct python_pressure_printer {
    /// @brief Вывод профилей в файл
    /// @tparam OX Тип данных по оси ОX
    /// @param time_moment Время моделирования
    /// @param x Ось ОХ
    /// @param params Вектор параметров для вывода
    /// @param params_name Названия параметров для вывода
    /// @param filename Название файла
    template <typename OX>
    static void print_p_profiles(const time_t& time_moment, const std::vector<OX>& x, const std::vector<std::vector<double>>& params, const std::string& params_name, const std::string& filename) {
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
        const std::string& folder,
        const time_t& t,
        const pipe_noniso_properties_t& pipe,
        const qsm_noniso_TP_layer& layer,
        const std::vector<double>& etalon_values = {}
    ) {
        print_p_profiles<double>(
            t,
            pipe.profile.coordinates,
            { layer.temperature, layer.pressure },
            "time,coordinates,temperature,pressure",
            folder + "results.csv"
        );

        // Вывод временного ряда отклонения расчётного температуры от фактического
        if (!etalon_values.empty())
        {
            double etalon = etalon_values[0];
            double calculated = layer.temperature.back();
            //double calculated_shukhov = layer.temp_shukhov.back();
            double calculated_p = layer.pressure.back();
            double pressure_delta = etalon_values[0] - layer.pressure.back();
            print_p_profiles<std::string>(static_cast<time_t>(0),
                std::vector<std::string>{ UnixToString(t) },
                std::vector<std::vector<double>>{ { calculated }, /*{calculated_shukhov},*/ {etalon}, {calculated_p}, {pressure_delta} },
                "time,time,calculated,etalon,calculated_p,diff_press",
                folder + "diff_temp.csv");
        }
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
    std::string path;
    // Временные ряды краевых условий
    std::vector<std::pair<std::vector<time_t>, std::vector<double>>> tag_data;
    // Временные ряды эталонных данных давления
    std::vector<std::pair<std::vector<time_t>, std::vector<double>>> etalon_tag_data;
    // Путь к реальным данным с трубопровода
    std::string folder = "../research_out/data/";

    /// @brief Подготовка к расчету для семейства тестов
    virtual void SetUp() override {
        using namespace std;

        // Указываем имя файла и желаемый шаг новой сетки
        std::string file_name = folder + "coord_heights.csv";
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

    }

    static std::tuple<
        std::vector<double>, 
        std::vector<std::vector<double>>, 
        std::vector<double>
    >
        prepare_real_data(double dt, const std::string& path_to_real_data)
    {
        using namespace std;
        std::vector<std::pair<std::string, std::string>>parameters =
        {
            { path_to_real_data + "Q_in", "m3/h-m3/s"s },
            { path_to_real_data + "t_in_nps", "C"s },
            { path_to_real_data + "p_in", "MPa"s },
            { path_to_real_data + "p_out", "MPa"s}

        };

        // Задаём период
        std::string start_period = "01.08.2021 00:00:00";
        std::string end_period = "01.09.2021 00:00:00";
        return prepare_timeseries_data(dt, start_period, end_period, parameters);
    }
};

/// @brief Тесты для расчёта на реальных данных
class qsmNonisoNoParties : public ::testing::Test {
protected:
    static std::tuple<std::vector<double>, std::vector<std::vector<double>>, std::vector<double>>
        prepare_real_data(double dt, const std::string& path_to_real_data)
    {
        using namespace std;
        std::vector<std::pair<std::string, std::string>>parameters =
        {
            { path_to_real_data + "Q_in", "m3/h-m3/s"s },
            { path_to_real_data + "t_in_nps", "C"s },
            { path_to_real_data + "t_out_n", "C"s}
        };

        // Задаём период
        std::string start_period = "01.08.2021 00:00:00";
        std::string end_period = "01.09.2021 00:00:00";

        return prepare_timeseries_data(dt, start_period, end_period, parameters);
    }
    /// @brief Зачитыват параметры трубы, которые лежат по заданному пути,
    /// и параметризует настроечный параметры Kто и тип модели (динамика, Шухов и проч.)
    /// @param pipe_data_path 
    /// @param ambient_heat_transfer 
    /// @return 
    static qsm_noniso_T_task_t prepare_task(const std::string& pipe_data_path,
        double ambient_heat_transfer,
        noniso_qsm_model_type model_type
        ) 
    {
        // По-хорошему, параметры нефти oil должны зависеть от трубы. 
        // Либо из папки зачитываться, либо прямо тут хардкодом прописываться в зависимости от трубы
        oil_parameters_t oil = get_noniso_research_oil();

        pipe_noniso_properties_t pipe = get_research_pipe_heatmodel(pipe_data_path);
        qsm_noniso_T_task_t task(pipe, oil, model_type);

        return task;
    }
};

/// @brief Проводит моделирование беспартийной тепловой модели 
/// на длительном интервале на холодном ЛУ (cold_lu)
/// Выводится в файл выходная температура: модельная, эталонная, а также разница между ними
/// (см. python_temperature_printer)
TEST_F(qsmNonisoNoParties, ColdLU_DynamicTemperature)
{
    // Путь к реальным данным с Линейного участка трубопровода
    const std::string pipe_data_path = get_pipe_data_path("cold_lu");

    double ambient_heat_transfer = 1.3786917741689342; // результат идентификации
    qsm_noniso_T_task_t task = prepare_task(pipe_data_path, 
        ambient_heat_transfer, noniso_qsm_model_type::Dynamic);

    double dt = 60;
    auto [times, control_data, etalon_temperature] = prepare_real_data(dt, pipe_data_path);

    // Путь к результатам ресёрча
    std::string result_file = prepare_research_folder_for_qsm_model2() + "temperature.csv";

    // переписать python_temperature_printer по аналогии с nonisothermal_qsm_batch_Tout_collector_t
    // Реализовать у него не статический метод process_data, где и делать всю работу
    // Раскомментить строчки ниже
    python_temperature_printer printer(result_file,
        times, control_data, etalon_temperature, task.get_pipe());
    quasistatic_batch(task, times, control_data, &printer);
}
