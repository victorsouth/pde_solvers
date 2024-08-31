#pragma once

class quasistatic_printer {
public:
    template <typename OX>
    static void print_profiles(time_t time_moment, vector<OX> x, vector<vector<double>> params, std::string params_name, std::string filename) {
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

    static void print_all(time_t t, pipe_properties_t pipe, const density_viscosity_quasi_layer<false>& layer, std::string folder) {
        print_profiles<double>(
            t,
            pipe.profile.coordinates,
            { layer.density, layer.viscosity, layer.pressure },
            "time,coordinates,density,viscosity,pressure",
            folder + "results.csv"
        );
    }
};

/// @brief Тесты для солвера
class RealPipeline : public ::testing::Test {
protected:
    /// @brief Параметры трубы
    pipe_properties_t pipe;

    std::string folder = "../research/2024-08-quasistationary-with-real-data/data/";

    /// @brief Подготовка к расчету для семейства тестов
    virtual void SetUp() override {
        // Указываем имя файла и желаемый шаг новой сетки
        string file_name = folder + "coord_heights.csv";
        //Желаемый шаг
        double desired_dx = 200;

        // Создаём новый профиль с постоянным шагом
        pipe.profile = pipe_profile_uniform::get_uniform_profile_from_csv(desired_dx, file_name);
        pipe.wall.diameter = 1;
    }

    /// @brief Перегрузка функции для возможности задания постоянного
    /// шага по времени без эталонных данных
    template <typename Solver, typename Printer>
    void perform_quasistatic_simulation(
        const string& path,
        const pipe_properties_t& pipe,
        const vector_timeseries_t& boundary_timeseries,
        double dt = std::numeric_limits<double>::quiet_NaN())
    {
        perform_quasistatic_simulation<Solver>(path, pipe, boundary_timeseries, vector_timeseries_t(), dt);
    }

    
     /// @brief Стационарный расчет (с помощью initial boundaries),
    /// а затем квазистационарный расчет по краевым условиям (boundary_timeseries)
    /// @tparam Solver Численный метод расчета движения партий
    /// @tparam Printer Класс для вывода результатов в файл
    /// @param path Путь к файлу с результатом
    /// @param pipe Модель трубы
    /// @param boundary_timeseries Краевые условия
    /// !!! Важно, чтобы вектор на заданный момент времени был совместим по порядку параметров с isothermal_quasistatic_task_boundaries_t !!!
    /// @param etalon_timeseries Эталонные данные давления в конце трубопровода 
    /// @param dt Шаг по времени либо задаётся постоянным, 
    /// либо рассчитывается на каждом шаге моделирования для Cr = 1
    template <typename Solver, typename Printer>
    void perform_quasistatic_simulation(
        const string& path, 
        const pipe_properties_t& pipe,
        const vector_timeseries_t& boundary_timeseries, 
        //const vector<pair<vector<time_t>, vector<double>>>& params_data,
        const vector_timeseries_t& etalon_timeseries = vector_timeseries_t(),
        double dt = std::numeric_limits<double>::quiet_NaN())
    {

        time_t t = boundary_timeseries.get_start_date(); // Момент времени начала моделирования
        isothermal_quasistatic_task_boundaries_t initial_boundaries(boundary_timeseries(t));

        isothermal_quasistatic_task_t<Solver> task(pipe);
        task.solve(initial_boundaries);


        // Печатаем профиль трубы и первый слой к нему
        // task.print_profile(path);
        write_profile(pipe.profile, path + "pipe_coord_heights");
        //task.print_all(t, path);
        Printer::print_all(t, pipe, task.buffer.current(), path);

        do
        {
            // Интерполируем значения параметров в заданный момент времени
            vector<double> values_in_time_model = boundary_timeseries(t);
            isothermal_quasistatic_task_boundaries_t boundaries(values_in_time_model);

            double time_step = dt;
            if (std::isnan(time_step)) {
                double v = boundaries.volumetric_flow / pipe.wall.getArea(); 
                time_step = task.get_time_step_assuming_max_speed(v);
            }
            t += static_cast<time_t>(time_step);

            task.step(time_step, boundaries);
            Printer::print_all(t, pipe, task.buffer.current(), path);

            if (!etalon_timeseries.data.empty())
            {
                double pressure_delta = etalon_timeseries(t)[0] - task.buffer.current().pressure.back();
                Printer::print_profiles<std::string>(static_cast<time_t>(0),
                    vector<string>{ UnixToString(t) },
                    vector<vector<double>>{ { pressure_delta  } },
                    "time,time,diff_press",
                    path + "diff_press.csv");
            }
        } while (t <= boundary_timeseries.get_end_date());
    }
};

TEST(Test, QuasiStationaryFull)
{
    vector<int> v1 = { 1, 2, 3, 4, 5 };
    vector<int> v2;
    v2 = vector<int>(v1.begin() + 1, v1.begin() + 5);

    cout << v2.size() << endl;
}


/// @brief Пример использования метода Quickest Ultimate с гидравлическим расчетом  
TEST_F(RealPipeline, QuasiStationaryFull)
{
    // Создаём папку с результатами и получаем путь к ней
    string path = prepare_research_folder_for_qsm_model();

    vector<pair<string, string>>parameters =
    {
        { folder + "Q_in", "m3/h-m3/s"s },
        { folder + "p_in", "MPa"s },
        { folder + "rho_in", "kg/m3"s },
        { folder + "visc_in", "mm^2/s-m^2/s"s },
        { folder + "p_out", "MPa"s}

    };

    vector<pair<string, string>> etalon_parameters =
    {
        { folder + "p_out", "MPa"s}
    };

    // Задаём период
    string start_period = "01.08.2021 00:00:00";
    string end_period = "01.09.2021 00:00:00";

    // Считываем временные ряды параметров
    csv_multiple_tag_reader tags(parameters);
    auto tag_data = tags.read_csvs(start_period, end_period);

    // Помещаем временные ряды в вектор
    vector_timeseries_t params(tag_data);

    // Считываем временные ряды параметров
    csv_multiple_tag_reader etalon_tags(etalon_parameters);
    auto etalon_tag_data = tags.read_csvs(start_period, end_period);

    // Помещаем временные ряды в вектор
    vector_timeseries_t etalon_params(etalon_tag_data);

    perform_quasistatic_simulation<advection_moc_solver, quasistatic_printer>(
        path, pipe, params, etalon_params
    );
}