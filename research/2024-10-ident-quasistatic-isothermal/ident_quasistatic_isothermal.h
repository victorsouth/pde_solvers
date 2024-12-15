#pragma once



template <typename LayerType>
class data_processor_t {
public:
    virtual void process_data(size_t step_index, const LayerType& layer) = 0;
};


class isothermal_qsm_simulation_result_collector_t
    : public data_processor_t<density_viscosity_quasi_layer<true>>
{
public:
    typedef density_viscosity_quasi_layer<true> layer_type;
protected:
    vector<double> pipe_pressure_out;
public:
    isothermal_qsm_simulation_result_collector_t(const vector<double>& times)
        : pipe_pressure_out(times.size(), std::numeric_limits<double>::quiet_NaN())
    {

    }
    virtual void process_data(size_t step_index, 
        const density_viscosity_quasi_layer<true>& layer) override
    {
        // at() - проверяет выход за границы массива
        //pipe_pressure_out.at(step_index) = layer.pressure.back();
        pipe_pressure_out[step_index] = layer.pressure.back();
    }
    const vector<double>& get_pressure_out_calculated() const {
        return pipe_pressure_out;
    }
};

template <typename Solver, typename LayerType>
inline void perform_quasistatic_simulation(
    isothermal_quasistatic_task_t<Solver>& task,
    const vector<double>& times,
    const vector<vector<double>>& boundary_timeseries,
    data_processor_t<LayerType>* data_processor
)
{
    isothermal_quasistatic_task_boundaries_t initial_boundaries(boundary_timeseries[0]);
    task.solve(initial_boundaries);
    data_processor->process_data(0, task.get_buffer().current());

    for (size_t step_index = 1; step_index < times.size(); step_index++)
    {
        double time_step = times[step_index] - times[step_index - 1];
        isothermal_quasistatic_task_boundaries_t boundaries(boundary_timeseries[step_index]);

        task.step(time_step, boundaries);

        data_processor->process_data(step_index, task.get_buffer().current());
    }
};

TEST(TimeSeries, PrepareTimeSeries)
{
    // Параметры трубы
    pipe_properties_t pipe;
    // Путь к результатам ресёрча
    string path;
    // Путь к реальным данным с трубопровода
    std::string folder = "../research/2024-08-quasistationary-with-real-data/data/";
    // Временные ряды краевых условий
    vector<pair<vector<time_t>, vector<double>>> control_tag_data;
    // Временные ряды эталонных данных
    vector<pair<vector<time_t>, vector<double>>> etalon_tag_data;
    // Задаём период
    string start_period = "02.08.2021 00:00:00";
    string end_period = "02.08.2021 02:00:00";

    // Указываем имя файла и желаемый шаг новой сетки
    string file_name = folder + "coord_heights.csv";
    //Желаемый шаг
    double desired_dx = 200;

    // Создаём новый профиль с постоянным шагом
    pipe.profile = pipe_profile_uniform::get_uniform_profile_from_csv(desired_dx, file_name);
    pipe.wall.diameter = 1;

    vector<pair<string, string>>parameters =
    {
        { folder + "Q_in", "m3/h-m3/s"s },
        { folder + "p_in", "MPa"s },
        { folder + "rho_in", "kg/m3"s },
        { folder + "visc_in", "mm^2/s-m^2/s"s },
        { folder + "p_out", "MPa"s}

    };

    // Считываем временные ряды параметров
    csv_multiple_tag_reader tags(parameters);
    etalon_tag_data = { control_tag_data.back() };
    control_tag_data.pop_back();

    // Помещаем временные ряды в вектор
    vector_timeseries_t control_parameters_time_series(control_tag_data);
    vector_timeseries_t etalon_parameters_time_series(etalon_tag_data);

    double step = 600;

    time_t start_period_time = max(control_parameters_time_series.get_start_date(), etalon_parameters_time_series.get_start_date());
    time_t end_period_time = min(control_parameters_time_series.get_end_date(), etalon_parameters_time_series.get_end_date());
    time_t duration = (end_period_time - start_period_time);

    size_t dots_count = static_cast<size_t>(ceil(duration / step) + 0.00001);

    vector<double> times(dots_count);
    vector<vector<double>> control_data(dots_count);
    vector<vector<double>> etalon_pressure(dots_count);

    for (size_t i = 0; i < dots_count; i++)
    {
        times[i] = step * i;
        time_t t = start_period_time + static_cast<time_t>(times[i] + 0.5);

        control_data[i] = control_parameters_time_series(t);
        etalon_pressure[i] = control_parameters_time_series(t);

    }
    
};

pipe_properties_t read_profile_data()
{
    pipe_properties_t pipe;
    // Путь к реальным данным с трубопровода
    std::string folder = "../research/2024-08-quasistationary-with-real-data/data/";

    // Указываем имя файла и желаемый шаг новой сетки
    string file_name = folder + "coord_heights.csv";
    //Желаемый шаг
    double desired_dx = 200;

    // Создаём новый профиль с постоянным шагом
    pipe.profile = pipe_profile_uniform::get_uniform_profile_from_csv(desired_dx, file_name);
    pipe.wall.diameter = 1;
    return pipe;
};

//vector<vector<double>> read_parameters_data()
//{
//    // Путь к реальным данным с трубопровода
//    std::string folder = "../research/2024-08-quasistationary-with-real-data/data/";
//    // Временные ряды краевых условий
//    vector<pair<vector<time_t>, vector<double>>> control_tag_data;
//    // Временные ряды эталонных данных
//    vector<pair<vector<time_t>, vector<double>>> etalon_tag_data;
//    // Задаём период
//    string start_period = "01.08.2021 00:00:00";
//    string end_period = "02.08.2021 00:00:00";
//
//    vector<pair<string, string>>parameters =
//    {
//        { folder + "Q_in", "m3/h-m3/s"s },
//        { folder + "p_in", "MPa"s },
//        { folder + "rho_in", "kg/m3"s },
//        { folder + "visc_in", "mm^2/s-m^2/s"s },
//        { folder + "p_out", "MPa"s}
//
//    };
//
//    // Считываем временные ряды параметров
//    csv_multiple_tag_reader tags(parameters);
//    control_tag_data = tags.read_csvs(start_period, end_period);
//    etalon_tag_data = { control_tag_data.back() };
//    control_tag_data.pop_back();
//
//    // Помещаем временные ряды в вектор
//    vector_timeseries_t control_parameters_time_series(control_tag_data);
//    vector_timeseries_t etalon_parameters_time_series(etalon_tag_data);
//
//    double step = 60;
//
//    time_t start_period_time = max(control_parameters_time_series.get_start_date(), etalon_parameters_time_series.get_start_date());
//    time_t end_period_time = min(control_parameters_time_series.get_end_date(), etalon_parameters_time_series.get_end_date());
//    time_t duration = (end_period_time - start_period_time);
//
//    size_t dots_count = static_cast<size_t>(ceil(duration / step) + 0.00001);
//    vector<double> times(dots_count);
//    vector<vector<double>> control_data(dots_count);
//    vector<double> etalon_pressure(dots_count);
//
//    for (size_t i = 0; i < dots_count; i++)
//    {
//        times[i] = step * i;
//        time_t t = start_period_time + static_cast<time_t>(times[i] + 0.5);
//
//        control_data[i] = control_parameters_time_series(t);
//        etalon_pressure[i] = etalon_parameters_time_series(t).front();
//
//    };
//}


class pipe_identification_t : public fixed_least_squares_function_t
{
    // Параметры трубы
    pipe_properties_t pipe;


    string folder = prepare_research_folder_for_qsm_model();

    vector<double> times;
    vector<vector<double>> control_data;
    vector<double> etalon_pressure;
public:
    pipe_identification_t()
    {
        // Путь к реальным данным с трубопровода
        std::string folder = "../research/2024-08-quasistationary-with-real-data/data/";
        // Временные ряды краевых условий
        vector<pair<vector<time_t>, vector<double>>> control_tag_data;
        // Временные ряды эталонных данных
        vector<pair<vector<time_t>, vector<double>>> etalon_tag_data;
        // Задаём период
        string start_period = "01.08.2021 00:00:00";
        string end_period = "01.09.2021 00:00:00";

        // Указываем имя файла и желаемый шаг новой сетки
        string file_name = folder + "coord_heights.csv";
        //Желаемый шаг
        double desired_dx = 200;

        // Создаём новый профиль с постоянным шагом
        pipe.profile = pipe_profile_uniform::get_uniform_profile_from_csv(desired_dx, file_name);
        pipe.wall.diameter = 1;

        vector<pair<string, string>>parameters =
        {
            { folder + "Q_in", "m3/h-m3/s"s },
            { folder + "p_in", "MPa"s },
            { folder + "rho_in", "kg/m3"s },
            { folder + "visc_in", "mm^2/s-m^2/s"s },
            { folder + "p_out", "MPa"s}

        };

        // Считываем временные ряды параметров
        csv_multiple_tag_reader tags(parameters);
        control_tag_data = tags.read_csvs(start_period, end_period);
        etalon_tag_data = { control_tag_data.back() };
        control_tag_data.pop_back();

        // Помещаем временные ряды в вектор
        vector_timeseries_t control_parameters_time_series(control_tag_data);
        vector_timeseries_t etalon_parameters_time_series(etalon_tag_data);

        double step = 60;

        time_t start_period_time = max(control_parameters_time_series.get_start_date(), etalon_parameters_time_series.get_start_date());
        time_t end_period_time = min(control_parameters_time_series.get_end_date(), etalon_parameters_time_series.get_end_date());
        time_t duration = (end_period_time - start_period_time);

        size_t dots_count = static_cast<size_t>(ceil(duration / step) + 0.00001);
        times = vector<double>(dots_count);
        control_data = vector<vector<double>>(dots_count);
        etalon_pressure = vector<double>(dots_count);

        for (size_t i = 0; i < dots_count; i++)
        {
            times[i] = step * i;
            time_t t = start_period_time + static_cast<time_t>(times[i] + 0.5);

            control_data[i] = control_parameters_time_series(t);
            etalon_pressure[i] = etalon_parameters_time_series(t).front();

        };
    }
public:

    void print_diff_pressure_before_after(const double d_before, const double d_after)
    {
        python_printer printer;

        printer.print_profiles<double>(static_cast<time_t>(0),
            times,
            vector<vector<double>>{ calc_vector_residuals(d_before), calc_vector_residuals(d_after) },
            "time,time,diff_press_before,diff_press_after",
            folder + "diff_press.csv");
    }

    vector<double> calc_vector_residuals(const double d)
    {
        isothermal_qsm_simulation_result_collector_t collector(times);
        pipe.wall.diameter = d;
        isothermal_quasistatic_task_t<quickest_ultimate_fv_solver> task(pipe);

        perform_quasistatic_simulation<quickest_ultimate_fv_solver, isothermal_qsm_simulation_result_collector_t::layer_type>(
            task,
            times,
            control_data,
            &collector
        );

        const vector<double>& calc_pressure = collector.get_pressure_out_calculated();
        vector<double> simulation_result(times.size());

        std::transform(calc_pressure.begin(), calc_pressure.end(), etalon_pressure.begin(), simulation_result.begin(),
            [](double etalon, double calc) { return etalon - calc; });


        return simulation_result;
    }

    void print_j_d(const double d, const double J)
    {

        python_printer printer;
        printer.print_profiles<double>(static_cast<time_t>(0),
            { d },
            vector<vector<double>>{ {J} },
            "time,d,J",
            folder + "j.csv");
    }

    VectorXd residuals(const VectorXd& d) {

        vector<double> simulation_result = calc_vector_residuals(d(0));

        Eigen::Map<VectorXd> result(simulation_result.data(), simulation_result.size());

        double J = result.squaredNorm() / 1e13;
        print_j_d(d(0), J);

        return result;
    }
};



TEST(OptimiseGaussNewton, PipeIdentification)
{
    VectorXd initial_d = VectorXd::Zero(1);
    initial_d(0) = 1;

    pipe_identification_t test_ident;
    //VectorXd residuals = test_ident.residuals(initial_d);

    fixed_optimizer_parameters_t parameters;
    fixed_optimizer_result_t result;

    fixed_optimize_gauss_newton::optimize(test_ident, initial_d, parameters, &result);

}

TEST(OptimiseGaussNewton, PipeIdentificationWithPrinter)
{
    std::string path = prepare_research_folder_for_qsm_model();
    VectorXd initial_d = VectorXd::Zero(1);
    initial_d(0) = 1;

    pipe_identification_t test_ident;
    //VectorXd residuals = test_ident.residuals(initial_d);

    fixed_optimizer_parameters_t parameters;
    fixed_optimizer_result_t result;

    fixed_optimize_gauss_newton::optimize(test_ident, initial_d, parameters, &result);

    test_ident.print_diff_pressure_before_after(initial_d(0), result.argument(0));

}