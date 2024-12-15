#pragma once



inline pipe_properties_t prepare_pipe(const std::string& path)
{
    // Указываем имя файла и желаемый шаг новой сетки
    //string file_name = folder + "coord_heights.csv";
    //Желаемый шаг
    std::string folder = path + "coord_heights.csv";
    double desired_dx = 200;
    pipe_properties_t pipe;

    // Создаём новый профиль с постоянным шагом
    pipe.profile = pipe_profile_uniform::get_uniform_profile_from_csv(desired_dx, folder);
    pipe.wall.diameter = 1;

    return pipe;
};

inline std::tuple<vector<double>, vector<vector<double>>, vector<double>> prepare_real_data(const std::string& path_to_real_data)
{
    // Временные ряды краевых условий
    vector<pair<vector<time_t>, vector<double>>> control_tag_data;
    // Временные ряды эталонных данных
    vector<pair<vector<time_t>, vector<double>>> etalon_tag_data;
    // Задаём период
    string start_period = "01.08.2021 00:00:00";
    string end_period = "01.09.2021 00:00:00";

    vector<pair<string, string>>parameters =
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

    double step = 60;

    time_t start_period_time = max(control_parameters_time_series.get_start_date(), etalon_parameters_time_series.get_start_date());
    time_t end_period_time = min(control_parameters_time_series.get_end_date(), etalon_parameters_time_series.get_end_date());
    time_t duration = (end_period_time - start_period_time);

    size_t dots_count = static_cast<size_t>(ceil(duration / step) + 0.00001);

    vector<double>  times = vector<double>(dots_count);
    vector<vector<double>> control_data = vector<vector<double>>(dots_count);
    vector<double> etalon_pressure = vector<double>(dots_count);

    for (size_t i = 0; i < dots_count; i++)
    {
        times[i] = step * i;
        time_t t = start_period_time + static_cast<time_t>(times[i] + 0.5);

        control_data[i] = control_parameters_time_series(t);
        etalon_pressure[i] = etalon_parameters_time_series(t).front();

    };

    return std::make_tuple(std::move(times), std::move(control_data), std::move(etalon_pressure));
}

//void print_diff_pressure_before_after(const double d_before, const double d_after)
//{
//    python_printer printer;
//
//    printer.print_profiles<double>(static_cast<time_t>(0),
//        times,
//        vector<vector<double>>{ calc_vector_residuals(d_before), calc_vector_residuals(d_after) },
//        "time,time,diff_press_before,diff_press_after",
//        folder + "diff_press.csv");
//}


struct ident_parameters_isothermal_qsm {
    double diameter_adaptation{ 1 };
    double friction_adaptation{ 1 };
    void set_adaptation(const pipe_properties_t& nominal_pipe, pipe_properties_t* pipe_to_ident) const 
    {
        pipe_to_ident->wall.diameter = diameter_adaptation * nominal_pipe.wall.diameter;
        pipe_to_ident->wall.resistance_function_adaptation = friction_adaptation;
    }
    
};

struct ident_isothermal_qsm_pipe_settings {
    bool ident_diameter{false};
    bool ident_friction{false};
    ident_parameters_isothermal_qsm unpack_ident_parameters(const VectorXd& packed_ident_parameters) {
        ident_parameters_isothermal_qsm result;
        size_t index = 0;
        if (ident_diameter) {
            result.diameter_adaptation = packed_ident_parameters(index++);
        }
        if (ident_friction) {
            result.friction_adaptation = packed_ident_parameters(index++);
        }
        return result;
    }
    void check_parameters() const {
        // если ничего не задано или наоборот, все задано, кинуть исключение
        // throw std::runtime_error("");
    }
};


class ident_isothermal_qsm_pipe_diameter_t : public fixed_least_squares_function_t
{
    const ident_isothermal_qsm_pipe_settings settings;
    const pipe_properties_t pipe_nominal;
    pipe_properties_t pipe_to_ident;

    const vector<double>& times;
    const vector<vector<double>>& control_data;
    const vector<double>& etalon_pressure;
public:
    ident_isothermal_qsm_pipe_diameter_t(
        const ident_isothermal_qsm_pipe_settings& settings,
        const pipe_properties_t& pipe, const vector<double>& times, const vector<vector<double>>& control_data, const vector<double>& etalon_pressure)
        : settings(settings)
        , pipe_nominal{ pipe }
        , pipe_to_ident{ pipe }
        , times{ times }
        , control_data{ control_data }
        , etalon_pressure{ etalon_pressure }
    {
        settings.check_parameters();
    }
protected:
    vector<double> calc_vector_residuals(const ident_parameters_isothermal_qsm& ident_parameters)
    {
        isothermal_qsm_batch_Pout_collector_t collector(times);
        ident_parameters.set_adaptation(pipe_nominal, &pipe_to_ident);

        isothermal_quasistatic_PQ_task_t<quickest_ultimate_fv_solver> task(pipe_to_ident);
        isothermal_quasistatic_batch<quickest_ultimate_fv_solver, isothermal_qsm_batch_Pout_collector_t::layer_type>(
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
public:
    virtual VectorXd residuals(const VectorXd& d) override {

        ident_parameters_isothermal_qsm ident_parameters;
        ident_parameters.diameter_adaptation = d(0);


        vector<double> simulation_result = calc_vector_residuals(ident_parameters);

        Eigen::Map<VectorXd> result(simulation_result.data(), simulation_result.size());

        return result;
    }
    double ident(fixed_solver_result_t<-1>* result = nullptr, fixed_solver_result_analysis_t<-1>* analysis = nullptr) {
        fixed_solver_parameters_t<-1, 0, golden_section_search> parameters;

        fixed_solver_result_t<-1> result_local;
        if (result == nullptr) {
            result = &result_local;
        }

        VectorXd initial_d(1); 
        initial_d(0) = 1;
        fixed_optimize_gauss_newton::optimize(*this, initial_d, parameters, result, analysis);

        if (result->result_code == numerical_result_code_t::Converged) {
            return result->argument(0);
        }
        else {
            return std::numeric_limits<double>::quiet_NaN();
        }
    }
};

void print_j_d(const double d, const double J, const std::string folder)
{
    python_printer printer;
    printer.print_profiles<double>(static_cast<time_t>(0),
        { d },
        vector<vector<double>>{ {J} },
        "time,d,J",
        folder + "j.csv");
}


TEST(IdentIsothermalQSM, Diameter)
{
    // Путь к реальным данным с трубопровода
    std::string data_path = "../research/2024-08-quasistationary-with-real-data/data/";

    pipe_properties_t pipe = prepare_pipe(data_path);
    auto [times, control_data, etalon_pressure] = prepare_real_data(data_path);

    ident_isothermal_qsm_pipe_settings ident_settings;
    ident_settings.ident_diameter = true;

    ident_isothermal_qsm_pipe_diameter_t test_ident(ident_settings, pipe, times, control_data, etalon_pressure);

    fixed_solver_result_t<-1> result;
    fixed_solver_result_analysis_t<-1> analysis;

    double result_d = test_ident.ident(&result, &analysis);
}

TEST(IdentIsothermalQSM, PipeIdentificationWithPrinter)
{
    std::string path = prepare_research_folder_for_qsm_model();


    //string folder = prepare_research_folder_for_qsm_model();

    //VectorXd initial_d = VectorXd::Zero(1);
    //initial_d(0) = 1;

    //// Путь к реальным данным с трубопровода
    //std::string data_path = "../research/2024-08-quasistationary-with-real-data/data/";
    //pipe_properties_t pipe = prepare_pipe(data_path);
    //auto [times, control_data, etalon_pressure] = prepare_real_data(data_path);

    //ident_isothermal_qsm_pipe_diameter_t test_ident;
    ////VectorXd residuals = test_ident.residuals(initial_d);



    //test_ident.print_diff_pressure_before_after(initial_d(0), result.argument(0));

}