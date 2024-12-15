#pragma once
#include <string>
#include <vector>
//#include "../timeseries/timeseries_helpers.h"

namespace pde_solvers {
;


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
    bool ident_diameter{ false };
    bool ident_friction{ false };
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

}
