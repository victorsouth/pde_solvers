#pragma once

namespace pde_solvers {

class pipe_solver_transport_interface_t {
public:
    virtual void transport_step(double dt, double volumetric_flow, const pde_solvers::endogenous_values_t& boundaries) = 0;
    virtual pde_solvers::endogenous_values_t get_endogenous_output(double volumetric_flow) const = 0;
};

class pipe_solver_hydro_interface_t {
public:
    virtual double hydro_solve_PP(double pressure_input, double pressure_output) = 0;
    virtual double hydro_solve_QP(double volumetric_flow, double pressure_output) = 0;
    virtual double hydro_solve_PQ(double volumetric_flow, double pressure_in) = 0;
    virtual std::array<double, 2> hydro_solve_PP_jacobian(double pressure_input, double pressure_output) = 0;
};

class pipe_solver_interface_t
    : public pipe_solver_hydro_interface_t
    , public pipe_solver_transport_interface_t
{
};

/// @brief Traits для солвера, реализующего гидравлический интерфейс
template<typename T>
constexpr bool has_hydro_interface_v = std::is_base_of_v<pde_solvers::pipe_solver_hydro_interface_t, T>;

/// @brief Traits для солвера, реализующего только транспортный интерфейс (без гидравлического)
template<typename T>
constexpr bool is_transport_only_solver_v = 
    std::is_base_of_v<pde_solvers::pipe_solver_transport_interface_t, T> &&
    !std::is_base_of_v<pde_solvers::pipe_solver_hydro_interface_t, T>;

}
