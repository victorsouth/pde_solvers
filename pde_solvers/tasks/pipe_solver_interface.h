#pragma once

namespace pde_solvers {

class pipe_solver_transport_interface_t {
public:
    virtual void transport_step(double dt, double volumetric_flow, const pde_solvers::endogenous_values_t& boundaries) = 0;
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

}

