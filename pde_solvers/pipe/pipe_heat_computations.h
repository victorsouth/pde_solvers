#pragma once


namespace pde_solvers {
;

inline double temperature_shukhov(
    const pipe_noniso_properties_t & pipe,
    const oil_parameters_t& oil,
    double Tleft, double mass_flow, double coordinate)
{
    // формула Шухова (Лурье-2012, формула 1.52),
    double T0 = Tleft;
    double S_0 = pipe.wall.getArea();


    double velocity = mass_flow / (S_0 * oil.density.nominal_density);
    double Re = velocity * pipe.wall.diameter / oil.viscosity();
    double hydraulic_resistance = hydraulic_resistance_shifrinson(Re, pipe.wall.relativeRoughness());


    double t_friction = 0.125 * hydraulic_resistance * oil.density.nominal_density * fabs(pow(velocity, 3)) / pipe.heat.ambient_heat_transfer;
    double t_stationary = pipe.heat.ambientTemperature + t_friction;

    return t_stationary + (T0 - t_stationary) *
        exp(-M_PI * pipe.wall.diameter * pipe.heat.ambient_heat_transfer * coordinate / (oil.heat.HeatCapacity * fabs(mass_flow)));
}


/// @brief Расчет распределения по Шухову по Лурье 2019, параграф 6.2 стр. 241
inline void compute_shuhov_temperature_distribution(
    const pipe_noniso_properties_t & pipe,
    const oil_parameters_t& oil,
    double mass_flow,
    double Tstart,
    const double& dx, std::vector<double>* _result)
{
    std::vector<double>& temperature = *_result;

    double d = pipe.wall.diameter;
    double S_0 = pipe.wall.getArea();
    double density = oil.density.nominal_density;
    double Kt = pipe.get_equivalent_heat_transfers(oil).first;
    double C = oil.heat.HeatCapacity;

    double v = mass_flow / (S_0 * density);
    double Re = v * pipe.wall.diameter / oil.viscosity();
    double lambda = hydraulic_resistance_shifrinson(Re, pipe.wall.relativeRoughness());

    // Расчет выделения тепла за счет внутреннего трения
    double Tx = lambda * density * pow(v, 3) / (8 * Kt);

    // Расчет распределения
    double Tstationary = pipe.heat.ambientTemperature + Tx;
    for (size_t i = 0; i < temperature.size(); i++)
    {
        //temperature[i] = Tstationary + (Tstart - Tstationary) * exp(-M_PI * d * Kt * i * dx / (mass_flow * C));
        temperature[i] = temperature_shukhov(pipe, oil, Tstart, mass_flow, i * dx);
    }
}

}