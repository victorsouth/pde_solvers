#pragma once

namespace pde_solvers {
;

/// @brief Ќефть по умолчанию дл€ тепловых задач
inline oil_parameters_t get_noniso_default_oil()
{
    array<viscosity_data_point, 2> viscosity_data{
    viscosity_data_point{60 + KELVIN_OFFSET, 0.0000006},
    viscosity_data_point{10 + KELVIN_OFFSET, 0.00006}
    };

    oil_parameters_t oil;
    oil.viscosity = oil_viscosity_parameters_t(viscosity_data);

    oil.density.nominal_density = 940;
    oil.heat.HeatCapacity = 2000;
    oil.heat.internalHeatTransferCoefficient = 257;

    oil.heat.pourPointTemperature = KELVIN_OFFSET + 12;
    return oil;
}

/// @brief “рубопровод по умолчанию дл€ тепловых задач
inline pipe_noniso_properties_t  get_noniso_default_pipe(double length = 12000, double dx = 1000)
{
    pipe_noniso_properties_t  pipe;

    //pipe.profile = PipeProfile::create(120, 0, 120000, 150, 30, 10e6);
    pipe.profile = pipe_profile_t::create(static_cast<size_t>(0.5 + length / dx), 0, length, 150, 30, 10e6);
    pipe.wall.equivalent_roughness = 0.0001;

    // это диаметр внутренний
    pipe.wall.diameter = 0.7;
    // это толщина одной стенки, к внешнему надо прибавл€ть удвоенную толщину
    pipe.wall.wallThickness = 0.01;

    // здесь бы структурировать
    // “еплопроводность материалов трубопровода, ¬т*м-1* -1
    // сталь, первый слой изол€ции, второй слой изол€ции, защитный слой
    pipe.heat.thermalConductivity = { 40, 0.2, 0.035, 0.3 };

    // “олщина стенки, слоев изол€ции и защитного сло€, м
    // первый слой изол€ции, второй слой изол€ции, защитный слой
    pipe.heat.layerThickness = { 0.020, 0.003, 0.050, 0.005 };

    return pipe;
}



/// @brief «онированный трубопровод дл€ тепловых задач с равной длиной участков с разным грунтом
inline zoned_pipe_properties get_zoned_pipe(
    const simple_pipe_properties& spipe,
    const vector<thermophysical_properties_t>& soils,
    HeatModelVer model_version = HeatModelVer::V2)
{
    zoned_pipe_properties pipe;
    pipe.model_version = model_version;

    //pipe.profile = PipeProfile::create(120, 0, 120000, 150, 30, 10e6);
    pipe.profile = pipe_profile_t::create(
        static_cast<size_t>(0.5 + spipe.length / spipe.dx),
        0, spipe.length,
        0, spipe.elevation, 10e6);
    pipe.wall.equivalent_roughness = 0.0001;

    // это диаметр внутренний
    pipe.wall.diameter = spipe.diameter;
    // это толщина одной стенки, к внешнему надо прибавл€ть удвоенную толщину
    pipe.wall.wallThickness = 0.01;

    size_t index = 0;
    size_t n = pipe.profile.get_point_count();

    double sensor_step = 20e3; // датчики каждые 20 км
    int sensor_count = std::max(2, static_cast<int>(spipe.length / sensor_step + 0.5) + 1);
    sensor_count -= 2; // датчики на границах будут в любом случае
    vector<double> sensor_coordinates(sensor_count);
    for (size_t index = 0; index < sensor_coordinates.size(); ++index) {
        sensor_coordinates[index] = sensor_step * (index + 1);
    }
    pipe.init_pressure_sensors(sensor_coordinates);

    for (const auto& soil : soils) {
        primary_heat_zone_parameters_t zone;
        zone.coordinate_begin = index;
        index += n / soils.size();

        // “олщина стенки, слоев изол€ции и защитного сло€, м
        // первый слой изол€ции, второй слой изол€ции, защитный слой
        zone.isolation_thickness = { 0.020, 0.003, 0.050, 0.005 };
        // “еплопроводность материалов трубопровода, ¬т*м-1* -1
        // сталь, первый слой изол€ции, второй слой изол€ции, защитный слой
        zone.isolation_conductivity = { 40, 0.2, 0.035, 0.3 };

        zone.soil = soil;

        pipe.heat_zones.emplace_back(zone);
    }

    //pipe.eqheat_zones = pipe.get_equivalent_zones();

    pipe.init_zone_coeffs();


    return pipe;
}

/// @brief «онированный трубопровод дл€ тепловых задач с равной длиной участков с разным грунтом
inline zoned_pipe_properties get_zoned_pipe(
    const vector<thermophysical_properties_t>& soils,
    HeatModelVer model_version = HeatModelVer::V2,
    double length = 12000, double dx = 1000, double diameter = 0.7)
{
    simple_pipe_properties spipe;
    spipe.elevation = 0;
    spipe.diameter = diameter;
    spipe.length = length;
    spipe.dx = dx;

    return get_zoned_pipe(spipe, soils, model_version);
}

}