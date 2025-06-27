#pragma once

namespace pde_solvers {
;


/// @brief Нефть по умолчанию для тепловых задач
inline oil_parameters_t get_default_oil_heatmodel()
{
    std::array<viscosity_data_point, 2> viscosity_data{
    viscosity_data_point{60 + KELVIN_OFFSET, 0.0000006},
    viscosity_data_point{10 + KELVIN_OFFSET, 0.00006}
    };

    oil_parameters_t oil;
    oil.viscosity = oil_viscosity_parameters_t(viscosity_data);

    oil.density.nominal_density = 940;
    oil.heat.HeatCapacity = 2000;
    oil.heat.internalHeatTransferCoefficient = 257;

    oil.heat.pour_point_temperature = KELVIN_OFFSET + 12;
    return oil;
}


// @brief Создание модели трубопровода по реальным данным 
/// @param path Путь к файлу с профилем реального ЛУ
/// @return Модель трубы с профилем на основе профиля реального участка трубы 
inline pipe_noniso_properties_t get_research_pipe_heatmodel(const std::string& path)
{
    // Указываем имя файла
    std::string folder = path + "coord_heights.csv";
    //Желаемый шаг
    double desired_dx = 200;

    pipe_noniso_properties_t pipe;

    // Создаём новый профиль с постоянным шагом
    pipe.profile = pipe_profile_uniform::get_uniform_profile_from_csv(desired_dx, folder);
    // Номинальный диаметр трубы
    pipe.wall.diameter = 1;
    //Температура грунта по таблице
    pipe.heat.ambientTemperature = 286.15;

    return pipe;
};

/// @brief Трубопровод по умолчанию для тепловых задач
inline pipe_noniso_properties_t  get_default_pipe_heatmodel(double length = 12000, double dx = 1000)
{
    pipe_noniso_properties_t  pipe;

    //pipe.profile = PipeProfile::create(120, 0, 120000, 150, 30, 10e6);
    pipe.profile = pipe_profile_t::create(static_cast<size_t>(0.5 + length / dx), 0, length, 150, 30, 10e6);
    pipe.wall.equivalent_roughness = 0.0001;

    // это диаметр внутренний
    pipe.wall.diameter = 0.7;
    // это толщина одной стенки, к внешнему надо прибавлять удвоенную толщину
    pipe.wall.wallThickness = 0.01;

    // здесь бы структурировать
    // Теплопроводность материалов трубопровода, Вт*м-1*К-1
    // сталь, первый слой изоляции, второй слой изоляции, защитный слой
    pipe.heat.thermalConductivity = { 40, 0.2, 0.035, 0.3 };

    // Толщина стенки, слоев изоляции и защитного слоя, м
    // первый слой изоляции, второй слой изоляции, защитный слой
    pipe.heat.layerThickness = { 0.020, 0.003, 0.050, 0.005 };

    return pipe;
}



/// @brief Зонированный трубопровод для тепловых задач с равной длиной участков с разным грунтом
inline zoned_pipe_properties get_zoned_pipe_heatmodel(
    const simple_pipe_properties& spipe,
    const std::vector<thermophysical_properties_t>& soils,
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
    // это толщина одной стенки, к внешнему надо прибавлять удвоенную толщину
    pipe.wall.wallThickness = 0.01;

    size_t index = 0;
    size_t n = pipe.profile.get_point_count();

    double sensor_step = 20e3; // датчики каждые 20 км
    int sensor_count = std::max(2, static_cast<int>(spipe.length / sensor_step + 0.5) + 1);
    sensor_count -= 2; // датчики на границах будут в любом случае
    std::vector<double> sensor_coordinates(sensor_count);
    for (size_t index = 0; index < sensor_coordinates.size(); ++index) {
        sensor_coordinates[index] = sensor_step * (index + 1);
    }
    pipe.init_pressure_sensors(sensor_coordinates);

    for (const auto& soil : soils) {
        primary_heat_zone_parameters_t zone;
        zone.coordinate_begin = index;
        index += n / soils.size();

        // Толщина стенки, слоев изоляции и защитного слоя, м
        // первый слой изоляции, второй слой изоляции, защитный слой
        zone.isolation_thickness = { 0.020, 0.003, 0.050, 0.005 };
        // Теплопроводность материалов трубопровода, Вт*м-1*К-1
        // сталь, первый слой изоляции, второй слой изоляции, защитный слой
        zone.isolation_conductivity = { 40, 0.2, 0.035, 0.3 };

        zone.soil = soil;

        pipe.heat_zones.emplace_back(zone);
    }

    //pipe.eqheat_zones = pipe.get_equivalent_zones();

    pipe.init_zone_coeffs();


    return pipe;
}

/// @brief Зонированный трубопровод для тепловых задач с равной длиной участков с разным грунтом
inline zoned_pipe_properties get_zoned_pipe_heatmodel(
    const std::vector<thermophysical_properties_t>& soils,
    HeatModelVer model_version = HeatModelVer::V2,
    double length = 12000, double dx = 1000, double diameter = 0.7)
{
    simple_pipe_properties spipe;
    spipe.elevation = 0;
    spipe.diameter = diameter;
    spipe.length = length;
    spipe.dx = dx;

    return get_zoned_pipe_heatmodel(spipe, soils, model_version);
}

}