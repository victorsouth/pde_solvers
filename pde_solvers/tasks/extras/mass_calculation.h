#pragma once

namespace pde_solvers {
;

/// @brief Вычисляет полную массу вещества в трубе по профилю плотности, кг
/// @details Стенка трубы абсолютно жесткая. Нефть несжимаемая, 
/// плотность не зависит от давления и температуры (т.е. переданная плотность трактуется как фактическая).
inline double calculate_isothermal_incompressible_fluid_mass_in_rigid_pipe(
    const std::vector<double>& coordinates,
    const std::vector<double>& density_profile,
    double cross_section_area)
{
    double fluid_mass_in_pipe = 0.0;
    for (size_t cell_index = 0; cell_index < coordinates.size() - 1; ++cell_index) {
        const double cell_density = density_profile[cell_index];
        if (!std::isfinite(cell_density)) {
            throw std::runtime_error(
                "Density is NaN at cell index " + std::to_string(cell_index) + " in pipe");
        }
        const double cell_length = coordinates[cell_index + 1] - coordinates[cell_index];
        fluid_mass_in_pipe += cell_density * cross_section_area * cell_length;
    }
    return fluid_mass_in_pipe;
}

/// @brief Вычисляет массу небаротропного флюида в жесткой трубе на основе сортовой плотности и давления
/// @details Формулу для расчета взята из [2025 - Южанин, Гришухин - Небаротроная модель] - формула (2)
/// @param nominal_pressure - номинальное давление, Па (по умолчанию 101325 Па - стандартные условия)
inline double calculate_iso_nonbaro_fluid_mass_in_rigid_pipe(
    const std::vector<double>& coordinates,
    const std::vector<double>& nominal_density_profile,
    const std::vector<double>& pressure_profile,
    double cross_section_area,
    double oil_compressibility,
    double nominal_pressure = 101325.0)
{  
    double fluid_mass_in_pipe = 0.0;
    for (size_t cell_index = 0; cell_index < coordinates.size() - 1; ++cell_index) {
        const double nominal_density = nominal_density_profile[cell_index];
        
        if (!std::isfinite(nominal_density)) {
            throw std::runtime_error(
                "Nominal density is NaN at cell index " + std::to_string(cell_index) + " in pipe");
        }

        const double pressure_in_cell = pressure_profile[cell_index];
        if (!std::isfinite(pressure_in_cell)) {
            throw std::runtime_error(
                "Pressure is NaN at cell index " + std::to_string(cell_index) + " in pipe");
        }
       
        const double working_density =
            nominal_density * (1.0 + oil_compressibility * (pressure_in_cell - nominal_pressure));
        const double cell_length = coordinates[cell_index + 1] - coordinates[cell_index];
        fluid_mass_in_pipe += working_density * cross_section_area * cell_length;
    }
    return fluid_mass_in_pipe;
}



}
