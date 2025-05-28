#pragma once

/// @brief Проверяет способность корректно инициализировать oil_viscosity_parameters_t 
/// по стандартной таблице вязкости
TEST(Viscosity, InitializesWithTable)
{
    using namespace pde_solvers;

    std::array<double, 3> visc_std_table{ 35.2166964842424e-6, 15.1959389818927e-6, 4.30720885400170e-6 };

    oil_viscosity_parameters_t viscosity(visc_std_table);

    auto coeffs = viscosity_table_model_t::reconstruct(visc_std_table);

    double v1 = viscosity_table_model_t::calc(coeffs, KELVIN_OFFSET + 30);
    double v2 = viscosity(KELVIN_OFFSET + 30);

    ASSERT_NEAR(v1, v2, 1e-12);
}