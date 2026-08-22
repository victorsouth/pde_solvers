#pragma once

/// @brief Проверяет, что масса в трубе увеличивается при увеличении плотности нефти
TEST(NonBaroPipeWithMassAccounting, MassInPipeIncreases_WhenDensityIncreases) {
    // Arrange: подготавливаем трубу с поддержкой расчета массы
    auto pipe = pde_solvers::iso_nonbaro_pipe_mass_accounting_properties_t::default_values();
    pipe.calculate_mass = true;
    const size_t point_count = pipe.profile.get_point_count();
    pde_solvers::iso_nonbaro_pipe_mass_accounting_solver_t::buffer_type buffer(2, point_count);
    pde_solvers::iso_nonbaro_pipe_mass_accounting_solver_t solver(pipe, buffer);

    // Act: считаем массу при начальной плотности и после локального увеличения плотности
    pde_solvers::endogenous_values_t boundaries;
    boundaries.density_std.value = 850.0;
    boundaries.density_std.confidence = true;
    solver.transport_solve(1.0, boundaries);
    const double mass_850 = solver.calculate_mass();

    boundaries.density_std.value = 900.0;
    solver.transport_solve(1.0, boundaries);
    const double mass_mixed = solver.calculate_mass();

    // Assert: масса после увеличения плотности должна возрасти
    ASSERT_GT(mass_mixed, mass_850);
}
