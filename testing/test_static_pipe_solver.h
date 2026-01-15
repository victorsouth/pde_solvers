#pragma once


TEST(Static_Hydraulic_Solver, UseCase)
{
    typedef profile_collection_t<2> layer_variables_type;
    typedef moc_solver<2>::specific_layer layer_moc_type;

    typedef composite_layer_t<layer_variables_type, layer_moc_type> composite_layer_type;

    ring_buffer_t<composite_layer_type> buffer(2, 3);

    pipe_properties_t pipe;
    pipe.profile.coordinates = { 0, 1000, 2000 };
    pipe.profile.heights = pipe.profile.capacity = std::vector<double>(pipe.profile.coordinates.size(), 0);

    oil_parameters_t oil;
    PipeModelPGConstArea pipeModel(pipe, oil);

    profile_wrapper<double, 2> start_layer(get_profiles_pointers(buffer.current().vars.point_double));

    double Pin = 5.2e5;
    double Pout = 5e5;

    double G = solve_pipe_PP(pipeModel, Pin, Pout, &start_layer);
}

/// @brief расчет разницы давления между началом и концом трубы в QP задаче
/// @param pipe 
/// @param conditions 
/// @return разница давления между началом и концом
inline double deltaP_for_solve_qP(
    const pde_solvers::iso_nonbarotropic_pipe_properties_t& pipe,
    const pde_solvers::iso_nonbarotropic_pipe_PQ_task_boundaries_t& conditions)
{

    iso_nonbarotropic_pipe_layer_t layer(pipe.profile.get_point_count());
    for (double& density : layer.density_std.value) {
        density = conditions.density;
    }
    layer.std_volumetric_flow = conditions.volumetric_flow;
    std::vector<double>& p_profile = layer.pressure;
    int euler_direction = +1; // Задаем направление для Эйлера

    iso_nonbaro_impulse_equation_t pipeModel(pipe, layer.density_std.value, conditions.volumetric_flow, euler_direction);

    solve_euler<1>(pipeModel, euler_direction, conditions.pressure_in, &p_profile);
    double pressure_drop = layer.pressure.front() - layer.pressure.back();

    return pressure_drop;
}

/// @brief расчет расхода в PP задаче
/// @param pipe 
/// @param conditions 
/// @param initial_Q_for_Newton начальное значение расхода при решении методом Ньютона
/// @return расход в PP задаче
inline double Q_for_solve_PP(
    const pde_solvers::iso_nonbarotropic_pipe_properties_t& pipe,
    const pde_solvers::iso_nonbarotropic_pipe_PP_task_boundaries_t& conditions,
    double initial_Q_for_Newton = 0.2)
{
    iso_nonbarotropic_pipe_layer_t layer(pipe.profile.get_point_count());
    for (double& density : layer.density_std.value) {
        density = conditions.density;
    }

    solve_condensate_PP<iso_nonbarotropic_pipe_PP_task_boundaries_t, iso_nonbarotropic_pipe_layer_t> test = solve_condensate_PP(pipe, conditions, layer);
    fixed_solver_parameters_t<1, 0, golden_section_search> parameters;
    parameters.residuals_norm = 0.1; // погрешность 0.1 Па
    parameters.argument_increment_norm = 0;
    parameters.residuals_norm_allow_early_exit = true;
    // Создание структуры для записи результатов расчета
    fixed_solver_result_t<1> result;
    fixed_newton_raphson<1>::solve_dense(test, { initial_Q_for_Newton }, parameters, &result);
    return result.argument;
}

/// @brief Проверяет способность системы производить нулевой перепад давления при нулевом расходе
TEST(CondensatePipeQPPde, ProducesZeroPressureDrop_WhenFlowRateIsZero) {

    //Arrange
    auto pipe = create_test_pipe_for_PQ();
    auto conditions = pde_solvers::iso_nonbarotropic_pipe_PQ_task_boundaries_t::default_values();
    conditions.volumetric_flow = 0; // нулевой расход

    // Act 
    double pressure_drop = deltaP_for_solve_qP(pipe, conditions);

    //// Assert
    ASSERT_NEAR(pressure_drop, 0, 1e-6);


}

/// @brief Проверяет способность системы увеличивать потери давления с ростом расхода
TEST(CondensatePipeQPPde, IncreasesPressureLoss_WithIncreasingFlowRate) {
    //Arrange
    auto pipe = create_test_pipe_for_PQ();
    pde_solvers::iso_nonbarotropic_pipe_PQ_task_t task(pipe);

    std::vector<double> flows = { 0.1, 0.3, 0.5, 0.7 }; // м³/с


    auto initial_conditions = pde_solvers::iso_nonbarotropic_pipe_PQ_task_boundaries_t::default_values();
    initial_conditions.pressure_in = 5e6;
    initial_conditions.density = 850.0;

    //Act
    std::vector<double> pressure_drops;
    for (double flow : flows) {
        initial_conditions.volumetric_flow = flow;
        double pressure_drop = deltaP_for_solve_qP(pipe, initial_conditions);
        pressure_drops.push_back(pressure_drop);
    }

    //Assert
    // Потери давления должны увеличиваться с ростом расхода
    for (size_t i = 1; i < pressure_drops.size(); ++i) {
        EXPECT_GT(pressure_drops[i], pressure_drops[i - 1]);
    }

}

/// @brief Проверяет способность системы производить нулевой расход при нулевом перепаде давления
TEST(CondensatePipePPPde, ProducesZeroFlowRate_WhenPressureDropIsZero) {
    // Arrange
    auto pipe = create_test_pipe_for_PP();
    auto conditions = pde_solvers::iso_nonbarotropic_pipe_PP_task_boundaries_t::default_values();
    conditions.pressure_in = 5e6;
    conditions.pressure_out = 5e6;
    conditions.density = 850.0;

    // Act
    double Q = Q_for_solve_PP(pipe, conditions);

    // Assert
    EXPECT_NEAR(Q, 0, 1e-6);

}


/// @brief Проверяет способность системы уменьшать расход с увеличением плотности при фиксированном перепаде давления
TEST(CondensatePipePPPde, DecreasesFlowRate_WithIncreasingDensity_AtFixedPressureDrop) {
    // Arrange
    auto pipe = create_test_pipe_for_PP();
    auto base_conditions = pde_solvers::iso_nonbarotropic_pipe_PP_task_boundaries_t::default_values();
    base_conditions.pressure_in = 5e6;
    base_conditions.pressure_out = 4e6;

    std::vector<double> densities = { 700.0, 850.0, 1000.0 };
    std::vector<double> calculated_flows;

    // Act
    for (double density : densities) {
        base_conditions.density = density;
        double Q = Q_for_solve_PP(pipe, base_conditions);
        calculated_flows.push_back(Q);
    }

    // Assert
    // Расход должен уменьшаться с увеличением плотности
    for (size_t i = 1; i < calculated_flows.size(); ++i) {
        EXPECT_LT(calculated_flows[i], calculated_flows[i - 1]);
    }
}

