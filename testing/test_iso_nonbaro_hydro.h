#pragma once


/// @brief расчет разницы давления между началом и концом трубы в QP задаче
/// @param pipe 
/// @param conditions 
/// @return разница давления между началом и концом
inline double deltaP_for_solve_qP(
    const pde_solvers::iso_nonbaro_pipe_properties_t& pipe,
    const pde_solvers::iso_nonbarotropic_pipe_PQ_task_boundaries_t& conditions)
{

    iso_nonbaro_pipe_layer_t layer(pipe.profile.get_point_count());
    for (double& density : layer.density_std.value) {
        density = conditions.density;
    }
    layer.std_volumetric_flow = conditions.volumetric_flow;
    std::vector<double>& p_profile = layer.pressure;
    int euler_direction = +1; // Задаем направление для Эйлера

    iso_nonbaro_impulse_equation_t pipeModel(pipe, layer, conditions.volumetric_flow, euler_direction);

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
    const pde_solvers::iso_nonbaro_pipe_properties_t& pipe,
    const pde_solvers::iso_nonbaro_pipe_PP_task_boundaries_t& conditions,
    double initial_Q_for_Newton = 0.2)
{
    iso_nonbaro_pipe_layer_t layer(pipe.profile.get_point_count());
    for (double& density : layer.density_std.value) {
        density = conditions.density;
    }

    rigorous_impulse_solver_PP<iso_nonbaro_impulse_equation_t>
        solver(pipe, layer, conditions.pressure_in, conditions.pressure_out);
    return solver.solve(initial_Q_for_Newton);
}

/// @brief Проверяет способность системы производить нулевой перепад давления при нулевом расходе
TEST(IsoNonbaroPipeQPPde, ProducesZeroPressureDrop_WhenFlowRateIsZero) {

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
TEST(IsoNonbaroPipeQPPde, IncreasesPressureLoss_WithIncreasingFlowRate) {
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
TEST(IsoNonbaroPipePPPde, ProducesZeroFlowRate_WhenPressureDropIsZero) {
    // Arrange
    auto pipe = create_test_pipe_for_PP();
    auto conditions = pde_solvers::iso_nonbaro_pipe_PP_task_boundaries_t::default_values();
    conditions.pressure_in = 5e6;
    conditions.pressure_out = 5e6;
    conditions.density = 850.0;

    // Act
    double Q = Q_for_solve_PP(pipe, conditions);

    // Assert
    EXPECT_NEAR(Q, 0, 1e-6);

}


/// @brief Проверяет способность системы уменьшать расход с увеличением плотности при фиксированном перепаде давления
TEST(IsoNonbaroPipePPPde, DecreasesFlowRate_WithIncreasingDensity_AtFixedPressureDrop) {
    // Arrange
    auto pipe = create_test_pipe_for_PP();
    auto base_conditions = pde_solvers::iso_nonbaro_pipe_PP_task_boundaries_t::default_values();
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

/// @brief Проверяет монотонное увеличение потерь давления с ростом кинематической вязкости при фиксированном расходе
TEST(IsoNonbaroPipeQPPde, IncreasesPressureLoss_WithIncreasingViscosity) {
    // Arrange: подготавливаем трубу, граничные условия и диапазон кинематических вязкостей
    auto pipe = create_test_pipe_for_PQ();
    auto conditions = pde_solvers::iso_nonbarotropic_pipe_PQ_task_boundaries_t::default_values();
    conditions.volumetric_flow = 0.3;
    conditions.pressure_in = 5e6;
    conditions.density = 850.0;

    std::vector<double> viscosities = { 1e-7, 4e-7, 8e-7, 2e-6 };
    std::vector<double> pressure_drops;

    // Act: рассчитываем потери давления для каждого значения вязкости
    for (double nu : viscosities) {
        pipe.kinematic_viscosity = nu;
        double pressure_drop = deltaP_for_solve_qP(pipe, conditions);
        pressure_drops.push_back(pressure_drop);
    }

    // Assert: потери давления должны строго увеличиваться с ростом вязкости
    for (size_t i = 1; i < pressure_drops.size(); ++i) {
        EXPECT_GT(pressure_drops[i], pressure_drops[i - 1]);
    }
}

/// @brief Проверяет монотонное уменьшение расхода с ростом кинематической вязкости при фиксированном перепаде давления
TEST(IsoNonbaroPipePPPde, DecreasesFlowRate_WithIncreasingViscosity_AtFixedPressureDrop) {
    // Arrange: подготавливаем трубу, фиксированный перепад давления и диапазон кинематических вязкостей
    auto pipe = create_test_pipe_for_PP();
    auto conditions = pde_solvers::iso_nonbaro_pipe_PP_task_boundaries_t::default_values();
    conditions.pressure_in = 5e6;
    conditions.pressure_out = 4e6;
    conditions.density = 850.0;

    std::vector<double> viscosities = { 1e-7, 4e-7, 8e-7, 2e-6 };
    std::vector<double> calculated_flows;

    // Act: рассчитываем расход методом Ньютона для каждого значения вязкости
    for (double nu : viscosities) {
        pipe.kinematic_viscosity = nu;
        double Q = Q_for_solve_PP(pipe, conditions);
        calculated_flows.push_back(Q);
    }

    // Assert: расход должен строго убывать с ростом вязкости
    for (size_t i = 1; i < calculated_flows.size(); ++i) {
        EXPECT_LT(calculated_flows[i], calculated_flows[i - 1]);
    }
}

/// @brief Проверяет малую чувствительность перепада давления к вариациям вязкости (+-30%) в турбулентном режиме
TEST(IsoNonbaroPipeQPPde, HasLowSensitivityToViscosityVariation_WhenTurbulent) {
    // Arrange: номинальная рабочая вязкость конденсата и вариации на +-30%
    auto pipe = create_test_pipe_for_PQ();
    auto conditions = pde_solvers::iso_nonbarotropic_pipe_PQ_task_boundaries_t::default_values();
    conditions.volumetric_flow = 0.3;
    conditions.pressure_in = 5e6;
    conditions.density = 750.0;

    double nu_nominal = 8e-7;
    double nu_low = nu_nominal * 0.7;
    double nu_high = nu_nominal * 1.3;

    // Act: рассчитываем перепад давления для номинала и отклонений
    pipe.kinematic_viscosity = nu_nominal;
    double dp_nominal = deltaP_for_solve_qP(pipe, conditions);

    pipe.kinematic_viscosity = nu_low;
    double dp_low = deltaP_for_solve_qP(pipe, conditions);

    pipe.kinematic_viscosity = nu_high;
    double dp_high = deltaP_for_solve_qP(pipe, conditions);

    // Assert: относительное изменение перепада давления не превышает 5%, абсолютная разница менее 0.1 атм
    double rel_diff_low = std::abs(dp_low - dp_nominal) / dp_nominal;
    double rel_diff_high = std::abs(dp_high - dp_nominal) / dp_nominal;

    EXPECT_LT(rel_diff_low, 0.05);
    EXPECT_LT(rel_diff_high, 0.05);
    EXPECT_LT(std::abs(dp_high - dp_low), 1e4);
}

