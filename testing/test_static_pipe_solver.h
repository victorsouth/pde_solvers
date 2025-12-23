#pragma once

// для использования функции генерации трубы create_test_pipe_for_PQ(), create_test_pipe_for_PP()
#include "test_condensate_pipe_PP.h"
#include "test_condensate_pipe_PQ.h"
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



    /// @brief Проверяет способность системы производить нулевой перепад давления при нулевом расходе
TEST(CondensatePipeQP, ProducesZeroPressureDrop_WhenFlowRateIsZero) {

    //Arrange
    auto pipe = create_test_pipe_for_PQ();
    pde_solvers::condensate_pipe_PQ_task_t task(pipe);
    auto conditions = pde_solvers::condensate_pipe_PQ_task_boundaries_t::default_values();
    conditions.volumetric_flow = 0; // нулевой расход

    // Act 
    task.solve(conditions);

    // Assert
    auto& layer = task.get_current_layer();
    double pressure_drop = layer.pressure.front() - layer.pressure.back();
    ASSERT_NEAR(pressure_drop, 0, 1e-6);

}

/// @brief Проверяет способность системы увеличивать потери давления с ростом расхода
TEST(CondensatePipeQP, IncreasesPressureLoss_WithIncreasingFlowRate) {
    //Arrange
    auto pipe = create_test_pipe_for_PQ();
    pde_solvers::condensate_pipe_PQ_task_t task(pipe);

    std::vector<double> flows = { 0.1, 0.3, 0.5, 0.7 }; // м³/с


    auto initial_conditions = pde_solvers::condensate_pipe_PQ_task_boundaries_t::default_values();
    initial_conditions.pressure_in = 5e6;
    initial_conditions.density = 850.0;

    //Act
    std::vector<double> pressure_drops;
    for (double flow : flows) {
        initial_conditions.volumetric_flow = flow;
        task.solve(initial_conditions);

        auto& layer = task.get_current_layer();
        double pressure_drop = layer.pressure.front() - layer.pressure.back();
        pressure_drops.push_back(pressure_drop);
    }

    //Assert
    // Потери давления должны увеличиваться с ростом расхода
    for (size_t i = 1; i < pressure_drops.size(); ++i) {
        EXPECT_GT(pressure_drops[i], pressure_drops[i - 1]);
    }

}

/// @brief Проверяет способность системы производить нулевой расход при нулевом перепаде давления
TEST(CondensatePipePP, ProducesZeroFlowRate_WhenPressureDropIsZero) {
    // Arrange
    auto pipe = pde_solvers::create_test_pipe_for_PP();
    pde_solvers::condensate_pipe_PP_task_t task(pipe);
    auto conditions = pde_solvers::condensate_pipe_PP_task_boundaries_t::default_values();
    conditions.pressure_in = 5e6;
    conditions.pressure_out = 5e6;
    conditions.density = 850.0;

    // Act
    task.solve(conditions);
    auto& layer = task.get_current_layer();

    // Assert
    EXPECT_NEAR(layer.std_volumetric_flow, 0, 1e-6);

}


/// @brief Проверяет способность системы уменьшать расход с увеличением плотности при фиксированном перепаде давления
TEST(CondensatePipePP, DecreasesFlowRate_WithIncreasingDensity_AtFixedPressureDrop) {
    // Arrange
    auto pipe = create_test_pipe_for_PP();
    pde_solvers::condensate_pipe_PP_task_t task(pipe);
    auto base_conditions = pde_solvers::condensate_pipe_PP_task_boundaries_t::default_values();
    base_conditions.pressure_in = 5e6;
    base_conditions.pressure_out = 4e6;

    std::vector<double> densities = { 700.0, 850.0, 1000.0 };
    std::vector<double> calculated_flows;

    // Act
    for (double density : densities) {
        base_conditions.density = density;
        task.solve(base_conditions);
        calculated_flows.push_back(task.get_current_layer().std_volumetric_flow);
    }

    // Assert
    // Расход должен уменьшаться с увеличением плотности
    for (size_t i = 1; i < calculated_flows.size(); ++i) {
        EXPECT_LT(calculated_flows[i], calculated_flows[i - 1]);
    }
}


/// @brief Проверяет способность системы PP задачи рассчитывать валидный расход для заданных граничных условий по давлению
TEST(CondensatePipePP, CalculatesValidFlowRate_ForGivenPressureBoundaries) {
    // Arrange
    auto pipe = create_test_pipe_for_PP();
    pde_solvers::condensate_pipe_PP_task_t task(pipe);
    auto initial_conditions = pde_solvers::condensate_pipe_PP_task_boundaries_t::default_values();
    initial_conditions.pressure_in = 5e6;
    initial_conditions.pressure_out = 4e6;
    initial_conditions.density = 850.0;

    // Act
    task.solve(initial_conditions);
    auto& layer = task.get_current_layer();

    // Assert
    // Проверяем инициализацию плотности
    for (const auto& density : layer.density) {
        EXPECT_NEAR(density, 850.0, 1e-6);
    }

    // Проверяем, что расход был рассчитан
    EXPECT_FALSE(std::isnan(layer.std_volumetric_flow));
    EXPECT_GT(layer.std_volumetric_flow, 0.0);

    // Проверяем расчет давления
    ASSERT_FALSE(layer.pressure.empty());
    EXPECT_NEAR(layer.pressure.front(), 5e6, 0.1);
    EXPECT_NEAR(layer.pressure.back(), 4e6, 0.1);

    // Давление должно уменьшаться вдоль трубы
    for (size_t i = 1; i < layer.pressure.size(); ++i) {
        EXPECT_LT(layer.pressure[i], layer.pressure[i - 1]);
    }
}