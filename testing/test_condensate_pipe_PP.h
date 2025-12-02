#pragma once


namespace pde_solvers {

class CondensatePipePP : public ::testing::Test {
protected:
    void SetUp() override {
        // Arrange: Создаем простую равномерную трубу для тестов
        pde_solvers::pipe_profile_t profile;
        const size_t points_count = 10;
        const double length = 1000.0; // 1 км

        // Создаем равномерную сетку
        for (size_t i = 0; i < points_count; ++i) {
            double x = i * length / (points_count - 1);
            double h = 0.0; // Нулевой уклон для упрощения
            profile.coordinates.push_back(x);
            profile.heights.push_back(h);
        }

        pipe.profile = profile;
        pipe.wall.diameter = 0.5; // 500 мм
        pipe.wall.wallThickness = 0.0001; // 0.1 мм
        pipe.kinematic_viscosity = 1e-6; // 1 сСт
    }

    pde_solvers::condensate_pipe_properties_t pipe;
};

/// @brief Проверяет способность системы PP задачи рассчитывать валидный расход для заданных граничных условий по давлению
TEST_F(CondensatePipePP, CalculatesValidFlowRate_ForGivenPressureBoundaries) {
    // Arrange
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
    EXPECT_NEAR(layer.pressure.front(), 5e6, 1e-6);
    EXPECT_NEAR(layer.pressure.back(), 4e6, 1e-6);

    // Давление должно уменьшаться вдоль трубы
    for (size_t i = 1; i < layer.pressure.size(); ++i) {
        EXPECT_LT(layer.pressure[i], layer.pressure[i - 1]);
    }
}

/// @brief Проверяет способность системы восстанавливать исходный расход при использовании давления из PQ задачи в PP задаче
TEST_F(CondensatePipePP, RecoversOriginalFlowRate_WhenPQPressureUsedInPP) {
    // Arrange
    pde_solvers::condensate_pipe_PQ_task_t pq_task(pipe);
    pde_solvers::condensate_pipe_PP_task_t pp_task(pipe);

    auto pq_conditions = pde_solvers::condensate_pipe_PQ_task_boundaries_t::default_values();
    pq_conditions.pressure_in = 5e6;
    pq_conditions.volumetric_flow = 0.3;
    pq_conditions.density = 850.0;

    // Act: 1. Рассчитываем PQ задачу
    pq_task.solve(pq_conditions);
    double calculated_pressure_out = pq_task.get_current_layer().pressure.back();

    // Act: 2. Используем найденное давление в PP задаче
    auto pp_conditions = pde_solvers::condensate_pipe_PP_task_boundaries_t::default_values();
    pp_conditions.pressure_in = pq_conditions.pressure_in;
    pp_conditions.pressure_out = calculated_pressure_out;
    pp_conditions.density = pq_conditions.density;

    pp_task.solve(pp_conditions);
    double pp_flow = pp_task.get_current_layer().std_volumetric_flow;
    double relative = abs(pp_flow - pq_conditions.volumetric_flow) / pq_conditions.volumetric_flow;
    // Assert: Расход должен совпадать
    EXPECT_NEAR(relative, 0, 1e-5);
}

/// @brief Проверяет способность системы PP задачи поддерживать стабильность расхода при множественных шагах по времени
TEST_F(CondensatePipePP, MaintainsFlowRateStability_OverMultipleTimeSteps) {
    // Arrange
    pde_solvers::condensate_pipe_PP_task_t task(pipe);
    auto initial_conditions = pde_solvers::condensate_pipe_PP_task_boundaries_t::default_values();
    initial_conditions.pressure_in = 5e6;
    initial_conditions.pressure_out = 4.5e6;
    initial_conditions.density = 850.0;

    const int steps = 10;
    const double dt = 10.0;
    std::vector<double> flow_history;

    // Act: Начальный расчет
    task.solve(initial_conditions);
    flow_history.push_back(task.get_current_layer().std_volumetric_flow);

    // Act: Много шагов по времени
    for (int i = 0; i < steps; ++i) {
        task.step(dt, initial_conditions);
        flow_history.push_back(task.get_current_layer().std_volumetric_flow);
    }

    // Assert
    // Проверяем, что расход стабилизировался
    double max_variation = 0.0;
    for (int i = 1; i < flow_history.size(); ++i) {
        max_variation = std::max(max_variation,
            std::abs(flow_history[i] - flow_history[i - 1]));
    }

    EXPECT_LT(max_variation, 1e-6);
}

/// @brief Проверяет способность системы уменьшать расход с увеличением плотности при фиксированном перепаде давления
TEST_F(CondensatePipePP, DecreasesFlowRate_WithIncreasingDensity_AtFixedPressureDrop) {
    // Arrange
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

/// @brief Проверяет способность системы производить нулевой расход при нулевом перепаде давления
TEST_F(CondensatePipePP, ProducesZeroFlowRate_WhenPressureDropIsZero) {
    // Arrange
    pde_solvers::condensate_pipe_PP_task_t task(pipe);
    auto conditions = pde_solvers::condensate_pipe_PP_task_boundaries_t::default_values();
    conditions.pressure_in = 5e6;
    conditions.pressure_out = 5e6; 
    conditions.density = 850.0;

    // Act
    task.solve(conditions);
    auto& layer = task.get_current_layer();

    // Assert
    EXPECT_EQ(layer.std_volumetric_flow, 0);
}

}

