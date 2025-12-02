#pragma once


namespace pde_solvers {

class CondensatePipeQP : public ::testing::Test {
protected:
    void SetUp() override {
        // Создаем простую равномерную трубу для тестов
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
        pipe.wall.wallThickness = 10e-3; // 0.1 мм
        pipe.kinematic_viscosity = 1e-6; // 1 сСт
    }

    pde_solvers::condensate_pipe_properties_t pipe;
};


/// @brief Проверяет способность системы PQ задачи поддерживать корректный профиль давления при заданных начальных условиях.
/// Тест выполняет следующие проверки:
/// 1. Корректность инициализации плотности - все значения плотности в профиле должны соответствовать заданному начальному значению (800.0 кг/м³) с точностью 1e-6
/// 2. Наличие рассчитанного профиля давления - профиль давления не должен быть пустым
/// 3. Сохранение граничного условия на входе - давление на входе должно точно совпадать с заданным значением (5 МПа) с точностью 1e-6
/// 4. Монотонное убывание давления вдоль трубы - давление должно строго уменьшаться от входа к выходу из-за гидравлических потерь
/// 5. Физическая корректность значений - все значения давления в профиле должны быть положительными
TEST_F(CondensatePipeQP, MaintainsCorrectPressureProfile_WhenGivenInitialConditions) {

    //Arrange
    pde_solvers::condensate_pipe_PQ_task_t task(pipe);

    auto initial_conditions = pde_solvers::condensate_pipe_PQ_task_boundaries_t::default_values();
    initial_conditions.volumetric_flow = 0.3; // м³/с
    initial_conditions.pressure_in = 5e6; // 5 МПа
    initial_conditions.density = 800.0; // кг/м³

    // Act
    task.solve(initial_conditions);
    auto& layer = task.get_current_layer();

    //Assert
    
    // Проверяем инициализацию плотности
    for (const auto& density : layer.density) {
        EXPECT_NEAR(density, 800.0, 1e-6);
    }

    // Проверяем расчет давления
    ASSERT_FALSE(layer.pressure.empty());

    // Давление на входе должно быть равно заданному
    EXPECT_NEAR(layer.pressure.front(), 5e6, 1e-6);

    // Давление должно уменьшаться вдоль трубы из-за гидравлических потерь
    for (size_t i = 1; i < layer.pressure.size(); ++i) {
        EXPECT_LT(layer.pressure[i], layer.pressure[i - 1]);
    }

    // Проверяем, что давление положительное
    for (const auto& pressure : layer.pressure) {
        EXPECT_GT(pressure, 0.0);
    }
}

/// @brief Проверяет способность системы увеличивать потери давления с ростом расхода
TEST_F(CondensatePipeQP, IncreasesPressureLoss_WithIncreasingFlowRate) {
    //Arrange
    pde_solvers::condensate_pipe_PQ_task_t task(pipe);

    std::vector<double> flows = { 0.1, 0.3, 0.5, 0.7 }; // м³/с
    std::vector<double> pressure_drops;

    auto initial_conditions = pde_solvers::condensate_pipe_PQ_task_boundaries_t::default_values();
    initial_conditions.pressure_in = 5e6;
    initial_conditions.density = 850.0;
    
    //Act
    for (double flow : flows) {
        initial_conditions.volumetric_flow = flow;
        task.solve(initial_conditions);

        auto& layer = task.get_current_layer();
        double pressure_drop = layer.pressure.front() - layer.pressure.back();
        pressure_drops.push_back(pressure_drop);

        // Для отладки
        std::cout << "Flow: " << flow << " m³/s, Pressure drop: "
            << pressure_drop << " Pa" << std::endl;
    }

    //Assert
    // Потери давления должны увеличиваться с ростом расхода
    for (size_t i = 1; i < pressure_drops.size(); ++i) {
        EXPECT_GT(pressure_drops[i], pressure_drops[i - 1]);
    }

}

/// @brief Проверяет способность системы поддерживать стабильность давления при множественных шагах по времени
TEST_F(CondensatePipeQP, MaintainsPressureStability_OverMultipleTimeSteps) {
    //Arrange
    pde_solvers::condensate_pipe_PQ_task_t task(pipe);

    auto initial_conditions = pde_solvers::condensate_pipe_PQ_task_boundaries_t::default_values();
    task.solve(initial_conditions);

    // Выполняем много шагов и проверяем стабильность
    const int steps = 100;
    double dt = 1.0; // 1 секунда

    std::vector<double> inlet_pressures;
    std::vector<double> outlet_pressures;

    for (int i = 0; i < steps; ++i) {
        //Act
        task.step(dt, initial_conditions);
        auto& layer = task.get_current_layer();

        inlet_pressures.push_back(layer.pressure.front());
        outlet_pressures.push_back(layer.pressure.back());

        //Assert
        // Проверяем, что давление остается в физических пределах
        EXPECT_GT(layer.pressure.front(), 0.0);
        EXPECT_GT(layer.pressure.back(), 0.0);
    }

    //Assert
    // При постоянных граничных условиях давления должны стабилизироваться
    double inlet_variation = 0.0;
    double outlet_variation = 0.0;

    for (int i = 1; i < steps; ++i) {
        inlet_variation += std::abs(inlet_pressures[i] - inlet_pressures[i - 1]);
        outlet_variation += std::abs(outlet_pressures[i] - outlet_pressures[i - 1]);
    }

    inlet_variation /= (steps - 1);
    outlet_variation /= (steps - 1);

    // Вариации должны быть малы
    EXPECT_LT(inlet_variation, 1.0); // Менее 1 Па в среднем
    EXPECT_LT(outlet_variation, 1.0);
}

/// @brief Проверяет способность системы производить нулевой перепад давления при нулевом расходе
TEST_F(CondensatePipeQP, ProducesZeroPressureDrop_WhenFlowRateIsZero) {
    // Очень маленький расход
    {
        //Arrange
        pde_solvers::condensate_pipe_PQ_task_t task(pipe);
        auto conditions = pde_solvers::condensate_pipe_PQ_task_boundaries_t::default_values();
        conditions.volumetric_flow = 0; // нулевой расход

        EXPECT_NO_THROW(task.solve(conditions));

        auto& layer = task.get_current_layer();
        // При нулевом расходе перепад давления должен быть минимальным
        double pressure_drop = layer.pressure.front() - layer.pressure.back();
        EXPECT_EQ(pressure_drop, 0); // 
    }
}


class CondensatePPTest : public ::testing::Test {
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
TEST_F(CondensatePPTest, CalculatesValidFlowRate_ForGivenPressureBoundaries) {
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
TEST_F(CondensatePPTest, RecoversOriginalFlowRate_WhenPQPressureUsedInPP) {
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
TEST_F(CondensatePPTest, MaintainsFlowRateStability_OverMultipleTimeSteps) {
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
TEST_F(CondensatePPTest, DecreasesFlowRate_WithIncreasingDensity_AtFixedPressureDrop) {
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
TEST_F(CondensatePPTest, ProducesZeroFlowRate_WhenPressureDropIsZero) {
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
