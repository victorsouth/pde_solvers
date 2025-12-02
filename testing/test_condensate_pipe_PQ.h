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

}

