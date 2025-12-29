#pragma once

    /// @brief Создает простую равномерную трубу для тестовых расчетов
    inline pde_solvers::condensate_pipe_properties_t create_test_pipe_for_PP() {
        pde_solvers::condensate_pipe_properties_t pipe;
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

        return pipe;
    }


    /// @brief Проверяет способность системы PP задачи рассчитывать валидный расход для заданных граничных условий по давлению
    TEST(CondensatePipePPTask, CalculatesValidFlowRate_ForGivenPressureBoundaries) {
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
            EXPECT_NEAR(density, initial_conditions.density, 1e-6);
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


    /// @brief Проверяет способность системы восстанавливать исходный расход при использовании давления из PQ задачи в PP задаче
    TEST(CondensatePipePPTask, RecoversOriginalFlowRate_WhenPQPressureUsedInPP) {
        // Arrange
        auto pipe = create_test_pipe_for_PP();
        pde_solvers::condensate_pipe_PQ_task_t pq_task(pipe);
        pde_solvers::condensate_pipe_PP_task_t pp_task(pipe);

        auto pq_conditions = pde_solvers::condensate_pipe_PQ_task_boundaries_t::default_values();
        pq_conditions.pressure_in = 5e6;
        pq_conditions.volumetric_flow = 0.3;
        pq_conditions.density = 850.0;

        // Act
        // Рассчитываем PQ задачу
        pq_task.solve(pq_conditions);
        double calculated_pressure_out = pq_task.get_current_layer().pressure.back();

        // Используем найденное давление в PP задаче
        auto pp_conditions = pde_solvers::condensate_pipe_PP_task_boundaries_t::default_values();
        pp_conditions.pressure_in = pq_conditions.pressure_in;
        pp_conditions.pressure_out = calculated_pressure_out;
        pp_conditions.density = pq_conditions.density;

        pp_task.solve(pp_conditions);

        // Assert: 
        // Расход должен совпадать
        double pp_flow = pp_task.get_current_layer().std_volumetric_flow;
        EXPECT_NEAR(pp_flow, pq_conditions.volumetric_flow, 1e-5);
    }

    /// @brief Проверяет способность системы PP задачи поддерживать стабильность расхода при множественных шагах по времени
    TEST(CondensatePipePPTask, MaintainsFlowRateStability_OverMultipleTimeSteps) {
        // Arrange
        auto pipe = create_test_pipe_for_PP();
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

        //Много шагов по времени
        for (int i = 0; i < steps; ++i) {
            task.step(dt, initial_conditions);
            flow_history.push_back(task.get_current_layer().std_volumetric_flow);
        }

        // Assert
        //проверяем, что расходы будут физичными
        EXPECT_GT(flow_history[0], 0.0);
        
        // Проверяем, что расход стабилизировался
        auto [min_it, max_it] = std::minmax_element(flow_history.begin(), flow_history.end());
        double diff = std::abs(*max_it - *min_it);
        EXPECT_LT(diff, 0.01); // погрешность по расходу меньше 0.01 м3/с
    }
