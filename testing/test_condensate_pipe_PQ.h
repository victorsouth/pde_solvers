#pragma once

/// @brief Создает простую равномерную трубу для тестовых расчетов
inline pde_solvers::iso_nonbaro_pipe_properties_t create_test_pipe_for_PQ() {
    pde_solvers::iso_nonbaro_pipe_properties_t pipe;
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
    pipe.wall.wallThickness = 10e-3; // 10 мм
    pipe.kinematic_viscosity = 1e-6; // 1 сСт

    return pipe;
}

/// @brief Возвращает знаки приращений функции, встертившиеся 
/// при увеличении аргумента от start_arg до end_arg с шагом step
/// @tparam OneArgumentFunction Функция double(double)
template <typename OneArgumentFunction>
inline std::set<int> get_increments_sign_for_func(OneArgumentFunction f, 
    double start_arg, double end_arg, double step = 1e-3) 
{
    if (start_arg > end_arg)
        throw std::runtime_error("Wrong start and end values of argument");
    
    std::set<int> increment_sign;
    bool has_previous = false;
    double previous_func = 0.0;
    for (double x = start_arg; x <= end_arg; x += step)
    {
        double new_func = f(x);

        if (!has_previous) {
            previous_func = new_func;
            has_previous = true;
            continue;
        }

        int sgn = fixed_solvers::sgn(new_func - previous_func);
        increment_sign.insert(sgn);
        
        previous_func = new_func;
    }

    return increment_sign;
}

/// @brief Проверяет способность системы PQ задачи поддерживать корректный профиль давления при заданных начальных условиях.
/// Тест выполняет следующие проверки:
/// 1. Корректность инициализации плотности - все значения плотности в профиле должны соответствовать заданному начальному значению (800.0 кг/м³) с точностью 1e-6
/// 2. Наличие рассчитанного профиля давления - профиль давления не должен быть пустым
/// 3. Сохранение граничного условия на входе - давление на входе должно точно совпадать с заданным значением (5 МПа) с точностью 1e-6
/// 4. Монотонное убывание давления вдоль трубы - давление должно строго уменьшаться от входа к выходу из-за гидравлических потерь
/// 5. Физическая корректность значений - все значения давления в профиле должны быть положительными
TEST(CondensatePipeQPTask, MaintainsCorrectPressureProfile_WhenGivenInitialConditions) {

    //Arrange
    auto pipe = create_test_pipe_for_PQ();
    pde_solvers::iso_nonbarotropic_pipe_PQ_task_t task(pipe);

    auto initial_conditions = pde_solvers::iso_nonbarotropic_pipe_PQ_task_boundaries_t::default_values();
    initial_conditions.volumetric_flow = 0.3; // м³/с
    initial_conditions.pressure_in = 5e6; // 5 МПа
    initial_conditions.density = 800.0; // кг/м³

    // Act
    task.solve(initial_conditions);
    auto& layer = task.get_current_layer();

    //Assert

    // Проверяем, что в профиль попала заданная плотность
    for (const auto& density : layer.density_std.value) {
        EXPECT_NEAR(density, initial_conditions.density, 1e-6);
    }

    // Давление на входе должно быть равно заданному
    ASSERT_FALSE(layer.pressure.empty());
    // Проверяем, что в профиль попало начальное давление
    EXPECT_NEAR(layer.pressure.front(), initial_conditions.pressure_in, 1e-6);

    // Давление должно уменьшаться вдоль трубы из-за гидравлических потерь
    for (size_t i = 1; i < layer.pressure.size(); ++i) {
        EXPECT_LT(layer.pressure[i], layer.pressure[i - 1]);
    }

    // Проверяем, что давление положительное
    for (const auto& pressure : layer.pressure) {
        EXPECT_GT(pressure, 0.0);
    }

    // Расход должен быть равен заданному
    ASSERT_NEAR(layer.std_volumetric_flow, initial_conditions.volumetric_flow, 1e-6);
}

    

/// @brief Проверяет способность системы поддерживать стабильность давления при множественных шагах по времени
TEST(CondensatePipeQPTask, MaintainsPressureStability_OverMultipleTimeSteps) {
    //Arrange
    auto pipe = create_test_pipe_for_PQ();
    pde_solvers::iso_nonbarotropic_pipe_PQ_task_t task(pipe);

    auto initial_conditions = pde_solvers::iso_nonbarotropic_pipe_PQ_task_boundaries_t::default_values();
    task.solve(initial_conditions);

    // Выполняем много шагов и проверяем стабильность
    const int steps = 100;
    double dt = 1.0; // 1 секунда

    std::vector<double> inlet_pressures;
    std::vector<double> outlet_pressures;
    std::vector<double> volumes;

    for (int i = 0; i < steps; ++i) {
        //Act
        task.step(dt, initial_conditions);
        auto& layer = task.get_current_layer();

        inlet_pressures.push_back(layer.pressure.front());
        outlet_pressures.push_back(layer.pressure.back());
        volumes.push_back(layer.std_volumetric_flow);
    }

    //Assert
    // При постоянных граничных условиях давления должны стабилизироваться

    auto [min_it, max_it] = std::minmax_element(inlet_pressures.begin(), inlet_pressures.end());
    double diff = std::abs(*max_it - *min_it);
    EXPECT_LT(diff, 1.0);

    auto [min_it1, max_it1] = std::minmax_element(outlet_pressures.begin(), outlet_pressures.end());
    double diff1 = std::abs(*max_it1 - *min_it1);
    EXPECT_LT(diff1, 1.0);

    // проверка соответсвия расхода заданному в условии
    for (int i = 0; i < steps; ++i) {
        ASSERT_NEAR(volumes[i], initial_conditions.volumetric_flow, 1e-6);
    }

}

/// @brief Проверяет монотонность замыкающего соотношения для трубы
/// При использовании для расчета гидравлического сопротивления 
/// при переходе между формулами Исаева и Шифринсона наблюдался скачок, нарушавший монотонность
TEST(CondensatePipeQPTask, CheckClosingRelationMonotonicity) {

    using hydro_solver = pde_solvers::iso_nonbaro_pipe_solver_t;

    // Arrange - создание солвера для трубы с известными параметрами
    double length = 500e3;
    double dx = 200;
    double diameter = 0.8;

    pde_solvers::pipe_json_data pipe_data{ diameter, 0, length };
    auto pipe_parameters = iso_nonbaro_pipe_properties_t(pipe_data);
    pipe_parameters.make_uniform_profile(dx);
 
    hydro_solver::layer_type layer(pipe_parameters.profile.get_point_count());

    double Pin = 872947.27809030749;
    double Pout = 100000.0;
    pde_solvers::rigorous_impulse_solver_PP<pde_solvers::iso_nonbaro_impulse_equation_t>
        pp_solver(pipe_parameters, layer, Pin, Pout);

    // Act - расчет знака приращения невязки в диапазоне расхода
    auto residuals = [&pp_solver](double flow) { return pp_solver.residuals(flow); };
    double start_flow = -0.03;
    double end_flow = 1.0;
    std::set<int> residuals_increment_sign = get_increments_sign_for_func(residuals, start_flow, end_flow);

    // Assert - Замыкающее соотношение монотонное - все приращения имеют одинаковый знак
    ASSERT_TRUE(residuals_increment_sign.size() == 1);
}


