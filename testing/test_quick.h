#pragma once

/// @brief Учитывает инверсию потока: при отрицательном расходе плотность изменяется с конца трубы
TEST(UpstreamDifferencing, ConsidersFlowSwap) {
    using namespace upstream_solver_types;

    // Arrange: труба, модель адвекции и буфер
    pipe_properties_t pipe = pipe_properties_t::build_simple_pipe(
        simple_pipe_properties::sample_district());
    double rho_initial = 850;
    double rho_in = 860;
    double rho_out = 870;
    double dt = 60;

    std::vector<double> Q(pipe.profile.get_point_count(), -0.5);
    PipeQAdvection advection_model(pipe, Q);
    ring_buffer_t<layer_t> buffer = build_buffer<layer_t>(pipe, rho_initial);

    // Act: один шаг адвекции
    upstream_fv_solver solver(advection_model, buffer.previous(), buffer.current());
    solver.step(dt, rho_in, rho_out);

    // Assert: плотность в конце выросла, в начале не изменилась
    const auto& rho_prev = buffer.previous().vars.cell_double[0];
    const auto& rho_curr = buffer.current().vars.cell_double[0];
    ASSERT_GT(rho_curr.back(), rho_prev.back());
    ASSERT_NEAR(rho_curr.front(), rho_prev.front(), 1e-8);
}

/// @brief Дымовой тест метода верхней разности
TEST(UpstreamDifferencing, ExecutesSingleStep)
{
    using namespace upstream_solver_types;

    // Arrange: труба, модель адвекции и буфер
    pipe_properties_t pipe = pipe_properties_t::build_simple_pipe(simple_pipe_properties::sample_district());
    double rho_initial = 850;
    double rho_in = 860;
    double rho_out = 870;
    double dt = 60;

    std::vector<double> Q(pipe.profile.get_point_count(), 0.5);
    PipeQAdvection advection_model(pipe, Q);
    ring_buffer_t<layer_t> buffer = build_buffer<layer_t>(pipe, rho_initial);

    // Act: один шаг адвекции
    upstream_fv_solver solver(advection_model, buffer.previous(), buffer.current());
    solver.step(dt, rho_in, rho_out);
}

/// @brief Солвер Quickest учитывает инверсию потока:
/// при перекачке в обратном направлении более тяжёлая нефть приходит с конца трубы
TEST(QuickestUltimate, ConsidersFlowSwap) {
    using namespace quickest_ultimate_solver_types;

    // Arrange: труба, модель адвекции и буфер
    pipe_properties_t pipe = pipe_properties_t::build_simple_pipe({.length = 700e3, .dx = 100, .diameter = 0.514});
    double rho_initial = 850;
    double rho_in = 860;
    double rho_out = 870;

    std::vector<double> Q(pipe.profile.get_point_count(), -0.5);
    PipeQAdvection advection_model(pipe, Q);
    ring_buffer_t<layer_t> buffer = build_buffer<layer_t>(pipe, rho_initial);
    double dt = calc_time_step_by_Courant(advection_model, 0.9);

    // Act: один шаг адвекции
    quickest_sequential_solver solver(advection_model, buffer.previous(), buffer.current());
    solver.step(dt, rho_in, rho_out);

    // Assert: плотность в конце выросла, в начале не изменилась
    const auto& rho_prev = buffer.previous().vars.cell_double[0];
    const auto& rho_curr = buffer.current().vars.cell_double[0];
    ASSERT_GT(rho_curr.back(), rho_prev.back());
    ASSERT_NEAR(rho_curr.front(), rho_prev.front(), 1e-8);
}

/// @brief Дымовой тест метода Quickest Ultimate на тепловой модели
TEST(QuickestUltimate, ExecutesTemperatureAdvection)
{
    using namespace quickest_ultimate_solver_types;

    // Arrange: труба, нефть, тепловая модель и буфер
    oil_parameters_t oil = pde_solvers::get_default_oil_heatmodel();
    pipe_noniso_properties_t pipe = pde_solvers::get_default_pipe_heatmodel();
    pipe.heat.ambient_heat_transfer = 3;

    double temp_initial = 300;
    double temp_in = 310;
    double temp_out = 290;
    double dt = 1000;
    double mass_flow = 300;

    std::vector<double> G(pipe.profile.get_point_count(), mass_flow);
    PipeHeatInflowConstArea heat_model(pipe, oil, G);
    ring_buffer_t<layer_t> buffer = build_buffer<layer_t>(pipe, temp_initial);

    // Act: 10 шагов адвекции температуры
    for (size_t index = 0; index < 10; ++index) {
        quickest_sequential_solver solver(heat_model, buffer);
        solver.step(dt, temp_in, temp_out);
        buffer.advance(+1);
    }

    // TODO: Assert
}

/// @brief Проверка возможности применения метода Quickest Ultimate для трубы из одной ячейки
/// При прямом потоке в ячейке мгновенно устанавливается значение из краевого условия на входе.
TEST(QuickestUltimate, HandlesShortPipe_Direct)
{
    using namespace quickest_ultimate_solver_types;

    // Arrange: труба из одной ячейки (length/dx → 1 сегмент), модель адвекции и буфер
    pipe_properties_t pipe = pipe_properties_t::build_simple_pipe({.length = 50.0, .dx = 50.0, .diameter = 0.7});
    double rho_initial = 840;
    double rho_in = 860;
    double rho_out = 850;
    double flow = +1.0;
    double Cr = 1.0;

    std::vector<double> Q(pipe.profile.get_point_count(), flow);
    PipeQAdvection advection_model(pipe, Q);
    ring_buffer_t<layer_t> buffer = build_buffer<layer_t>(pipe, rho_initial);
    double dt = calc_time_step_by_Courant(advection_model, Cr);

    // Act: один шаг адвекции
    quickest_sequential_solver solver(advection_model, buffer.previous(), buffer.current());
    solver.step(dt, rho_in, rho_out);

    // Assert: плотность в ячейке стремится к граничному значению на входе
    double rho_prev = buffer.previous().vars.cell_double[0].front();
    double rho_curr = buffer.current().vars.cell_double[0].front();
    ASSERT_NEAR(rho_curr, rho_in, 1e-6);
    ASSERT_NE(rho_curr, rho_prev);
}

/// @brief Проверка возможности применения метода Quickest Ultimate для трубы из одной ячейки
/// При обратном потоке в ячейке мгновенно устанавливается значение из краевого условия на выходе.
TEST(QuickestUltimate, HandlesShortPipe_Reverse)
{
    using namespace quickest_ultimate_solver_types;

    // Arrange: труба из одной ячейки, модель адвекции и буфер
    pipe_properties_t pipe = pipe_properties_t::build_simple_pipe({.length = 50.0, .dx = 50.0, .diameter = 0.7});
    double rho_initial = 840;
    double rho_in = 860;
    double rho_out = 850;
    double flow = -1.0;

    std::vector<double> Q(pipe.profile.get_point_count(), flow);
    PipeQAdvection advection_model(pipe, Q);
    ring_buffer_t<layer_t> buffer = build_buffer<layer_t>(pipe, rho_initial);

    double Cr = 1.0;
    double dt = calc_time_step_by_Courant(advection_model, Cr);

    // Act: один шаг адвекции
    quickest_sequential_solver solver(advection_model, buffer.previous(), buffer.current());
    solver.step(dt, rho_in, rho_out);

    // Assert: плотность в ячейке стремится к граничному значению на выходе
    double rho_prev = buffer.previous().vars.cell_double[0].front();
    double rho_curr = buffer.current().vars.cell_double[0].front();
    ASSERT_NEAR(rho_curr, rho_out, 1e-6);
    ASSERT_NE(rho_curr, rho_prev);
}

/// @brief Результаты параллельного и последовательного расчетов совпадают.
/// В начальном состоянии параметр в каждой ячейке различен.
/// Вероятно, если при одном шаге расчета результаты совпадают, то и при длительном расчете они будут совпадать.
TEST(MultiThreadedQuickestUltimate, IsCoherentWithSingleThreaded)
{
    // Arrange: Подготовка трубы с большим количеством ячеек и линейным профилем плотности
    pipe_properties_t pipe = pipe_properties_t::build_simple_pipe(
        {.length = 1'000'000.0, .dx = 1.0, .diameter = 0.7});

    double rho_initial_left = 850;
    double rho_initial_right = 860;
    double rho_new = 870;
    double flow = 0.5;

    std::vector<double> Q(pipe.profile.get_point_count(), flow);
    PipeQAdvection advection_model(pipe, Q);
    double dt = calc_time_step_by_Courant(advection_model, 1.0);

    // Act: Расчет адвекции в последовательном и параллельном режимах
    // TODO: явно сгенерировать слой в виде вектора через STL и передать в конструктор
    auto seqential_calculation = build_linear_buffer<quickest_ultimate_solver_types::layer_t>(pipe, rho_initial_left, rho_initial_right);
    quickest_ultimate_solver_types::quickest_sequential_solver(advection_model, seqential_calculation.previous(), seqential_calculation.current())
        .step(dt, rho_new, rho_initial_right);
    std::vector<double> sequential_profile = seqential_calculation.current().vars.cell_double[0];

    auto parallel_calculation = build_linear_buffer<quickest_ultimate_solver_types::layer_t>(pipe, rho_initial_left, rho_initial_right);
    quickest_ultimate_solver_types::quickest_parallel_solver(advection_model, parallel_calculation.previous(), parallel_calculation.current())
        .step(dt, rho_new, rho_initial_right);
    std::vector<double> parallel_profile = parallel_calculation.current().vars.cell_double[0];

    // Assert: Расчетные профили параметра в каждой ячейке совпадают
    ASSERT_EQ(sequential_profile, parallel_profile); // GTest корректно сравнивает векторы
}

/// @brief На длинной трубе параллельный расчет дает выигрыш во времени по сравнению с последовательным.
/// Ожидаемое ускорение пропорционально числу аппаратных потоков.
TEST(MultiThreadedQuickestUltimate, IncreasesPerformance)
{
    using namespace quickest_ultimate_solver_types;

    unsigned int threads = std::thread::hardware_concurrency();
    if (threads <= 2) {
        // При двух ядрах случайный разбег времени выполнения может быть близок к ожидаемому ускорению 
        GTEST_SKIP() << "для проверки ускорения нужно не менее 2 аппаратных потоков";
    }
    
    // Arrange: Подготовка трубы с большим количеством ячеек
    pipe_properties_t pipe = pipe_properties_t::build_simple_pipe(
        {.length = 500, .dx = 1.0, .diameter = 0.7});

    double rho_initial_left = 850;
    double rho_initial_right = 860;
    double rho_new = 870;
    double flow = 0.5;

    std::vector<double> Q(pipe.profile.get_point_count(), flow);
    PipeQAdvection advection_model(pipe, Q);
    double dt = calc_time_step_by_Courant(advection_model, 0.5);

    constexpr size_t step_count = 100;

    // Act: Расчет адвекции в последовательном и параллельном режимах
    double seqential_calculation_time = run_timed_step<quickest_sequential_solver, layer_t>(
        advection_model, pipe, dt, rho_initial_left, rho_initial_right, rho_new, step_count).second;
    double parallel_calculation_time = run_timed_step<quickest_parallel_solver, layer_t>(
        advection_model, pipe, dt, rho_initial_left, rho_initial_right, rho_new, step_count).second;

    // Assert: Фактическое ускорение не меньше ожидаемого
    double actual_speedup = seqential_calculation_time / parallel_calculation_time;
    double expected_speedup = threads / 2.0;
    EXPECT_GE(actual_speedup, expected_speedup) << "Потоков: " << threads << ", ожидаемое ускорение: " 
        << expected_speedup << "x, фактическое: " << actual_speedup << "x";
}
