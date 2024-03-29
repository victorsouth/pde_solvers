#pragma once

/// @brief Пример квазистационарного расчёта трубопровода для плотности
TEST(TransportMocSolver, UseCaseDensity)
{
    // Зададимся исходными данными
    pipe_properties_t pipe;
    double x0 = 0;
    double xl = 1e5;
    double d = 1;
    size_t n = 100;
    double rho_init = 850;
    double rho_left = 870;
    double rho_right = 830;
    double volumetric_flow = 0.5;
    double T = 1000;

    // Модель трубопровода
    pipe.profile = PipeProfile::create(n, x0, xl, 0, 0, 10e6);
    pipe.wall.diameter = d;

    // Буфер для хранения слоёв
    ring_buffer_t<vector<double>> buffer(2, pipe.profile.getPointCount());
    buffer.previous() = vector<double>(pipe.profile.getPointCount(), rho_init);

    double modeling_time = 0;


    while (modeling_time <= T)
    {
        // Расход в трубопроводе может меняться во времени
        advection_moc_solver solver(pipe, volumetric_flow, buffer.previous(), buffer.current());
        // При изменение расхода меняется шаг по времени, так как Курант всегда должен равняться единице
        modeling_time += solver.prepare_step();

        solver.step(rho_left, rho_right);

        buffer.advance(+1);
    }

}

/// @brief Проверка на то, правильно ли учитывается инверсия потока 
TEST(TransportMocSolver, CanConsiderFlowSwap)
{
    // Зададимся исходными данными
    pipe_properties_t pipe;
    double x0 = 0;
    double xl = 1e5;
    double d = 1;
    size_t n = 100;
    double rho_init = 850;
    double rho_left = 870;
    double rho_right = 830;
    double volumetric_flow = -0.5;

    // Буфер для хранения слоёв
    pipe.profile = PipeProfile::create(n, x0, xl, 0, 0, 10e6);
    pipe.wall.diameter = d;

    // Модель трубопровода
    ring_buffer_t<vector<double>> buffer(2, pipe.profile.getPointCount());
    buffer.previous() = vector<double>(pipe.profile.getPointCount(), rho_init);
    

    advection_moc_solver solver(pipe, volumetric_flow, buffer.previous(), buffer.current());

    solver.step(rho_left, rho_right);

    // Проверим, что давление изменилось в конце трубопровода, но не изменилось в начале
    ASSERT_NEAR(buffer.current().back(), rho_right, 0.05);
    ASSERT_NEAR(buffer.current().front(), rho_init, 0.05);
}