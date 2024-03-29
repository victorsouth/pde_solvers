#pragma once

/// @brief Пример квазистационарного расчёта трубопровода
TEST(TransportMocSolver, UseCase)
{
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

    pipe.profile = PipeProfile::create(n, x0, xl, 0, 0, 10e6);
    pipe.wall.diameter = d;

    ring_buffer_t<vector<double>> buffer(2, pipe.profile.getPointCount());
    buffer.previous() = vector<double>(pipe.profile.getPointCount(), rho_init);

    double modeling_time = 0;


    while (modeling_time <= T)
    {
        transport_moc_solver solver(pipe, volumetric_flow, buffer.previous(), buffer.current());
        modeling_time += solver.prepare_step();

        solver.step(rho_left, rho_right);

        buffer.advance(+1);
    }

}

/// @brief Проверка на то, правильно ли учитывается инверсия потока 
TEST(TransportMocSolver, CanConsiderFlowSwap)
{
    pipe_properties_t pipe;
    double x0 = 0;
    double xl = 1e5;
    double d = 1;
    size_t n = 100;
    double rho_init = 850;
    double rho_left = 870;
    double rho_right = 830;
    double volumetric_flow = -0.5;

    pipe.profile = PipeProfile::create(n, x0, xl, 0, 0, 10e6);
    pipe.wall.diameter = d;

    ring_buffer_t<vector<double>> buffer(2, pipe.profile.getPointCount());
    buffer.previous() = vector<double>(pipe.profile.getPointCount(), rho_init);

    transport_moc_solver solver(pipe, volumetric_flow, buffer.previous(), buffer.current());

    solver.step(rho_left, rho_right);

    ASSERT_NEAR(buffer.current().back(), rho_right, 0.05);
    ASSERT_NEAR(buffer.current().front(), rho_init, 0.05);
}