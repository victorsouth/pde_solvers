#pragma once

/// @brief Класс для тестов advection_moc_solver
class AdvectionMocSolver : public ::testing::Test
{
public:
    /// @brief Подготовка модели трубопровода 
    static pipe_properties_t PrepareTestPipe()
    {
        // Зададимся исходными данными
        pipe_properties_t pipe;
        double x0 = 0;
        double xl = 1e5;
        double d = 1;
        size_t n = 100;

        pipe.profile = PipeProfile::create(n, x0, xl, 0, 0, 10e6);
        pipe.wall.diameter = d;

        return pipe;
    }

protected:
    // Модель трубопровода
    pipe_properties_t pipe;
    // Исходные данные
    double rho_init = 850;
    double rho_left = 870;
    double rho_right = 830;
    double T = 1000;

    virtual void SetUp() override {
        pipe = PrepareTestPipe();
    };
};

/// @brief Пример квазистационарного расчёта трубопровода для плотности
TEST_F(AdvectionMocSolver, UseCaseDensity)
{
    // Зададимся расходом внутри трубопровода
    double volumetric_flow = 0.5;

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
TEST_F(AdvectionMocSolver, ConsiderFlowSwap)
{
    // Зададимся расходом внутри трубопровода для инверсии потока
    double volumetric_flow = -0.5;

    // Буфер для хранения слоёв
    ring_buffer_t<vector<double>> buffer(2, pipe.profile.getPointCount());
    buffer.previous() = vector<double>(pipe.profile.getPointCount(), rho_init);
    
    advection_moc_solver solver(pipe, volumetric_flow, buffer.previous(), buffer.current());

    solver.step(rho_left, rho_right);

    // Проверим, что давление изменилось в конце трубопровода, но не изменилось в начале
    ASSERT_NEAR(buffer.current().back(), rho_right, 0.05);
    ASSERT_NEAR(buffer.current().front(), rho_init, 0.05);
}