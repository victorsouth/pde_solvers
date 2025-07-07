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

        pipe.profile = pipe_profile_t::create(n, x0, xl, 0, 0, 10e6);
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

    size_t point_count = pipe.profile.get_point_count();

    // Буфер для хранения слоёв
    ring_buffer_t<std::vector<double>> buffer(2, point_count);
    buffer.previous() = std::vector<double>(point_count, rho_init);

    double modeling_time = 0;

    while (modeling_time <= T)
    {
        // Расход в трубопроводе может меняться во времени
        advection_moc_solver solver(pipe, volumetric_flow, buffer.previous(), buffer.current());
        // При изменение расхода меняется шаг по времени, так как Курант всегда должен равняться единице
        double time_step = solver.prepare_step();
        modeling_time += time_step;

        solver.step(time_step, rho_left, rho_right);

        buffer.advance(+1);
    }

}

/// @brief Проверка на то, правильно ли учитывается инверсия потока 
TEST_F(AdvectionMocSolver, ConsiderFlowSwap)
{
    // Зададимся расходом внутри трубопровода для инверсии потока
    double volumetric_flow = -0.5;

    size_t point_count = pipe.profile.get_point_count();

    // Буфер для хранения слоёв
    ring_buffer_t<std::vector<double>> buffer(2, point_count);
    buffer.previous() = std::vector<double>(point_count, rho_init);
    
    advection_moc_solver solver(pipe, volumetric_flow, buffer.previous(), buffer.current());
    double time_step = solver.prepare_step();

    solver.step(time_step, rho_left, rho_right);

    // Проверим, что давление изменилось в конце трубопровода, но не изменилось в начале
    ASSERT_NEAR(buffer.current().back(), rho_right, 0.05);
    ASSERT_NEAR(buffer.current().front(), rho_init, 0.05);
}

/// @brief Проверка случая, когда скорость потока в ходе моделирования меняется и Курант становится меньше единицы
TEST_F(AdvectionMocSolver, ConsiderCrLessOne)
{
    // Зададимся максимальным расходом
    double volumetric_flow = 0.5;
    // Посчитаем шаг по времени, при котором Курант равен единице
    double time_step = (pipe.profile.coordinates[1] - pipe.profile.coordinates[0]) / (volumetric_flow / pipe.wall.getArea());

    size_t point_count = pipe.profile.get_point_count();

    // Буфер для хранения слоёв
    ring_buffer_t<std::vector<double>> buffer(2, point_count);
    buffer.current() = std::vector<double>(point_count, rho_init);

    // Проведём моделирование, в котором на первой итерации расход будет равен начальному
    // а на второй - расход станет вдвое меньше
    // Шаг по времени при этом остаётся таким же, поэтому Курант станет равен 0.5
    for (int i = 1; i <= 2; i++)
    {
        buffer.advance(+1);
        advection_moc_solver solver(pipe, volumetric_flow / i, buffer.previous(), buffer.current());
        solver.step(time_step, rho_left, rho_right);
    }

    // В такой ситуации после второй итерации значение плотности во второй точке текущего профиля
    // станет равна значению из середины первой и второй точки предыдущего профиля
    ASSERT_NEAR(buffer.current()[1], (buffer.previous()[0] + buffer.previous()[1]) / 2, 0.05);
}

/// @brief Проверка случая, когда скорость потока в ходе моделирования меняется и Курант становится меньше единицы при инверсном потоке
TEST_F(AdvectionMocSolver, ConsiderCrLessOneInverseFlow)
{
    // Зададимся максимальным расходом в инверсном направлении
    double volumetric_flow = -0.5;
    // Посчитаем шаг по времени, при котором Курант равен единице
    double time_step = (pipe.profile.coordinates[1] - pipe.profile.coordinates[0]) / (-volumetric_flow / pipe.wall.getArea());

    size_t point_count = pipe.profile.get_point_count();

    // Буфер для хранения слоёв
    ring_buffer_t<std::vector<double>> buffer(2, point_count);
    buffer.current() = std::vector<double>(point_count, rho_init);

    // Проведём моделирование, в котором на первой итерации расход будет равен начальному
    // а на второй - расход станет вдвое меньше
    // Шаг по времени и по координате при этом остаётся таким же, поэтому Курант станет равен 0.5
    for (int i = 1; i <= 2; i++)
    {
        buffer.advance(+1);
        advection_moc_solver solver(pipe, volumetric_flow / i, buffer.previous(), buffer.current());
        solver.step(time_step, rho_left, rho_right);
    }

    // В такой ситуации после второй итерации значение плотности в предпоследней точке текущего профиля
    // станет равна значению из середины последней и предпоследней точки предыдущего профиля
    ASSERT_NEAR(buffer.current()[point_count - 2], (buffer.previous()[point_count - 1] + buffer.previous()[point_count - 2]) / 2, 0.05);
}