#pragma once






/// @brief Тесты для солвера upstream_fv_solver
class UpstreamDifferencing : public ::testing::Test {
protected:
    // Профиль переменных
    typedef upstream_fv_solver_traits<1>::var_layer_data target_var_t;
    typedef upstream_fv_solver_traits<1>::specific_layer specific_data_t;

    // Слой: переменных Vars + сколько угодно служебных Specific
    typedef composite_layer_t<target_var_t, specific_data_t> layer_t;
protected:
    /// @brief Параметры трубы
    pipe_properties_t pipe;
    /// @brief Профиль расхода
    std::vector<double> Q;
    std::unique_ptr<PipeQAdvection> advection_model;
    std::unique_ptr<ring_buffer_t<layer_t>> buffer;
protected:
    
    /// @brief Подготовка к расчету для семейства тестов
    virtual void SetUp() override {
        // Упрощенное задание трубы - 50км, с шагом разбиения для расчтной сетки 1км, диаметром 700мм
        simple_pipe_properties simple_pipe = simple_pipe_properties::sample_district();
        pipe = pipe_properties_t::build_simple_pipe(simple_pipe);

        Q = std::vector<double> (pipe.profile.get_point_count(), 0.5);
        advection_model = std::make_unique<PipeQAdvection>(pipe, Q);
        buffer = std::make_unique<ring_buffer_t<layer_t>>(2, pipe.profile.get_point_count());

        layer_t& prev = buffer->previous();
        prev.vars.cell_double[0] = std::vector<double>(prev.vars.cell_double[0].size(), 850);
    }
};

/// @brief Тесты для солвера quick_fv_solver
class QUICK : public ::testing::Test {
protected:
    // Профиль переменных
    typedef quick_fv_solver_traits<1>::var_layer_data target_var_t;
    typedef quick_fv_solver_traits<1>::specific_layer specific_data_t;

    // Слой: переменных Vars + сколько угодно служебных Specific
    typedef composite_layer_t<target_var_t, specific_data_t> layer_t;
protected:
    /// @brief Параметры трубы
    pipe_properties_t pipe;
    /// @brief Профиль расхода
    std::vector<double> Q;
    std::unique_ptr<PipeQAdvection> advection_model;
    std::unique_ptr<ring_buffer_t<layer_t>> buffer;
protected:

    /// @brief Подготовка к расчету для семейства тестов
    virtual void SetUp() override {
        // Упрощенное задание трубы - 50км, с шагом разбиения для расчтной сетки 1км, диаметром 700мм
        simple_pipe_properties simple_pipe;
        //simple_pipe.length = 50e3;
        simple_pipe.length = 700e3; // тест трубы 700км
        //simple_pipe.diameter = 0.7;
        simple_pipe.diameter = 0.514; // тест трубы 700км
        //simple_pipe.dx = 100;
        simple_pipe.dx = 100; // тест трубы 700км
        pipe = pipe_properties_t::build_simple_pipe(simple_pipe);

        Q = std::vector<double>(pipe.profile.get_point_count(), 0.5);
        advection_model = std::make_unique<PipeQAdvection>(pipe, Q);
        buffer = std::make_unique<ring_buffer_t<layer_t>>(2, pipe.profile.get_point_count());

        layer_t& prev = buffer->previous();
        prev.vars.cell_double[0] = std::vector<double>(prev.vars.cell_double[0].size(), 850);
    }
};

/// @brief Тесты для солвера quickest_fv_solver
class QUICKEST : public ::testing::Test {
protected:
    // Профиль переменных
    typedef quickest_fv_solver_traits<1>::var_layer_data target_var_t;
    typedef quickest_fv_solver_traits<1>::specific_layer specific_data_t;

    // Слой: переменных Vars + сколько угодно служебных Specific
    typedef composite_layer_t<target_var_t, specific_data_t> layer_t;
protected:
    /// @brief Параметры трубы
    pipe_properties_t pipe;
    /// @brief Профиль расхода
    std::vector<double> Q;
    std::unique_ptr<PipeQAdvection> advection_model;
    std::unique_ptr<ring_buffer_t<layer_t>> buffer;
protected:

    /// @brief Подготовка к расчету для семейства тестов
    virtual void SetUp() override {
        // Упрощенное задание трубы - 50км, с шагом разбиения для расчтной сетки 1км, диаметром 700мм
        simple_pipe_properties simple_pipe;
        //simple_pipe.length = 50e3;
        simple_pipe.length = 700e3; // тест трубы 700км
        //simple_pipe.diameter = 0.7;
        simple_pipe.diameter = 0.514; // тест трубы 700км
        //simple_pipe.dx = 100;
        simple_pipe.dx = 100; // тест трубы 700км
        pipe = pipe_properties_t::build_simple_pipe(simple_pipe);

        Q = std::vector<double>(pipe.profile.get_point_count(), 0.5);
        advection_model = std::make_unique<PipeQAdvection>(pipe, Q);
        buffer = std::make_unique<ring_buffer_t<layer_t>>(2, pipe.profile.get_point_count());

        layer_t& prev = buffer->previous();
        prev.vars.cell_double[0] = std::vector<double>(prev.vars.cell_double[0].size(), 850);
    }
};

/// @brief Тесты для солвера quickest_ultimate_fv_solver
class QUICKEST_ULTIMATE : public ::testing::Test {
protected:
    // Профиль переменных
    typedef quickest_ultimate_fv_solver_traits<1>::var_layer_data target_var_t;
    typedef quickest_ultimate_fv_solver_traits<1>::specific_layer specific_data_t;

    // Слой: переменных Vars + сколько угодно служебных Specific
    typedef composite_layer_t<target_var_t, specific_data_t> layer_t;
protected:
    /// @brief Параметры трубы
    pipe_properties_t pipe;
    /// @brief Профиль расхода
    std::vector<double> Q;
    std::unique_ptr<PipeQAdvection> advection_model;
    std::unique_ptr<ring_buffer_t<layer_t>> buffer;
protected:

    /// @brief Подготовка к расчету для семейства тестов
    virtual void SetUp() override {
        // Упрощенное задание трубы - 50км, с шагом разбиения для расчтной сетки 1км, диаметром 700мм
        simple_pipe_properties simple_pipe;
        //simple_pipe.length = 50e3;
        simple_pipe.length = 700e3; // тест трубы 700км
        //simple_pipe.diameter = 0.7;
        simple_pipe.diameter = 0.514; // тест трубы 700км
        //simple_pipe.dx = 100;
        simple_pipe.dx = 100; // тест трубы 700км
        pipe = pipe_properties_t::build_simple_pipe(simple_pipe);

        Q = std::vector<double>(pipe.profile.get_point_count(), 0.5);
        advection_model = std::make_unique<PipeQAdvection>(pipe, Q);
        buffer = std::make_unique<ring_buffer_t<layer_t>>(2, pipe.profile.get_point_count());

        layer_t& prev = buffer->previous();
        prev.vars.cell_double[0] = std::vector<double>(prev.vars.cell_double[0].size(), 850);
    }
};

/// @brief Проверка, правильно ли учитывается инверсия потока
TEST_F(UpstreamDifferencing, ConsidersFlowSwap) {
    Q = std::vector<double>(pipe.profile.get_point_count(), -0.5);

    //Получение текущего/предыдущего слоя
    layer_t& prev = buffer->previous();
    layer_t& next = buffer->current();

    double rho_in = 860;
    double rho_out = 870;
    double dt = 60; // 1 минута

    upstream_fv_solver solver(*advection_model, prev, next);
    solver.step(dt, rho_in, rho_out);

    const auto& rho_prev = prev.vars.cell_double[0];
    const auto& rho_curr = next.vars.cell_double[0];
    ASSERT_GT(rho_curr.back(), rho_prev.back()); // плотность в конце выросла
    ASSERT_NEAR(rho_curr.front(), rho_prev.front(), 1e-8); // плотность в начале не изменилась
}

/// @brief Разработка метода прямых разностей по [Leonard 1979]
TEST_F(UpstreamDifferencing, UseCaseSingleStep)
{
    //Получение текущего/предыдущего слоя
    layer_t& prev = buffer->previous();
    layer_t& next = buffer->current();

    double rho_in = 860;
    double rho_out = 870;
    double dt = 60; // 1 минута

    upstream_fv_solver solver(*advection_model, prev, next);
    solver.step(dt, rho_in, rho_out);

}

/// @brief Пример вывода в файл через
TEST_F(UpstreamDifferencing, DISABLED_UseCaseStepDensity)
{
    std::string path = prepare_test_folder();
    std::ofstream output(path + "output.csv");

    double rho_in = 860;
    double rho_out = 870;
    double T = 350000; // период моделирования

    const auto& x = advection_model->get_grid();
    double dx = x[1] - x[0];
    double v = advection_model->getEquationsCoeffs(0, 0);
    double dt_ideal = dx / v;

    for (double Cr = 0.5; Cr < 0.51; Cr += 0.05) {
        advection_model = std::make_unique<PipeQAdvection>(pipe, Q);
        buffer = std::make_unique<ring_buffer_t<layer_t>>(2, pipe.profile.get_point_count());

        layer_t& prev = buffer->previous();
        prev.vars.cell_double[0] = std::vector<double>(prev.vars.cell_double[0].size(), 850);

        double t = 0; // текущее время
        //double dt = 60; // 1 минута
        double dt = Cr * dt_ideal; // время в долях от Куранта

        std::stringstream filename;
        filename << path << "output Cr=" << Cr << ".csv";
        std::ofstream output(filename.str());


        size_t N = static_cast<int>(T / dt);
        for (size_t index = 0; index < N; ++index) {
            if (index == 0) {
                layer_t& prev = buffer->previous();
                prev.vars.print(t, output);
            }

            t += dt;

            upstream_fv_solver solver(*advection_model, *buffer);
            solver.step(dt, rho_in, rho_out);

            layer_t& next = buffer->current();
            next.vars.print(t, output);

            buffer->advance(+1);

        }
        output.flush();
        output.close();
    }
}

/// @brief Пример вывода в файл через
TEST_F(QUICK, DISABLED_UseCaseStepDensity)
{
    std::string path = prepare_test_folder();


    double rho_in = 860;
    double rho_out = 870;
    double T = 350000; // период моделирования

    const auto& x = advection_model->get_grid();
    double dx = x[1] - x[0];
    double v = advection_model->getEquationsCoeffs(0, 0);
    double dt_ideal = abs(dx / v);


    for (double Cr = 0.5; Cr < 0.51; Cr += 0.05) {
        advection_model = std::make_unique<PipeQAdvection>(pipe, Q);
        buffer = std::make_unique<ring_buffer_t<layer_t>>(2, pipe.profile.get_point_count());

        layer_t& prev = buffer->previous();
        prev.vars.cell_double[0] = std::vector<double>(prev.vars.cell_double[0].size(), 850);

        double t = 0; // текущее время
        //double dt = 60; // 1 минута
        double dt = Cr * dt_ideal; // время в долях от Куранта

        std::stringstream filename;
        filename << path << "output Cr=" << Cr << ".csv";
        std::ofstream output(filename.str());


        size_t N = static_cast<int>(T / dt);
        for (size_t index = 0; index < N; ++index) {
            if (index == 0) {
                layer_t& prev = buffer->previous();
                prev.vars.print(t, output);
            }

            t += dt;

            quick_fv_solver solver(*advection_model, *buffer);
            solver.step(dt, rho_in, rho_out);

            layer_t& next = buffer->current();
            next.vars.print(t, output);

            buffer->advance(+1);

        }
        output.flush();
        output.close();
    }

}

/// @brief Пример вывода в файл через
TEST_F(QUICKEST, DISABLED_UseCaseStepDensity)
{
    std::string path = prepare_test_folder();


    double rho_in = 860;
    double rho_out = 870;
    double T = 350000; // период моделирования

    const auto& x = advection_model->get_grid();
    double dx = x[1] - x[0];
    double v = advection_model->getEquationsCoeffs(0, 0);
    double dt_ideal = abs(dx / v);


    for (double Cr = 0.5; Cr < 0.51; Cr += 0.05) {
        advection_model = std::make_unique<PipeQAdvection>(pipe, Q);
        buffer = std::make_unique<ring_buffer_t<layer_t>>(2, pipe.profile.get_point_count());

        layer_t& prev = buffer->previous();
        prev.vars.cell_double[0] = std::vector<double>(prev.vars.cell_double[0].size(), 850);

        double t = 0; // текущее время
        //double dt = 60; // 1 минута
        double dt = Cr * dt_ideal; // время в долях от Куранта

        std::stringstream filename;
        filename << path << "output Cr=" << Cr << ".csv";
        std::ofstream output(filename.str());


        size_t N = static_cast<int>(T / dt);
        for (size_t index = 0; index < N; ++index) {
            if (index == 0) {
                layer_t& prev = buffer->previous();
                prev.vars.print(t, output);
            }

            t += dt;

            quickest_fv_solver solver(*advection_model, *buffer);
            solver.step(dt, rho_in, rho_out);

            layer_t& next = buffer->current();
            next.vars.print(t, output);

            buffer->advance(+1);

        }
        output.flush();
        output.close();
    }

}

//TEST_MOC
/// @brief Базовый пример использования метода характеристик для уравнения адвекции
TEST(MOC_Solver, DISABLED_MOC_Compare_With_QUICK)
{
    
    // Упрощенное задание трубы - 50км, с шагом разбиения для расчтной сетки 1км, диаметром 700мм
    simple_pipe_properties simple_pipe;
    //simple_pipe.length = 50e3;
    simple_pipe.length = 700e3; // тест трубы 700км
    //simple_pipe.diameter = 0.7;
    simple_pipe.diameter = 0.514; // тест трубы 700км
    //simple_pipe.dx = 100;
    simple_pipe.dx = 100; // тест трубы 700км

    std::string path = prepare_test_folder();

    double rho_in = 860;
    double rho_out = 870;
    double T = 350000; // период моделирования

    pipe_properties_t pipe = pipe_properties_t::build_simple_pipe(simple_pipe);

    // Одна переменная, и структуры метода характеристик для нееm
    typedef composite_layer_t<profile_collection_t<1>,
        moc_solver<1>::specific_layer> single_var_moc_t;


    std::vector<double> Q(pipe.profile.get_point_count(), 0.5); // задаем по трубе расход 0.5 м3/с
    PipeQAdvection advection_model(pipe, Q);

    const auto& x = advection_model.get_grid();
    double dx = x[1] - x[0];
    double v = advection_model.getEquationsCoeffs(0, 0);
    double dt_ideal = abs(dx / v);

    for (double Cr = 0.05; Cr < 0.51; Cr += 0.05) {
        PipeQAdvection advection_model(pipe, Q);
        ring_buffer_t<single_var_moc_t> buffer(2, pipe.profile.get_point_count());
        buffer.advance(+1);
        single_var_moc_t& prev = buffer.previous();
        single_var_moc_t& next = buffer.current();
        auto& rho_initial = prev.vars.point_double[0];
        rho_initial = std::vector<double>(rho_initial.size(), 850);
        double t = 0; // текущее время
        double dt = Cr * dt_ideal; // время в долях от Куранта
        std::stringstream filename;
        filename << path << "output Cr=" << Cr << ".csv";
        std::ofstream output(filename.str());
        size_t N = static_cast<int>(T / dt);
        for (size_t index = 0; index < N; ++index) {
            
            if (index == 0) {
                single_var_moc_t& prev = buffer.previous();
            }

            t += dt;
            moc_solver<1> solver(advection_model, buffer.previous(), buffer.current());
            solver.step_optional_boundaries(dt, rho_in, rho_out);

            single_var_moc_t& next = buffer.current();
            next.vars.print(t, output);

            buffer.advance(+1);

        }
        output.flush();
        output.close();
    }
}

/// @brief Пример вывода в файл через
TEST_F(QUICKEST_ULTIMATE, DISABLED_UseCaseStepDensity)
{
    std::string path = prepare_test_folder();


    double rho_in = 860;
    double rho_out = 870;
    double T = 350000; // период моделирования
    //double T = 800000; // период моделирования (тест трубы 700км)

    const auto& x = advection_model->get_grid();
    double dx = x[1] - x[0];
    double v = advection_model->getEquationsCoeffs(0, 0);
    double dt_ideal = abs(dx / v);


    for (double Cr = 0.05; Cr < 1.01; Cr += 0.05) {
        advection_model = std::make_unique<PipeQAdvection>(pipe, Q);
        buffer = std::make_unique<ring_buffer_t<layer_t>>(2, pipe.profile.get_point_count());

        layer_t& prev = buffer->previous();
        prev.vars.cell_double[0] = std::vector<double>(prev.vars.cell_double[0].size(), 850);

        double t = 0; // текущее время
        //double dt = 60; // 1 минута
        double dt = Cr * dt_ideal; // время в долях от Куранта

        std::stringstream filename;
        filename << path << "output Cr=" << Cr << ".csv";
        std::ofstream output(filename.str());


        size_t N = static_cast<int>(T / dt);
        for (size_t index = 0; index < N; ++index) {
            if (index == 0) {
                layer_t& prev = buffer->previous();
                prev.vars.print(t, output);
            }

            t += dt;

            quickest_ultimate_fv_solver solver(*advection_model, *buffer);
            solver.step(dt, rho_in, rho_out);

            layer_t& next = buffer->current();
            next.vars.print(t, output);


            buffer->advance(+1);

        }
        output.flush();
        output.close();
    }

}

/// @brief Проверка QUICKEST-ULTIMATE, правильно ли учитывается инверсия потока
TEST_F(QUICKEST_ULTIMATE, ConsidersFlowSwap) {
    Q = std::vector<double>(pipe.profile.get_point_count(), -0.5);

    //Получение текущего/предыдущего слоя
    layer_t& prev = buffer->previous();
    layer_t& next = buffer->current();

    double rho_in = 860;
    double rho_out = 870;

    // Рассчитываем dt для Cr = 0.9, чтобы выполнялось условие Куранта
    const auto& x = advection_model->get_grid();
    double dx = x[1] - x[0];
    double v = advection_model->getEquationsCoeffs(0, 0);
    double dt = 0.9 * dx / abs(v);

    quickest_ultimate_fv_solver solver(*advection_model, prev, next);
    solver.step(dt, rho_in, rho_out);

    const auto& rho_prev = prev.vars.cell_double[0];
    const auto& rho_curr = next.vars.cell_double[0];
    ASSERT_GT(rho_curr.back(), rho_prev.back()); // плотность в конце выросла
    ASSERT_NEAR(rho_curr.front(), rho_prev.front(), 1e-8); // плотность в начале не изменилась
}



/// @brief Тесты для солвера quickest_ultimate_fv_solver
class QUICKEST_ULTIMATE2 : public ::testing::Test {
protected:
    // Профиль переменных
    typedef quickest_ultimate_fv_solver_traits<1>::var_layer_data target_var_t;
    typedef quickest_ultimate_fv_solver_traits<1>::specific_layer specific_data_t;

    // Слой: переменных Vars + сколько угодно служебных Specific
    typedef composite_layer_t<target_var_t, specific_data_t> layer_t;
protected:
    /// @brief Параметры трубы
    pipe_noniso_properties_t pipe;
    oil_parameters_t oil;
    /// @brief Профиль расхода
    std::vector<double> G;
    std::unique_ptr<PipeHeatInflowConstArea> heatModel;
    std::unique_ptr<ring_buffer_t<layer_t>> buffer;
protected:

    /// @brief Подготовка к расчету для семейства тестов
    virtual void SetUp() override {
        oil_parameters_t oil = get_default_oil_heatmodel();
        pipe_noniso_properties_t  pipe = get_default_pipe_heatmodel();
        auto model = pipe.get_heat_eqivalent_model(oil);
        pipe.heat.ambient_heat_transfer = -model.A;    // Использовать пока Кт константу 1-5
        

        std::vector<double> G(pipe.profile.get_point_count(), 300);
        PipeHeatInflowConstArea heatModel(pipe, oil, G);

        buffer = std::make_unique<ring_buffer_t<layer_t>>(2, pipe.profile.get_point_count());

        layer_t& prev = buffer->previous();
        prev.vars.cell_double[0] = std::vector<double>(prev.vars.cell_double[0].size(), 300);
    }
};

/// @brief Базовый пример использования метода характеристик для уравнения адвекции
TEST_F(QUICKEST_ULTIMATE2, Quick_UseCase_Advection_Temperature)
{
    std::vector<double> x = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 };

    double temp_in = 310; // темп нефти, закачиваемой на входе трубы при положительном расходе
    double temp_out = 290; // темп нефти, закачиваемой с выхода трубы при отрицательном расходе

    double t = 0;
    //const auto& x = heatModel->get_grid();
    //double dx = x[1] - x[0];
    //double v = heatModel->getEquationsCoeffs(0, 0);
    //double dt_ideal = abs(dx / v);
    //double Cr = 0.9;
    double dt = 1000;
    oil_parameters_t oil = pde_solvers::get_default_oil_heatmodel();
    pipe_noniso_properties_t  pipe = pde_solvers::get_default_pipe_heatmodel();
    auto model = pipe.get_heat_eqivalent_model(oil);
    //pipe.heat.ambient_heat_transfer = -model.A;
    pipe.heat.ambient_heat_transfer = 3;

    std::vector<double> G(pipe.profile.get_point_count(), 300);
    heatModel = std::make_unique<PipeHeatInflowConstArea>(pipe, oil, G);
    buffer = std::make_unique<ring_buffer_t<layer_t>>(2, pipe.profile.get_point_count());
    layer_t& prev = buffer->previous();
    prev.vars.cell_double[0] = std::vector<double>(prev.vars.cell_double[0].size(), 300);

    for (size_t index = 0; index < 10; ++index) {

        quickest_ultimate_fv_solver solver(*heatModel, *buffer);
        t += dt;
        solver.step(dt, temp_in, temp_out);

        //auto& curr = buffer[0].layers;
        layer_t& next = buffer->current();
        buffer->advance(+1);
    }
}

/// @brief Расчет шага на короткой трубе
/// @return Значение параметра в единственной ячейки на предыдущем слое и на новом
inline std::pair<double, double> short_pipe_step(double vol_flow, 
    double rho_in, double rho_out, double rho_initial)
{
    simple_pipe_properties simple_pipe;
    simple_pipe.length = 100.0;   // длина 100 м
    simple_pipe.dx = 50.0;       // шаг сетки = длина (1 ячейка)
    simple_pipe.diameter = 0.7;
    pipe_properties_t pipe = pipe_properties_t::build_simple_pipe(simple_pipe);

    std::vector<double> Q(pipe.profile.get_point_count(), 1.0);
    PipeQAdvection advection_model(pipe, Q);

    typedef quickest_ultimate_fv_solver_traits<1>::var_layer_data target_var_t;
    typedef quickest_ultimate_fv_solver_traits<1>::specific_layer specific_data_t;
    typedef composite_layer_t<target_var_t, specific_data_t> layer_t;

    ring_buffer_t<layer_t> buffer(2, pipe.profile.get_point_count());
    ring_buffer_t<layer_t> buffer_reverse(2, pipe.profile.get_point_count());

    layer_t& prev = buffer.previous();
    layer_t& next = buffer.current();

    prev.vars.cell_double[0] = std::vector<double>(1, rho_initial);

    double v = advection_model.getEquationsCoeffs(0, Q[0]);
    double Cr = 1;
    double dt = abs((Cr * simple_pipe.dx) / v);       // шаг по времени

    // Запускаем шаг по QUICKEST-ULTIMATE
    quickest_ultimate_fv_solver solver(advection_model, prev, next);
    solver.step(dt, rho_in, rho_out);

    double rho_prev = prev.vars.cell_double[0].front();
    double rho_curr = next.vars.cell_double[0].front();

    return std::make_pair(rho_prev, rho_curr);

}

/// @brief Проверка QUICKEST-ULTIMATE на короткой трубе (1 ячейка)
/// Прямой поток
TEST(QuickestUltimate, HandlesShortPipe_Direct)
{
    double vol_flow = +1.0;
    double rho_in = 860;    // плотность на входе
    double rho_out = 850; // плотность на выходе
    double rho_initial = 840;
    auto [rho_prev, rho_curr] = short_pipe_step(vol_flow, rho_in, rho_out, rho_initial);

    // Проверка
    ASSERT_NEAR(rho_curr, rho_in, 1e-6)
        << "Плотность в текущей ячейке должна стремиться к граничному значению на входе";
    ASSERT_NE(rho_curr, rho_prev)
        << "Плотность не должна оставаться равной начальному значению";
}


/// @brief Проверка QUICKEST-ULTIMATE на короткой трубе (1 ячейка)
/// Обратный поток
TEST(QuickestUltimate, HandlesShortPipe_Reverse)
{
    double vol_flow = -1.0;
    double rho_in = 860;    // плотность на входе
    double rho_out = 850; // плотность на выходе
    double rho_initial = 840;
    auto [rho_prev, rho_curr] = short_pipe_step(vol_flow, rho_in, rho_out, rho_initial);

    // Проверка
    ASSERT_NEAR(rho_curr, rho_in, 1e-6)
        << "Плотность в текущей ячейке должна стремиться к граничному значению на входе";
    ASSERT_NE(rho_curr, rho_prev)
        << "Плотность не должна оставаться равной начальному значению";

}


/// @brief Тест для проверки консервативности численной схемы QUICKEST-ULTIMATE (для одной трубы) 
TEST(QuickestUltimate, SinglePipeConservativityTest)
{
    // === ARRANGE ===
    const double pipe_length = 50e3;/// Длина трубы [м]
    const double dx = 200.0;/// Шаг пространственной сетки [м]
    const double diameter = 0.7;/// Диаметр трубы [м]
    const double total_time = 120 * 3600;///Общее время моделирования [с]


    const double rho_first = 600.0;///Плотность первой партии [кг/м³]
    const double rho_second = 650.0;///Плотность второй партии [кг/м³]

    const double flow_rate = 0.5;///Скорость потока [м³/с]

    const int n_cells_int = static_cast<int>(pipe_length / dx);///Число ячеек сетки

    const double switch_time = 6 * 3600;///Время через которое приходит вторая партия
    double dt = 60;///Шаг по времени [с]
    int n_steps = static_cast<int>(total_time / dt);///Число шагов по времени

    const double area = 3.141592653589793 * diameter * diameter / 4.0;///Площадь сечения [м²]
    const double cell_volume = area * dx; ///Объем одной ячейки[м³]

    double initial_mass_in_pipe = 0.0;

    simple_pipe_properties simple_pipe;
    simple_pipe.length = pipe_length;
    simple_pipe.dx = dx;
    simple_pipe.diameter = diameter;
    pipe_properties_t pipe = pipe_properties_t::build_simple_pipe(simple_pipe);
    
    std::vector<double> Q(pipe.profile.get_point_count(), flow_rate);///Вектор расходов
    PipeQAdvection advection_model(pipe, Q);///Модель адвекции

    typedef quickest_ultimate_fv_solver_traits<1>::var_layer_data target_var_t;
    typedef quickest_ultimate_fv_solver_traits<1>::specific_layer specific_data_t;
    typedef composite_layer_t<target_var_t, specific_data_t> layer_t;

    ring_buffer_t<layer_t> buffer(2, pipe.profile.get_point_count());
    layer_t& prev = buffer.previous();
    prev.vars.cell_double[0] = std::vector<double>(prev.vars.cell_double[0].size(), rho_first);

    double v = 0.0;
    if (!Q.empty()) {
        v = advection_model.getEquationsCoeffs(0, Q[0]);
    }    
    const auto& initial_densities = prev.vars.cell_double[0];
    for (int i = 0; i < n_cells_int; ++i) {
        if (static_cast<size_t>(i) < initial_densities.size()) {
            initial_mass_in_pipe += initial_densities[i] * cell_volume;
        }
    }

    double cumulative_mass_in = 0.0;
    double cumulative_mass_out = 0.0;

    std::vector<double> mass_in_pipe_history, mass_balance_history, time_history, balance_error_history, rho_outlet_history;

    mass_in_pipe_history.push_back(initial_mass_in_pipe);
    mass_balance_history.push_back(initial_mass_in_pipe);
    time_history.push_back(0.0);
    balance_error_history.push_back(0.0);
    rho_outlet_history.push_back(rho_first);

    double current_time = 0.0;
    int switch_step = static_cast<int>(switch_time / dt);
    double current_mass_in_pipe = initial_mass_in_pipe;

    // === ACT ===
    /// Основной цикл расчетов
    for (int step = 1; step <= n_steps; ++step) {
        current_time += dt;

        double rho_in = (step <= switch_step) ? rho_first : rho_second;

        double rho_outlet_current = 0.0;
        if (!prev.vars.cell_double[0].empty()) {
            rho_outlet_current = prev.vars.cell_double[0].back();
        }
        else {
            rho_outlet_current = rho_first;
        }

        double mass_flow_in = rho_in * flow_rate * dt;
        double mass_flow_out = rho_outlet_current * flow_rate * dt;

        cumulative_mass_in += mass_flow_in;
        cumulative_mass_out += mass_flow_out;

        double mass_by_balance = initial_mass_in_pipe + cumulative_mass_in - cumulative_mass_out;

        try {
            layer_t& next = buffer.current();
            quickest_ultimate_fv_solver solver(advection_model, prev, next);

            solver.step(dt, rho_in, rho_outlet_current);

            current_mass_in_pipe = 0.0;
            const auto& densities = next.vars.cell_double[0];
            for (int i = 0; i < n_cells_int; ++i) {
                if (static_cast<size_t>(i) < densities.size()) {
                    current_mass_in_pipe += densities[i] * cell_volume;
                }
            }

            double balance_error = current_mass_in_pipe - mass_by_balance;

            mass_in_pipe_history.push_back(current_mass_in_pipe);
            mass_balance_history.push_back(mass_by_balance);
            time_history.push_back(current_time);
            balance_error_history.push_back(balance_error);
            rho_outlet_history.push_back(rho_outlet_current);

            buffer.advance(+1);
            prev = buffer.previous();
        }
        catch (const std::exception& e) {
            FAIL() << "Ошибка на шаге " << step << ": " << e.what();
        }
    }

    // === ASSERT ===
    double final_mass_by_scheme = current_mass_in_pipe;
    double final_mass_by_balance = initial_mass_in_pipe + cumulative_mass_in - cumulative_mass_out;
    double absolute_error = final_mass_by_scheme - final_mass_by_balance;

    double ref1 = std::abs(initial_mass_in_pipe);
    double ref2 = std::abs(final_mass_by_scheme);
    double ref3 = std::abs(cumulative_mass_in);
    double ref4 = std::abs(cumulative_mass_out);

    double max_mass = ref1;
    if (ref2 > max_mass) max_mass = ref2;
    if (ref3 > max_mass) max_mass = ref3;
    if (ref4 > max_mass) max_mass = ref4;

    if (max_mass < 1.0) max_mass = 1.0;
    double relative_error = std::abs(absolute_error) / max_mass;

    const double absolute_tolerance = 1e-6;
    const double relative_tolerance = 1e-10;

    double tolerance1 = absolute_tolerance;
    double tolerance2 = relative_tolerance * max_mass;
    double tolerance = (tolerance1 > tolerance2) ? tolerance1 : tolerance2;

    EXPECT_NEAR(final_mass_by_scheme, final_mass_by_balance, tolerance)
        << "Схема QUICKEST-ULTIMATE не консервативна. "
        << "M_расчетная = " << final_mass_by_scheme << " кг, "
        << "M_балансная = " << final_mass_by_balance << " кг, "
        << "разница = " << absolute_error << " кг";

    // Сохранение в CSV файл
    if (auto csv = std::ofstream("conservativity_results.csv")) {
        csv << "step,time_h,mass_scheme_kg,mass_balance_kg,error_kg,rho_out_kg_m3\n";
        for (size_t i = 0; i < time_history.size(); ++i)
            csv << i << ';' << time_history[i] / 3600 << ';' << mass_in_pipe_history[i] << ';'
            << mass_balance_history[i] << ';' << balance_error_history[i] << ';'
            << rho_outlet_history[i] << '\n';
    }
}


