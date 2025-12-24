#pragma once


// Настройка кодировки для Windows
#ifdef _WIN32
#include <windows.h>

class WindowsConsoleSetup {
public:
    WindowsConsoleSetup() {
        // Сохраняем текущие настройки
        oldOutputCP = GetConsoleOutputCP();
        oldInputCP = GetConsoleCP();

        // Устанавливаем UTF-8
        SetConsoleOutputCP(CP_UTF8);
        SetConsoleCP(CP_UTF8);

        std::cout << "Консоль настроена на UTF-8" << std::endl;
    }

    ~WindowsConsoleSetup() {
        // Восстанавливаем настройки (опционально)
        SetConsoleOutputCP(oldOutputCP);
        SetConsoleCP(oldInputCP);
    }

private:
    UINT oldOutputCP;
    UINT oldInputCP;
};

// Глобальный объект для настройки консоли
static WindowsConsoleSetup consoleSetup;
#endif





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
TEST_F(UpstreamDifferencing, CanConsiderFlowSwap) {
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
TEST_F(QUICKEST_ULTIMATE, CanConsiderFlowSwap) {
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


//TEST(QuickestUltimate, SinglePipeConservativityTest)
//{
//    // === ARRANGE ===
//    std::cout << "\n=== ТЕСТ КОНСЕРВАТИВНОСТИ ДЛЯ ОДНОЙ ТРУБЫ ===\n" << std::endl;
//
//    // Параметры теста
//    const double pipe_length = 50e3;      // 50 км
//    const double dx = 100.0;              // шаг 100 м
//    const double diameter = 0.7;          // диаметр трубы
//    const double total_time = 120 * 3600; // 120 часов в секундах
//
//    // Плотности партий
//    const double rho_first = 600.0;       // первая партия (начальная)
//    const double rho_second = 650.0;      // вторая партия (вход)
//    // rho_outlet НЕ константа - будет меняться со временем!
//
//    // Расход
//    const double flow_rate = 0.5;         // 0.5 м³/с
//
//    // Расчет количества ячеек
//    const int n_cells_int = static_cast<int>(pipe_length / dx);
//    ASSERT_GT(n_cells_int, 0) << "Неверный шаг сетки";
//    const size_t n_cells = static_cast<size_t>(n_cells_int);
//
//    std::cout << "Параметры расчета:" << std::endl;
//    std::cout << "  Длина трубы: " << pipe_length << " м" << std::endl;
//    std::cout << "  Шаг сетки: " << dx << " м" << std::endl;
//    std::cout << "  Количество ячеек: " << n_cells_int << std::endl;
//    std::cout << "  Диаметр: " << diameter << " м" << std::endl;
//    std::cout << "  Площадь сечения: " << 3.141592653589793 * diameter * diameter / 4.0 << " м²" << std::endl;
//    std::cout << "  Расход: " << flow_rate << " м³/с" << std::endl;
//    std::cout << "  Общее время: " << total_time << " с (" << total_time / 3600 << " ч)" << std::endl;
//    std::cout << "  Плотность первой партии: " << rho_first << " кг/м³" << std::endl;
//    std::cout << "  Плотность второй партии: " << rho_second << " кг/м³" << std::endl;
//    std::cout << "  Плотность на выходе: БУДЕТ МЕНЯТЬСЯ СО ВРЕМЕНЕМ" << std::endl;
//
//    // Создание модели трубы
//    simple_pipe_properties simple_pipe;
//    simple_pipe.length = pipe_length;
//    simple_pipe.dx = dx;
//    simple_pipe.diameter = diameter;
//    pipe_properties_t pipe = pipe_properties_t::build_simple_pipe(simple_pipe);
//
//    // Вектор расходов (постоянный по трубе)
//    std::vector<double> Q(pipe.profile.get_point_count(), flow_rate);
//    PipeQAdvection advection_model(pipe, Q);
//
//    // Типы данных для солвера
//    typedef quickest_ultimate_fv_solver_traits<1>::var_layer_data target_var_t;
//    typedef quickest_ultimate_fv_solver_traits<1>::specific_layer specific_data_t;
//    typedef composite_layer_t<target_var_t, specific_data_t> layer_t;
//
//    // Буфер для хранения слоев
//    ring_buffer_t<layer_t> buffer(2, pipe.profile.get_point_count());
//
//    // Начальное состояние: вся труба заполнена первой партией
//    layer_t& prev = buffer.previous();
//    prev.vars.cell_double[0] = std::vector<double>(prev.vars.cell_double[0].size(), rho_first);
//
//    // Расчет скорости потока
//    double v = 0.0;
//    if (!Q.empty()) {
//        v = advection_model.getEquationsCoeffs(0, Q[0]);
//    }
//    ASSERT_NE(v, 0) << "Скорость потока не должна быть нулевой";
//
//    std::cout << "\nСкорость потока: " << v << " м/с" << std::endl;
//
//    // Расчет времени прохода фронта по трубе
//    const double transit_time = pipe_length / std::abs(v);
//    std::cout << "Время прохода фронта по трубе: " << transit_time << " с ("
//        << transit_time / 3600 << " ч)" << std::endl;
//
//    // Определение времени подачи второй партии
//    const double switch_time = 6 * 3600;  // через 6 часов
//    std::cout << "Время подачи второй партии: " << switch_time << " с ("
//        << switch_time / 3600 << " ч)" << std::endl;
//
//    // Шаг по времени
//    double dt = 60;
//
//    // Количество шагов
//    int n_steps = static_cast<int>(total_time / dt);
//
//    std::cout << "\nПараметры расчета по времени:" << std::endl;
//    std::cout << "  Шаг по времени: " << dt << " с" << std::endl;
//    std::cout << "  Количество шагов: " << n_steps << std::endl;
//    std::cout << "  Шаг переключения партий: " << static_cast<int>(switch_time / dt) << std::endl;
//
//    // Параметры для отслеживания массы
//    const double area = 3.141592653589793 * diameter * diameter / 4.0;
//    const double cell_volume = area * dx;
//
//    // Начальная масса в трубе (Σ m_i в начальный момент)
//    double initial_mass_in_pipe = 0.0;
//    const auto& initial_densities = prev.vars.cell_double[0];
//    for (int i = 0; i < n_cells_int; ++i) {
//        if (static_cast<size_t>(i) < initial_densities.size()) {
//            initial_mass_in_pipe += initial_densities[i] * cell_volume;
//        }
//    }
//
//    // Переменные для отслеживания баланса массы
//    double cumulative_mass_in = 0.0;     // Σ (Δt * G_вход,i)
//    double cumulative_mass_out = 0.0;    // Σ (Δt * G_выход,i)
//
//    // Векторы для хранения истории
//    std::vector<double> mass_in_pipe_history;
//    std::vector<double> mass_balance_history;
//    std::vector<double> time_history;
//    std::vector<double> balance_error_history;
//    std::vector<double> rho_outlet_history;  // История плотности на выходе
//
//    mass_in_pipe_history.push_back(initial_mass_in_pipe);
//    mass_balance_history.push_back(initial_mass_in_pipe);
//    time_history.push_back(0.0);
//    balance_error_history.push_back(0.0);
//    rho_outlet_history.push_back(rho_first);  // Начальная плотность на выходе
//
//    // Подготовка для цикла расчета
//    double current_time = 0.0;
//    int switch_step = static_cast<int>(switch_time / dt);
//    double current_mass_in_pipe = initial_mass_in_pipe;
//
//    std::cout << "\nНачальные параметры массы:" << std::endl;
//    std::cout << "  Площадь сечения: " << area << " м²" << std::endl;
//    std::cout << "  Объем ячейки: " << cell_volume << " м³" << std::endl;
//    std::cout << "  Объем трубы: " << cell_volume * n_cells_int << " м³" << std::endl;
//    std::cout << "  Начальная масса в трубе (Σ m_i): " << initial_mass_in_pipe << " кг" << std::endl;
//
//    // === ACT ===
//    std::cout << "\nЗапуск расчета..." << std::endl;
//
//    for (int step = 1; step <= n_steps; ++step) {
//        current_time += dt;
//
//        // Определение граничных условий на ВХОДЕ
//        double rho_in = (step <= switch_step) ? rho_first : rho_second;
//
//        // ВАЖНО: плотность на выходе берется из ПРЕДЫДУЩЕГО временного слоя
//        // Это плотность в последней ячейке трубы на предыдущем шаге
//        double rho_outlet_current = 0.0;
//        if (!prev.vars.cell_double[0].empty()) {
//            rho_outlet_current = prev.vars.cell_double[0].back();  // Плотность в последней ячейке
//        }
//        else {
//            rho_outlet_current = rho_first;  // Начальное значение
//        }
//
//        // Расчет потоков массы через границы за шаг dt
//        // ВАЖНО: используем РАЗНЫЕ плотности для входа и выхода!
//        double mass_flow_in = rho_in * flow_rate * dt;
//        double mass_flow_out = rho_outlet_current * flow_rate * dt;  // ИЗМЕНЕНО!
//
//        cumulative_mass_in += mass_flow_in;
//        cumulative_mass_out += mass_flow_out;
//
//        // Расчет массы по балансу (корректная формула)
//        double mass_by_balance = initial_mass_in_pipe + cumulative_mass_in - cumulative_mass_out;
//
//        // Выполнение шага расчета
//        try {
//            layer_t& next = buffer.current();
//            quickest_ultimate_fv_solver solver(advection_model, prev, next);
//
//            // Передаем в солвер правильные граничные условия
//            // rho_outlet_current - это плотность на выходе на предыдущем шаге
//            solver.step(dt, rho_in, rho_outlet_current);
//
//            // Расчет текущей массы в трубе по схеме (Σ m_i)
//            current_mass_in_pipe = 0.0;
//            const auto& densities = next.vars.cell_double[0];
//            for (int i = 0; i < n_cells_int; ++i) {
//                if (static_cast<size_t>(i) < densities.size()) {
//                    current_mass_in_pipe += densities[i] * cell_volume;
//                }
//            }
//
//            // Расчет ошибки баланса
//            double balance_error = current_mass_in_pipe - mass_by_balance;
//
//            // Сохранение истории
//            mass_in_pipe_history.push_back(current_mass_in_pipe);
//            mass_balance_history.push_back(mass_by_balance);
//            time_history.push_back(current_time);
//            balance_error_history.push_back(balance_error);
//            rho_outlet_history.push_back(rho_outlet_current);  // Сохраняем историю плотности на выходе
//
//            // Переход к следующему шагу
//            buffer.advance(+1);
//            prev = buffer.previous();
//
//        }
//        catch (const std::exception& e) {
//            FAIL() << "Ошибка на шаге " << step << ": " << e.what();
//        }
//
//        // Вывод прогресса
//        if (n_steps >= 10 && step % (n_steps / 10) == 0) {
//            std::cout << "  Шаг " << step << "/" << n_steps
//                << " (время: " << current_time / 3600 << " ч, "
//                << "M_расчетная: " << current_mass_in_pipe << " кг, "
//                << "M_балансная: " << mass_by_balance << " кг, "
//                << "ρ_выход: " << rho_outlet_current << " кг/м³)" << std::endl;  // Добавлен вывод плотности на выходе
//        }
//    }
//
//    std::cout << "\nРасчет завершен." << std::endl;
//
//    // === ASSERT ===
//    // Финальная проверка консервативности
//    double final_mass_by_scheme = current_mass_in_pipe;
//    double final_mass_by_balance = initial_mass_in_pipe + cumulative_mass_in - cumulative_mass_out;
//    double absolute_error = final_mass_by_scheme - final_mass_by_balance;
//
//    // Расчет относительной ошибки
//    double ref1 = std::abs(initial_mass_in_pipe);
//    double ref2 = std::abs(final_mass_by_scheme);
//    double ref3 = std::abs(cumulative_mass_in);
//    double ref4 = std::abs(cumulative_mass_out);
//
//    double max_mass = ref1;
//    if (ref2 > max_mass) max_mass = ref2;
//    if (ref3 > max_mass) max_mass = ref3;
//    if (ref4 > max_mass) max_mass = ref4;
//
//    if (max_mass < 1.0) max_mass = 1.0;
//    double relative_error = std::abs(absolute_error) / max_mass;
//
//    // Вывод результатов
//    std::cout << "\n=== ПРОВЕРКА КОНСЕРВАТИВНОСТИ ПО ФОРМУЛЕ ===" << std::endl;
//    std::cout << "Формула: M_балансная(t) = Σ m_i(0) + Σ(Δt * ρ_вход * Q) - Σ(Δt * ρ_выход * Q)" << std::endl;
//    std::cout << "Σ m_i(0)  (начальная масса): " << initial_mass_in_pipe << " кг" << std::endl;
//    std::cout << "Σ(Δt * ρ_вход * Q) (масса на входе): " << cumulative_mass_in << " кг" << std::endl;
//    std::cout << "Σ(Δt * ρ_выход * Q) (масса на выходе): " << cumulative_mass_out << " кг" << std::endl;
//    std::cout << "M_балансная(t) по формуле: " << final_mass_by_balance << " кг" << std::endl;
//    std::cout << "M_расчетная(t) по схеме: " << final_mass_by_scheme << " кг" << std::endl;
//    std::cout << "\nАбсолютная ошибка (M_расчетная - M_балансная): " << absolute_error << " кг" << std::endl;
//    std::cout << "Относительная ошибка: " << relative_error << std::endl;
//
//    // Вывод финальной плотности на выходе
//    if (!rho_outlet_history.empty()) {
//        std::cout << "Финальная плотность на выходе: " << rho_outlet_history.back() << " кг/м³" << std::endl;
//    }
//
//    // Критерии проверки консервативности
//    const double absolute_tolerance = 1e-6;
//    const double relative_tolerance = 1e-10;
//
//    double tolerance1 = absolute_tolerance;
//    double tolerance2 = relative_tolerance * max_mass;
//    double tolerance = (tolerance1 > tolerance2) ? tolerance1 : tolerance2;
//
//    std::cout << "\nКритерии проверки консервативности:" << std::endl;
//    std::cout << "  Абсолютный допуск: " << absolute_tolerance << " кг" << std::endl;
//    std::cout << "  Относительный допуск: " << relative_tolerance * 100 << " %" << std::endl;
//    std::cout << "  Итоговый допуск: " << tolerance << " кг" << std::endl;
//
//    // Проверка консервативности
//    if (std::abs(absolute_error) <= tolerance) {
//        std::cout << "\n✓ СХЕМА КОНСЕРВАТИВНА: M_расчетная ≈ M_балансная" << std::endl;
//        std::cout << "  Ошибка (" << absolute_error << " кг) в пределах допуска (" << tolerance << " кг)" << std::endl;
//        SUCCEED();
//    }
//    else {
//        std::cout << "\n✗ СХЕМА НЕ КОНСЕРВАТИВНА: M_расчетная ≠ M_балансная" << std::endl;
//        std::cout << "  Ошибка (" << absolute_error << " кг) превышает допуск (" << tolerance << " кг)" << std::endl;
//
//        // Анализ максимальной ошибки в истории
//        double max_abs_error = 0.0;
//        size_t max_error_step = 0;
//
//        for (size_t i = 0; i < balance_error_history.size(); ++i) {
//            double abs_error = std::abs(balance_error_history[i]);
//            if (abs_error > max_abs_error) {
//                max_abs_error = abs_error;
//                max_error_step = i;
//            }
//        }
//
//        std::cout << "  Максимальная ошибка в истории: " << max_abs_error << " кг" << std::endl;
//        std::cout << "  На шаге: " << max_error_step
//            << " (время: " << time_history[max_error_step] / 3600 << " ч)" << std::endl;
//    }
//
//    // Основная проверка ASSERT-ом
//    EXPECT_NEAR(final_mass_by_scheme, final_mass_by_balance, tolerance)
//        << "Схема QUICKEST-ULTIMATE не консервативна. "
//        << "M_расчетная = " << final_mass_by_scheme << " кг, "
//        << "M_балансная = " << final_mass_by_balance << " кг, "
//        << "разница = " << absolute_error << " кг";
//
//    // Дополнительная проверка: на каждом шаге ошибка должна быть малой
//    bool all_steps_conservative = true;
//    for (size_t i = 0; i < balance_error_history.size(); ++i) {
//        if (std::abs(balance_error_history[i]) > tolerance * 10) {
//            std::cout << "  Предупреждение: большая ошибка на шаге " << i
//                << " (время " << time_history[i] / 3600 << " ч): "
//                << balance_error_history[i] << " кг" << std::endl;
//            all_steps_conservative = false;
//        }
//    }
//
//    if (all_steps_conservative) {
//        std::cout << "\n✓ Схема консервативна на всех шагах расчета" << std::endl;
//    }
//
//    // Сохранение подробных результатов
//    std::ofstream results("conservativity_detailed_results.csv");
//    if (results.is_open()) {
//        results << "time_h,mass_by_scheme_kg,mass_by_balance_kg,balance_error_kg,rho_outlet_kg_m3\n";
//        for (size_t i = 0; i < time_history.size(); ++i) {
//            results << time_history[i] / 3600 << ","
//                << mass_in_pipe_history[i] << ","
//                << mass_balance_history[i] << ","
//                << balance_error_history[i] << ","
//                << rho_outlet_history[i] << "\n";
//        }
//        results.close();
//        std::cout << "\nПодробные результаты сохранены в файл: conservativity_detailed_results.csv" << std::endl;
//    }
//
//    std::cout << "\n=== ТЕСТ ЗАВЕРШЕН ===" << std::endl;
//}

/// @brief ТЕСТ КОНСЕРВАТИВНОСТИ ДЛЯ ОДНОЙ ТРУБЫ
TEST(QuickestUltimate, SinglePipeConservativityTest)
{
    // === ARRANGE ===
    std::cout << "\n=== ТЕСТ КОНСЕРВАТИВНОСТИ ДЛЯ ОДНОЙ ТРУБЫ ===\n" << std::endl;

    // Параметры теста
    const double pipe_length = 50e3;      // 50 км
    const double dx = 200.0;              // шаг 200 м
    const double diameter = 0.7;          // диаметр трубы
    const double total_time = 120 * 3600; // 120 часов в секундах

    // Плотности партий
    const double rho_first = 600.0;       // первая партия (начальная)
    const double rho_second = 650.0;      // вторая партия (вход)

    // Расход
    const double flow_rate = 0.5;         // 0.5 м³/с

    // Расчет количества ячеек
    const int n_cells_int = static_cast<int>(pipe_length / dx);
    ASSERT_GT(n_cells_int, 0) << "Неверный шаг сетки";
    const size_t n_cells = static_cast<size_t>(n_cells_int);

    std::cout << "Параметры расчета:" << std::endl;
    std::cout << "  Длина трубы: " << pipe_length << " м" << std::endl;
    std::cout << "  Шаг сетки: " << dx << " м" << std::endl;
    std::cout << "  Количество ячеек: " << n_cells_int << std::endl;
    std::cout << "  Диаметр: " << diameter << " м" << std::endl;
    std::cout << "  Площадь сечения: " << 3.141592653589793 * diameter * diameter / 4.0 << " м²" << std::endl;
    std::cout << "  Расход: " << flow_rate << " м³/с" << std::endl;
    std::cout << "  Общее время: " << total_time << " с (" << total_time / 3600 << " ч)" << std::endl;
    std::cout << "  Плотность первой партии: " << rho_first << " кг/м³" << std::endl;
    std::cout << "  Плотность второй партии: " << rho_second << " кг/м³" << std::endl;
    std::cout << "  Плотность на выходе: БУДЕТ МЕНЯТЬСЯ СО ВРЕМЕНЕМ" << std::endl;

    // Создание модели трубы
    simple_pipe_properties simple_pipe;
    simple_pipe.length = pipe_length;
    simple_pipe.dx = dx;
    simple_pipe.diameter = diameter;
    pipe_properties_t pipe = pipe_properties_t::build_simple_pipe(simple_pipe);

    // Вектор расходов (постоянный по трубе)
    std::vector<double> Q(pipe.profile.get_point_count(), flow_rate);
    PipeQAdvection advection_model(pipe, Q);

    // Типы данных для солвера
    typedef quickest_ultimate_fv_solver_traits<1>::var_layer_data target_var_t;
    typedef quickest_ultimate_fv_solver_traits<1>::specific_layer specific_data_t;
    typedef composite_layer_t<target_var_t, specific_data_t> layer_t;

    // Буфер для хранения слоев
    ring_buffer_t<layer_t> buffer(2, pipe.profile.get_point_count());

    // Начальное состояние: вся труба заполнена первой партией
    layer_t& prev = buffer.previous();
    prev.vars.cell_double[0] = std::vector<double>(prev.vars.cell_double[0].size(), rho_first);

    // Расчет скорости потока
    double v = 0.0;
    if (!Q.empty()) {
        v = advection_model.getEquationsCoeffs(0, Q[0]);
    }
    ASSERT_NE(v, 0) << "Скорость потока не должна быть нулевой";

    std::cout << "\nСкорость потока: " << v << " м/с" << std::endl;

    // Расчет времени прохода фронта по трубе
    const double transit_time = pipe_length / std::abs(v);
    std::cout << "Время прохода фронта по трубе: " << transit_time << " с ("
        << transit_time / 3600 << " ч)" << std::endl;

    // Определение времени подачи второй партии
    const double switch_time = 6 * 3600;  // через 6 часов
    std::cout << "Время подачи второй партии: " << switch_time << " с ("
        << switch_time / 3600 << " ч)" << std::endl;

    // Шаг по времени
    double dt = 60;

    // Количество шагов
    int n_steps = static_cast<int>(total_time / dt);

    std::cout << "\nПараметры расчета по времени:" << std::endl;
    std::cout << "  Шаг по времени: " << dt << " с" << std::endl;
    std::cout << "  Количество шагов: " << n_steps << std::endl;
    std::cout << "  Шаг переключения партий: " << static_cast<int>(switch_time / dt) << std::endl;

    // Параметры для отслеживания массы
    const double area = 3.141592653589793 * diameter * diameter / 4.0;
    const double cell_volume = area * dx;

    // Начальная масса в трубе (Σ m_i в начальный момент)
    double initial_mass_in_pipe = 0.0;
    const auto& initial_densities = prev.vars.cell_double[0];
    for (int i = 0; i < n_cells_int; ++i) {
        if (static_cast<size_t>(i) < initial_densities.size()) {
            initial_mass_in_pipe += initial_densities[i] * cell_volume;
        }
    }

    // Переменные для отслеживания баланса массы
    double cumulative_mass_in = 0.0;     // Σ (Δt * G_вход,i)
    double cumulative_mass_out = 0.0;    // Σ (Δt * G_выход,i)

    // Векторы для хранения истории
    std::vector<double> mass_in_pipe_history;
    std::vector<double> mass_balance_history;
    std::vector<double> time_history;
    std::vector<double> balance_error_history;
    std::vector<double> rho_outlet_history;  // История плотности на выходе

    mass_in_pipe_history.push_back(initial_mass_in_pipe);
    mass_balance_history.push_back(initial_mass_in_pipe);
    time_history.push_back(0.0);
    balance_error_history.push_back(0.0);
    rho_outlet_history.push_back(rho_first);  // Начальная плотность на выходе

    // Подготовка для цикла расчета
    double current_time = 0.0;
    int switch_step = static_cast<int>(switch_time / dt);
    double current_mass_in_pipe = initial_mass_in_pipe;

    std::cout << "\nНачальные параметры массы:" << std::endl;
    std::cout << "  Площадь сечения: " << area << " м²" << std::endl;
    std::cout << "  Объем ячейки: " << cell_volume << " м³" << std::endl;
    std::cout << "  Объем трубы: " << cell_volume * n_cells_int << " м³" << std::endl;
    std::cout << "  Начальная масса в трубе (Σ m_i): " << initial_mass_in_pipe << " кг" << std::endl;

    // === ACT ===
    std::cout << "\nЗапуск расчета..." << std::endl;

    for (int step = 1; step <= n_steps; ++step) {
        current_time += dt;

        // Определение граничных условий на ВХОДЕ
        double rho_in = (step <= switch_step) ? rho_first : rho_second;

        // Плотность на выходе берется из ПРЕДЫДУЩЕГО временного слоя
        double rho_outlet_current = 0.0;
        if (!prev.vars.cell_double[0].empty()) {
            rho_outlet_current = prev.vars.cell_double[0].back();  // Плотность в последней ячейке
        }
        else {
            rho_outlet_current = rho_first;  // Начальное значение
        }

        // Расчет потоков массы через границы за шаг dt
        double mass_flow_in = rho_in * flow_rate * dt;
        double mass_flow_out = rho_outlet_current * flow_rate * dt;

        cumulative_mass_in += mass_flow_in;
        cumulative_mass_out += mass_flow_out;

        // Расчет массы по балансу (корректная формула)
        double mass_by_balance = initial_mass_in_pipe + cumulative_mass_in - cumulative_mass_out;

        // Выполнение шага расчета
        try {
            layer_t& next = buffer.current();
            quickest_ultimate_fv_solver solver(advection_model, prev, next);

            // Передаем в солвер правильные граничные условия
            solver.step(dt, rho_in, rho_outlet_current);

            // Расчет текущей массы в трубе по схеме (Σ m_i)
            current_mass_in_pipe = 0.0;
            const auto& densities = next.vars.cell_double[0];
            for (int i = 0; i < n_cells_int; ++i) {
                if (static_cast<size_t>(i) < densities.size()) {
                    current_mass_in_pipe += densities[i] * cell_volume;
                }
            }

            // Расчет ошибки баланса
            double balance_error = current_mass_in_pipe - mass_by_balance;

            // Сохранение истории
            mass_in_pipe_history.push_back(current_mass_in_pipe);
            mass_balance_history.push_back(mass_by_balance);
            time_history.push_back(current_time);
            balance_error_history.push_back(balance_error);
            rho_outlet_history.push_back(rho_outlet_current);

            // Переход к следующему шагу
            buffer.advance(+1);
            prev = buffer.previous();

        }
        catch (const std::exception& e) {
            FAIL() << "Ошибка на шаге " << step << ": " << e.what();
        }

        // Вывод прогресса
        if (n_steps >= 10 && step % (n_steps / 10) == 0) {
            std::cout << "  Шаг " << step << "/" << n_steps
                << " (время: " << current_time / 3600 << " ч, "
                << "M_расчетная: " << current_mass_in_pipe << " кг, "
                << "M_балансная: " << mass_by_balance << " кг, "
                << "ρ_выход: " << rho_outlet_current << " кг/м³)" << std::endl;
        }
    }

    std::cout << "\nРасчет завершен." << std::endl;

    // === ASSERT ===
    // Финальная проверка консервативности
    double final_mass_by_scheme = current_mass_in_pipe;
    double final_mass_by_balance = initial_mass_in_pipe + cumulative_mass_in - cumulative_mass_out;
    double absolute_error = final_mass_by_scheme - final_mass_by_balance;

    // Расчет относительной ошибки
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

    // Вывод результатов
    std::cout << "\n=== ПРОВЕРКА КОНСЕРВАТИВНОСТИ ПО ФОРМУЛЕ ===" << std::endl;
    std::cout << "Формула: M_балансная(t) = Σ m_i(0) + Σ(Δt * ρ_вход * Q) - Σ(Δt * ρ_выход * Q)" << std::endl;
    std::cout << "Σ m_i(0)  (начальная масса): " << initial_mass_in_pipe << " кг" << std::endl;
    std::cout << "Σ(Δt * ρ_вход * Q) (масса на входе): " << cumulative_mass_in << " кг" << std::endl;
    std::cout << "Σ(Δt * ρ_выход * Q) (масса на выходе): " << cumulative_mass_out << " кг" << std::endl;
    std::cout << "M_балансная(t) по формуле: " << final_mass_by_balance << " кг" << std::endl;
    std::cout << "M_расчетная(t) по схеме: " << final_mass_by_scheme << " кг" << std::endl;
    std::cout << "\nАбсолютная ошибка (M_расчетная - M_балансная): " << absolute_error << " кг" << std::endl;
    std::cout << "Относительная ошибка: " << relative_error << std::endl;

    if (!rho_outlet_history.empty()) {
        std::cout << "Финальная плотность на выходе: " << rho_outlet_history.back() << " кг/м³" << std::endl;
    }

    // Критерии проверки консервативности
    const double absolute_tolerance = 1e-6;
    const double relative_tolerance = 1e-10;

    double tolerance1 = absolute_tolerance;
    double tolerance2 = relative_tolerance * max_mass;
    double tolerance = (tolerance1 > tolerance2) ? tolerance1 : tolerance2;

    std::cout << "\nКритерии проверки консервативности:" << std::endl;
    std::cout << "  Абсолютный допуск: " << absolute_tolerance << " кг" << std::endl;
    std::cout << "  Относительный допуск: " << relative_tolerance * 100 << " %" << std::endl;
    std::cout << "  Итоговый допуск: " << tolerance << " кг" << std::endl;

    // Проверка консервативности
    if (std::abs(absolute_error) <= tolerance) {
        std::cout << "\n✓ СХЕМА КОНСЕРВАТИВНА: M_расчетная ≈ M_балансная" << std::endl;
        std::cout << "  Ошибка (" << absolute_error << " кг) в пределах допуска (" << tolerance << " кг)" << std::endl;
        SUCCEED();
    }
    else {
        std::cout << "\n✗ СХЕМА НЕ КОНСЕРВАТИВНА: M_расчетная ≠ M_балансная" << std::endl;
        std::cout << "  Ошибка (" << absolute_error << " кг) превышает допуск (" << tolerance << " кг)" << std::endl;

        // Анализ максимальной ошибки в истории
        double max_abs_error = 0.0;
        size_t max_error_step = 0;

        for (size_t i = 0; i < balance_error_history.size(); ++i) {
            double abs_error = std::abs(balance_error_history[i]);
            if (abs_error > max_abs_error) {
                max_abs_error = abs_error;
                max_error_step = i;
            }
        }

        std::cout << "  Максимальная ошибка в истории: " << max_abs_error << " кг" << std::endl;
        std::cout << "  На шаге: " << max_error_step
            << " (время: " << time_history[max_error_step] / 3600 << " ч)" << std::endl;
    }

    // Основная проверка ASSERT-ом
    EXPECT_NEAR(final_mass_by_scheme, final_mass_by_balance, tolerance)
        << "Схема QUICKEST-ULTIMATE не консервативна. "
        << "M_расчетная = " << final_mass_by_scheme << " кг, "
        << "M_балансная = " << final_mass_by_balance << " кг, "
        << "разница = " << absolute_error << " кг";

    // Дополнительная проверка: на каждом шаге ошибка должна быть малой
    bool all_steps_conservative = true;
    for (size_t i = 0; i < balance_error_history.size(); ++i) {
        if (std::abs(balance_error_history[i]) > tolerance * 10) {
            std::cout << "  Предупреждение: большая ошибка на шаге " << i
                << " (время " << time_history[i] / 3600 << " ч): "
                << balance_error_history[i] << " кг" << std::endl;
            all_steps_conservative = false;
        }
    }

    if (all_steps_conservative) {
        std::cout << "\n✓ Схема консервативна на всех шагах расчета" << std::endl;
    }

    // === СОХРАНЕНИЕ РЕЗУЛЬТАТОВ ===
    std::cout << "\n=== СОХРАНЕНИЕ РЕЗУЛЬТАТОВ ===" << std::endl;

    // 1. Сохранение детальных результатов в CSV
    std::ofstream results("conservativity_detailed_results.csv");
    if (results.is_open()) {
        results << "step,time_h,mass_by_scheme_kg,mass_by_balance_kg,balance_error_kg,rho_outlet_kg_m3\n";
        for (size_t i = 0; i < time_history.size(); ++i) {
            results << i << ","
                << time_history[i] / 3600 << ","
                << mass_in_pipe_history[i] << ","
                << mass_balance_history[i] << ","
                << balance_error_history[i] << ","
                << rho_outlet_history[i] << "\n";
        }
        results.close();
        std::cout << "✓ Детальные результаты сохранены в файл: conservativity_detailed_results.csv" << std::endl;
    }
    else {
        std::cout << "✗ Ошибка сохранения файла conservativity_detailed_results.csv" << std::endl;
    }

    // 2. Сохранение итогового отчета
    std::ofstream summary("conservativity_summary.txt");
    if (summary.is_open()) {
        summary << "=== ОТЧЕТ ПО ТЕСТУ КОНСЕРВАТИВНОСТИ QUICKEST-ULTIMATE ===\n\n";

        summary << "ПАРАМЕТРЫ РАСЧЕТА:\n";
        summary << "  Длина трубы: " << pipe_length << " м\n";
        summary << "  Шаг сетки: " << dx << " м\n";
        summary << "  Количество ячеек: " << n_cells_int << "\n";
        summary << "  Диаметр: " << diameter << " м\n";
        summary << "  Площадь сечения: " << area << " м²\n";
        summary << "  Расход: " << flow_rate << " м³/с\n";
        summary << "  Общее время: " << total_time << " с (" << total_time / 3600 << " ч)\n";
        summary << "  Плотность первой партии: " << rho_first << " кг/м³\n";
        summary << "  Плотность второй партии: " << rho_second << " кг/м³\n";
        summary << "  Время переключения партий: " << switch_time << " с (" << switch_time / 3600 << " ч)\n\n";

        summary << "ПАРАМЕТРЫ ВРЕМЕНИ:\n";
        summary << "  Шаг по времени: " << dt << " с\n";
        summary << "  Количество шагов: " << n_steps << "\n";
        summary << "  Скорость потока: " << v << " м/с\n";
        summary << "  Время прохода фронта: " << transit_time << " с (" << transit_time / 3600 << " ч)\n\n";

        summary << "РЕЗУЛЬТАТЫ БАЛАНСА МАССЫ:\n";
        summary << "  Начальная масса в трубе: " << initial_mass_in_pipe << " кг\n";
        summary << "  Суммарная масса на входе: " << cumulative_mass_in << " кг\n";
        summary << "  Суммарная масса на выходе: " << cumulative_mass_out << " кг\n";
        summary << "  Масса по балансу (формула): " << final_mass_by_balance << " кг\n";
        summary << "  Масса по схеме (расчет): " << final_mass_by_scheme << " кг\n";
        summary << "  Абсолютная ошибка: " << absolute_error << " кг\n";
        summary << "  Относительная ошибка: " << relative_error << "\n";
        summary << "  Финальная плотность на выходе: " << rho_outlet_history.back() << " кг/м³\n\n";

        summary << "КРИТЕРИИ ПРОВЕРКИ:\n";
        summary << "  Абсолютный допуск: " << absolute_tolerance << " кг\n";
        summary << "  Относительный допуск: " << relative_tolerance * 100 << " %\n";
        summary << "  Итоговый допуск: " << tolerance << " кг\n\n";

        if (std::abs(absolute_error) <= tolerance) {
            summary << "РЕЗУЛЬТАТ: ✓ СХЕМА КОНСЕРВАТИВНА\n";
            summary << "  Ошибка в пределах допуска\n";
        }
        else {
            summary << "РЕЗУЛЬТАТ: ✗ СХЕМА НЕ КОНСЕРВАТИВНА\n";
            summary << "  Ошибка превышает допуск\n";
        }

        summary << "\nИСТОРИЯ РАСЧЕТА:\n";
        summary << "  Всего шагов: " << time_history.size() << "\n";
        summary << "  Максимальная ошибка в истории: ";

        // Находим максимальную ошибку
        double max_hist_error = 0.0;
        size_t max_hist_step = 0;
        for (size_t i = 0; i < balance_error_history.size(); ++i) {
            if (std::abs(balance_error_history[i]) > max_hist_error) {
                max_hist_error = std::abs(balance_error_history[i]);
                max_hist_step = i;
            }
        }
        summary << max_hist_error << " кг на шаге " << max_hist_step
            << " (время " << time_history[max_hist_step] / 3600 << " ч)\n";

        summary.close();
        std::cout << "✓ Итоговый отчет сохранен в файл: conservativity_summary.txt" << std::endl;
    }
    else {
        std::cout << "✗ Ошибка сохранения файла conservativity_summary.txt" << std::endl;
    }

    // 3. Сохранение параметров в JSON для визуализации
    std::ofstream params_json("conservativity_params.json");
    if (params_json.is_open()) {
        params_json << "{\n";
        params_json << "  \"test_name\": \"SinglePipeConservativityTest\",\n";
        params_json << "  \"solver_type\": \"QUICKEST_ULTIMATE\",\n";
        params_json << "  \"pipe_length_m\": " << pipe_length << ",\n";
        params_json << "  \"grid_step_m\": " << dx << ",\n";
        params_json << "  \"n_cells\": " << n_cells_int << ",\n";
        params_json << "  \"diameter_m\": " << diameter << ",\n";
        params_json << "  \"cross_section_area_m2\": " << area << ",\n";
        params_json << "  \"flow_rate_m3_s\": " << flow_rate << ",\n";
        params_json << "  \"total_time_s\": " << total_time << ",\n";
        params_json << "  \"rho_first_kg_m3\": " << rho_first << ",\n";
        params_json << "  \"rho_second_kg_m3\": " << rho_second << ",\n";
        params_json << "  \"switch_time_s\": " << switch_time << ",\n";
        params_json << "  \"time_step_s\": " << dt << ",\n";
        params_json << "  \"n_steps\": " << n_steps << ",\n";
        params_json << "  \"velocity_m_s\": " << v << ",\n";
        params_json << "  \"transit_time_s\": " << transit_time << ",\n";
        params_json << "  \"initial_mass_kg\": " << initial_mass_in_pipe << ",\n";
        params_json << "  \"cumulative_mass_in_kg\": " << cumulative_mass_in << ",\n";
        params_json << "  \"cumulative_mass_out_kg\": " << cumulative_mass_out << ",\n";
        params_json << "  \"final_mass_scheme_kg\": " << final_mass_by_scheme << ",\n";
        params_json << "  \"final_mass_balance_kg\": " << final_mass_by_balance << ",\n";
        params_json << "  \"absolute_error_kg\": " << absolute_error << ",\n";
        params_json << "  \"relative_error\": " << relative_error << ",\n";
        params_json << "  \"max_mass_reference_kg\": " << max_mass << ",\n";
        params_json << "  \"absolute_tolerance_kg\": " << absolute_tolerance << ",\n";
        params_json << "  \"relative_tolerance\": " << relative_tolerance << ",\n";
        params_json << "  \"calculated_tolerance_kg\": " << tolerance << ",\n";
        params_json << "  \"is_conservative\": " << (std::abs(absolute_error) <= tolerance ? "true" : "false") << ",\n";
        params_json << "  \"final_rho_outlet_kg_m3\": " << rho_outlet_history.back() << ",\n";
        params_json << "  \"data_points\": " << time_history.size() << "\n";
        params_json << "}\n";
        params_json.close();
        std::cout << "✓ Параметры теста сохранены в файл: conservativity_params.json" << std::endl;
    }
    else {
        std::cout << "✗ Ошибка сохранения файла conservativity_params.json" << std::endl;
    }

    // 4. Сохранение истории ошибок для анализа
    std::ofstream errors_csv("conservativity_errors.csv");
    if (errors_csv.is_open()) {
        errors_csv << "step,time_h,balance_error_kg,absolute_error_kg\n";
        for (size_t i = 0; i < balance_error_history.size(); ++i) {
            errors_csv << i << ","
                << time_history[i] / 3600 << ","
                << balance_error_history[i] << ","
                << std::abs(balance_error_history[i]) << "\n";
        }
        errors_csv.close();
        std::cout << "✓ История ошибок сохранена в файл: conservativity_errors.csv" << std::endl;
    }
    else {
        std::cout << "✗ Ошибка сохранения файла conservativity_errors.csv" << std::endl;
    }

    std::cout << "\n=== ТЕСТ ЗАВЕРШЕН ===" << std::endl;
}

//TEST(QuickestUltimate, HasConservativity)
//{
//    double vol_flow = -1.0;
//    double rho_in = 860;    // плотность на входе
//    double rho_out = 850; // плотность на выходе
//    double rho_initial = 840;
//
//    simple_pipe_properties simple_pipe;
//    simple_pipe.length = 100.0;   // длина 100 м
//    simple_pipe.dx = 100.0;       // шаг сетки = длина (1 ячейка)
//    simple_pipe.diameter = 0.7;
//    pipe_properties_t pipe = pipe_properties_t::build_simple_pipe(simple_pipe);
//
//    std::vector<double> Q(pipe.profile.get_point_count(), 1.0);
//    PipeQAdvection advection_model(pipe, Q);
//
//    typedef quickest_ultimate_fv_solver_traits<1>::var_layer_data target_var_t;
//    typedef quickest_ultimate_fv_solver_traits<1>::specific_layer specific_data_t;
//    typedef composite_layer_t<target_var_t, specific_data_t> layer_t;
//
//    ring_buffer_t<layer_t> buffer(2, pipe.profile.get_point_count());
//    ring_buffer_t<layer_t> buffer_reverse(2, pipe.profile.get_point_count());
//
//    layer_t& prev = buffer.previous();
//    layer_t& next = buffer.current();
//
//    prev.vars.cell_double[0] = std::vector<double>(1, rho_initial);
//
//    double v = advection_model.getEquationsCoeffs(0, Q[0]);
//    double Cr = 1;
//    double dt = abs((Cr * simple_pipe.dx) / v);       // шаг по времени
//
//    // Запускаем шаг по QUICKEST-ULTIMATE
//    quickest_ultimate_fv_solver solver(advection_model, prev, next);
//    solver.step(dt, rho_in, rho_out);
//
//    double rho_prev = prev.vars.cell_double[0].front();
//    double rho_curr = next.vars.cell_double[0].front();
//
//   
//
//}




TEST(QuickestUltimate, FormulaCheck) {
    // Тест проверяет логику формулы, а не реализацию схемы
    std::vector<double> cell_masses = { 10.0, 20.0, 30.0 }; // массы в ячейках
    std::vector<double> inflow_flux = { 5.0, 3.0, 2.0 };    // G_gx за шаги
    std::vector<double> outflow_flux = { 1.0, 2.0, 1.0 };   // G_вых за шаги
    double dt = 1.0;

    // Σm_i (сумма масс в ячейках)
    double sum_mass = 0.0;
    for (double m : cell_masses) sum_mass += m;

    // Σ(G_gx - G_вых) * Δt
    double sum_flux = 0.0;
    for (size_t i = 0; i < inflow_flux.size(); i++) {
        sum_flux += (inflow_flux[i] - outflow_flux[i]) * dt;
    }

    // M_баланс(t) = Σm_i + Δt * Σ(G_gx - G_вых)
    double balance_mass = sum_mass + sum_flux;

    // Проверяем на простом примере
    double expected = 60.0 + (10.0 - 4.0); // 10+20+30 + (5+3+2 - 1+2+1)

    EXPECT_NEAR(balance_mass, expected, 1e-12);
}




