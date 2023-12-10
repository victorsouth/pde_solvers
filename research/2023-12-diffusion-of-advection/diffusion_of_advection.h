#pragma once

/// @brief Тесты для солвера quickest_ultimate_fv_solver
class DiffusionOfAdvection : public ::testing::Test {
protected:
    // Профиль переменных
    typedef quickest_ultimate_fv_solver_traits<1>::var_layer_data target_var_t;
    typedef quickest_ultimate_fv_solver_traits<1>::specific_layer specific_data_t;

    // Слой: переменных Vars + сколько угодно служебных Specific
    typedef composite_layer_t<target_var_t, specific_data_t> layer_t;
protected: // здесь пока копим параметры эксперимента
    const double density_initial{ 850 };
    const double density_final{ 860 };
    const double volumetric_flow{ 0.5 };
    const double experiment_time = 350000; 

protected:
    /// @brief Параметры трубы
    pipe_properties_t pipe;

    oil_parameters_t oil;

    /// @brief Профиль расхода
    vector<double> Q;
    std::unique_ptr<PipeQAdvection> advection_model;
    std::unique_ptr<ring_buffer_t<layer_t>> buffer;
protected:

    /// @brief Подготовка к расчету для семейства тестов
    virtual void SetUp() override {
        simple_pipe_properties simple_pipe = simple_pipe_properties::sample_district();
        pipe = pipe_properties_t::build_simple_pipe(simple_pipe);
        
        pipe.wall.equivalent_roughness = 15e-5; // нужно только для моделирования физической диффузии!
        oil.viscosity.nominal_viscosity = 6e-7; // аналогично

        Q = vector<double>(pipe.profile.getPointCount(), volumetric_flow);
        advection_model = std::make_unique<PipeQAdvection>(pipe, Q);
        buffer = std::make_unique<ring_buffer_t<layer_t>>(2, pipe.profile.getPointCount());

        layer_t& prev = buffer->previous();
        prev.vars.cell_double[0] = vector<double>(prev.vars.cell_double[0].size(), density_initial);
    }
    
    /// @brief Расчет вытеснения первой партии нефтью с другой плотностью 
    /// Расчет выполняетс яметодом Quickest Ultimate для заданного числа Куранта
    /// Результат расчетных профилей записывается в файл вида "output Cr.csv"
    /// @param rho_initial Плотность исходной нефти (вытесняемой)
    /// @param rho_final Плотность вытесняющей нефти
    /// @param Cr Число Куранта
    /// @param T Период моделирования
    /// @param path Путь, куда пишется результат расчета
    void calc_quickest_with_cr(
        double rho_initial, double rho_final, double v,
        double Cr, double T, const string& path)
    {
        // Фиктивное граничное условие на выходе. Реально в эксперименте не задействуется
        double rho_out = 870; 

        const auto& x = advection_model->get_grid();
        double dx = x[1] - x[0];
        double dt_ideal = abs(dx / v);

        advection_model = std::make_unique<PipeQAdvection>(pipe, Q);
        buffer = std::make_unique<ring_buffer_t<layer_t>>(2, pipe.profile.getPointCount());

        layer_t& prev = buffer->previous();
        prev.vars.cell_double[0] = vector<double>(prev.vars.cell_double[0].size(), rho_initial);

        double t = 0; // текущее время

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
            solver.step(dt, rho_final, rho_out);

            layer_t& next = buffer->current();
            next.vars.print(t, output);


            buffer->advance(+1);

        }
        output.flush();
        output.close();
    }

    void calc_physical_diffusion_model(
        double rho_initial, double rho_final,
        const vector<double>& t, double v, const string& path)
    {
        // точность расчета интеграла (шаг в секундах)
        // и одновременно период дискретизации для задания входных граничных условий
        double delta_t = 1;

        double t_change = 60;
        size_t n_change = static_cast<size_t>(t_change / delta_t + 0.5) + 1;
        size_t input_size = static_cast<size_t>(t.back() / delta_t + 0.5);
        vector<double> input = 
            diffusion_transport_solver::create_boundary(rho_initial, rho_final, 
                input_size, n_change, n_change);

        diffusion_transport_solver solver(pipe, oil);
        vector<double> density_output = solver.solve(t, delta_t, input, v, true);

        std::stringstream filename;
        filename << path << "output physical.csv";
        std::ofstream output(filename.str());

        for (size_t index = 0; index < density_output.size(); ++index) {
            output << t[index] << ";" << density_output[index] << std::endl;
        }
        output.close();
        output.flush();
    }

    vector<double> build_density_out_time_of_interest(double dt_out)
    {
        // чтобы не считать лишнего, делаем расчет для самых интересных двух часов
        size_t N = (3600 * 2) / dt_out; 
        vector<double> t(N);
        for (size_t i = 0; i < N; i++)
        {
            t[i] = (i + 1 + 4800) * dt_out; //Конец трубы
        }
        return t;
    }
};


/// @brief Расчет, демонстрирующий, что QUICKEST-ULTIMATE имеет диффузию 
/// не хуже физической
/// Строятся графики QUICKEST-ULTIMATE для разных Cr
/// Для сравнения строится также график с физической диффузией
/// !! Возможно, здесь же построить границы для физической диффузии 
/// по модели смесеобразования
TEST_F(DiffusionOfAdvection, CompareQuickestDiffusion)
{
    string path = prepare_research_folder();

    double v = advection_model->getEquationsCoeffs(0, 0);

    for (double Cr = 0.05; Cr < 1.01; Cr += 0.05) {
        calc_quickest_with_cr(density_initial, density_final, v,
            Cr, experiment_time, path);
    }

    // Задание массива моментов времени для расчета выходного параметра
    double dt_out = 60;
    vector<double> t_out = build_density_out_time_of_interest(dt_out);
    calc_physical_diffusion_model(density_initial, density_final, t_out, v, path);


}
