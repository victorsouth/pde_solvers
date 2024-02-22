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
    // Момент времени когда область смеси наполовину своей длины пройдёт через конец трубопровода
    double physical_diffusion_center_time; 
protected:
    /// @brief Параметры трубы
    pipe_properties_t pipe;
    /// @brief Параметры нефти
    oil_parameters_t oil;

    /// @brief Профиль расхода
    vector<double> Q;
    /// @brief Уравнение адвекции
    std::unique_ptr<PipeQAdvection> advection_model;
    /// @brief Буфер, содержащий в себе слои
    std::unique_ptr<ring_buffer_t<layer_t>> buffer;
protected:

    /// @brief Подготовка к расчету для семейства тестов
    virtual void SetUp() override {
        // Инициализация модели трубы
        simple_pipe_properties simple_pipe = simple_pipe_properties::sample_district();
        pipe = pipe_properties_t::build_simple_pipe(simple_pipe);
        
        pipe.wall.equivalent_roughness = 15e-5; // нужно только для моделирования физической диффузии!
        oil.viscosity.nominal_viscosity = 6e-7; // и для расчёта границ физической диффузии

        // Инициализация профиля расхода
        Q = vector<double>(pipe.profile.getPointCount(), volumetric_flow);
        advection_model = std::make_unique<PipeQAdvection>(pipe, Q);
        buffer = std::make_unique<ring_buffer_t<layer_t>>(2, pipe.profile.getPointCount());
    }
    
    /// @brief Формирует имя файл для результатов исследования разных численных метов
    /// @tparam Solver Класс солвера
    /// @param path Путь, в котором формируется файл
    /// @param Cr Число куранта
    /// @return Путь к файлу в заивисимости от указанного класса солвера
    template <typename Solver>
    static std::string get_courant_research_filename(const string& path, double Cr)
    {
        std::stringstream filename;
        if constexpr (std::is_same<Solver, quickest_ultimate_fv_solver>::value) {
            filename << path << "output ultimate Cr=" << Cr << ".csv";
        }
        else if constexpr (std::is_same<Solver, upstream_fv_solver>::value) {
            filename << path << "output upstream Cr=" << Cr << ".csv";
        }
        else if constexpr (std::is_same<Solver, quickest_fv_solver>::value) {
            filename << path << "output quickest Cr = " << Cr << ".csv";
        }
        else if constexpr (std::is_same<Solver, quick_fv_solver>::value) {
            filename << path << "output quick Cr = " << Cr << ".csv";
        }
        else if constexpr (std::is_same<Solver, moc_solver<1>>::value)
        {
            filename << path << "output MOC Cr=" << Cr << ".csv";
        }
        else {
            throw std::runtime_error("Please specify filename for solver type");
        }
        return filename.str();
    }

    /// @brief Расчет вытеснения первой партии нефтью с другой плотностью 
    /// Расчет выполняется методом QUICKEST или QUICKEST-ULTIMATE для заданного числа Куранта
    /// Результат расчетных профилей записывается в файл вида "output Cr.csv"
    /// @tparam Solver Шаблон солвера для расчёта диффузии
    /// @param rho_initial Плотность исходной нефти (вытесняемой)
    /// @param rho_final Плотность вытесняющей нефти
    /// @param v Скорость движения нефти в трубопроводе
    /// @param Cr Число Куранта
    /// @param T Период моделирования
    /// @param path Путь, куда пишется результат расчета
    template <typename Solver>
    void calc_quickest_with_cr(double rho_initial, double rho_final, double v,
        double Cr, double T, const string& path)
    {
        // Фиктивное граничное условие на выходе. Реально в эксперименте не задействуется
        double rho_out = 870;

        // Задаём исходные значения
        layer_t& prev = buffer->previous();
        prev.vars.cell_double[0] = vector<double>(prev.vars.cell_double[0].size(), density_initial);

        const auto& x = advection_model->get_grid();
        double dx = x[1] - x[0];
        double dt_ideal = abs(dx / v); // Расчёт идеального шага по времени

        double t = 0; // текущее время

        double dt = Cr * dt_ideal; // время в долях от Куранта

        string filename = get_courant_research_filename<Solver>(path, Cr);

        std::ofstream output(filename); // Открытие файла для записи
        output << "time;Density" << std::endl;

        size_t N = static_cast<int>(T / dt);
        for (size_t index = 0; index < N; ++index) {
            if (index == 0) {
                layer_t& prev = buffer->previous();
                // Вывод значения плотности в конце трубы в начале моделирования
                output << t << ';' << prev.vars.cell_double[0].back() << std::endl;
                
            }
            t += dt;
            Solver solver(*advection_model, *buffer);
            solver.step(dt, rho_final, rho_out); // Шаг расчёта
            layer_t& next = buffer->current();
            // Вывод значения плотности в конце трубы на каждом шаге моделирования
            output << t << ';' << next.vars.cell_double[0].back() << std::endl;
            buffer->advance(+1);
        }
        output.flush();
        output.close();
    }

    /// @brief Осуществляет поиск момента времени
    /// когда область смеси наполовину своей длины пройдёт через конец трубопровода
    /// @param time Временной ряд, по которому была рассчитана модель физической диффузии
    /// @param parameter Рассчитанная модель физической диффузии на выходе трубопровода
    double calc_physical_diffusion_center(const vector<double>& time, const vector<double>& parameter) const
    {
        double eps = 0.0001; // Точность для определения начала и конца области диффузии
        double start_diff = 0; // момент времени начала прохождения области смеси через конец трубопровода
        double end_diff; // момент времени конца прохождения области смеси через конец трубопровода
        
        for (size_t index = 0; index < parameter.size(); ++index)
        {
            if ((abs(parameter[index] - density_initial) < eps) && (abs(parameter[index + 1] - density_initial) >= eps))
            {
                start_diff = time[index];
                break;
            }
        }

        for (size_t index = parameter.size() - 1; index > 0; --index)
        {
            if ((abs(parameter[index] - density_final) < eps) && (abs(parameter[index - 1] - density_final) >= eps))
            {
                end_diff = time[index];
                break;
            }
        }

        if (start_diff == 0) {
            // Задаем значение start_diff по умолчанию чтобы не возвращался мусор
            start_diff = time[0];
        }
        
        return (start_diff + end_diff) / 2;
    }

    /// @brief Получения временного диапазона, в котором область смеси двух партий нефти
    /// проходит через конец трубопровода
    /// @param dt_out Шаг по времени
    /// @return возвращает временной ряд интересующего нас периода
    vector<double> build_density_out_time_of_interest(double dt_out)
    {
        // чтобы не считать лишнего, делаем расчет для самых интересных двух часов
        size_t N = static_cast<size_t>((3600 * 2) / dt_out);
        vector<double> t(N);
        for (size_t i = 0; i < N; i++)
        {
            t[i] = (i + 1 + 4800) * dt_out; //Конец трубы
        }
        return t;
    }

    /// @brief Расчёт модели физической диффузии
    /// В результате получаем изменение плотности на выходе трубопровода за два часа,
    /// за которые область смешения проходит через конец трубы
    /// Результат записывается в файл output "output physical.csv"
    /// @param rho_initial Плотность исходной нефти (вытесняемой)
    /// @param rho_final Плотность вытесняющей нефти
    /// @param t Выбранный отрезок времени - два часа
    /// @param v Скорость движения нефти в трубопроводе
    /// @param path Путь, куда пишется результат расчета
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

        physical_diffusion_center_time = calc_physical_diffusion_center(t, density_output); // Поиск временного центра области смешения
    }

    /// @brief Расчёт длины области смеси двух партий нефти
    /// @param speed Скорость движения среды в трубопроводе
    /// @return Возвращает длину области смеси в метрах
    double calc_physical_diffusion_length(const double speed)
    {
        double L = pipe.profile.getLength();
        double S = M_PI * pow(pipe.wall.diameter, 2) / 4; 
        double kc_v = diffusion_transport_solver::calc_diffusion_coefficient(pipe, oil, speed) / speed;
        double vc = 6.58 * S * sqrt(kc_v) * sqrt(L); // Расчёт объёма смеси, формула из учебника Лурье 2012, с. 409
        return vc / S;
    }

    /// @brief Расчёт временных границ области смеси двух партий нефти
    /// Моменты, когда область смешения начинает и заканчивает проходить через конец трубопровода
    /// записываются в файл "physical_diffusion.txt"
    /// @param speed Скорость движения среды в трубопроводе
    /// @param path Путь, куда сохраняются результы
    void calc_physical_diffusion_boundaries(const double& speed, const string& path)
    {
        double boundaries[2];
        double diff_length = calc_physical_diffusion_length(speed); // расчёт длины области смеси

        boundaries[0] = physical_diffusion_center_time - diff_length / (speed * 2);
        boundaries[1] = physical_diffusion_center_time + diff_length / (speed * 2);
        std::stringstream filename;
        filename << path << "physical_diffusion.txt";
        std::ofstream output(filename.str());

        output << boundaries[0] << "," << boundaries[1] << std::endl;
        
        output.close();
        output.flush();
    }

    /// @brief Расчет вытеснения первой партии нефтью с другой плотностью 
    /// Расчет выполняется методом характеристик для заданного числа Куранта
    /// Результат расчетных профилей записывается в файл вида "output Cr.csv"
    /// @tparam Solver Шаблон солвера для расчёта диффузии
    /// @param rho_initial Плотность исходной нефти (вытесняемой)
    /// @param rho_final Плотность вытесняющей нефти
    /// @param v Скорость движения нефти в трубопроводе
    /// @param Cr Число Куранта
    /// @param T Период моделирования
    /// @param path Путь, куда пишется результат расчета
    template <typename Solver>
    void calc_moc_with_cr(double rho_initial, double rho_final, double v,
        double Cr, double T, const string& path)
    {
        // Фиктивное граничное условие на выходе. Реально в эксперименте не задействуется
        double rho_out = 870;

        // Одна переменная, и структуры метода характеристик для нее
        typedef composite_layer_t<profile_collection_t<1>,
            moc_solver<1>::specific_layer> single_var_moc_t;
        // буфер, содержащий в себе слои
        ring_buffer_t<single_var_moc_t> buffer(2, pipe.profile.getPointCount());
        // объявляем переменную
        single_var_moc_t& prev = buffer.previous();
        // задаем исходные значения
        prev.vars.point_double[0] = vector<double>(prev.vars.point_double[0].size(), density_initial);
        // профиль расхода
        vector<double> Q(pipe.profile.getPointCount(), 0.5);
        // уравнение адвекции
        PipeQAdvection advection_model(pipe, Q);

        const auto& x = advection_model.get_grid();
        double dx = x[1] - x[0];
        double dt_ideal = abs(dx / v);

        double t = 0; // текущее время

        double dt = Cr * dt_ideal; // время в долях от Куранта

        string filename = get_courant_research_filename<Solver>(path, Cr);

        std::ofstream output(filename); // Открытие файла для записи
        output << "time;Density" << std::endl;

        size_t N = static_cast<int>(T / dt);
        for (size_t index = 0; index < N; ++index) {
            if (index == 0) {
                single_var_moc_t& prev = buffer.previous();
                // Вывод значения плотности в конце трубы на каждом шаге моделирования
                output << t << ';' << prev.vars.point_double[0].back() << std::endl;
            }

            t += dt;

            Solver solver(advection_model, buffer.previous(), buffer.current());
            solver.step_optional_boundaries(dt, rho_final, rho_out); // Шаг расчёта

            single_var_moc_t& next = buffer.current();
            // Вывод значения плотности в конце трубы на каждом шаге моделирования
            output << t << ';' << next.vars.point_double[0].back() << std::endl;

            buffer.advance(+1);

        }
        output.flush();
        output.close();
    }
};


/// @brief Расчет, демонстрирующий, что QUICKEST-ULTIMATE имеет диффузию 
/// не хуже физической
/// Строятся графики QUICKEST-ULTIMATE для разных Cr
/// Для сравнения строится также график с физической диффузией
/// Производится расчёт границ для физической диффузии 
/// по модели смесеобразования
TEST_F(DiffusionOfAdvection, CompareQuickestDiffusion)
{
    // Создаём папку с результатами и получаем путь к ней
    string path = prepare_research_folder();

    // Получаем значение скорости для эксперимента
    double v = advection_model->getEquationsCoeffs(0, 0);

    // Производим моделирование движения партий методом QUICKEST-ULTIMATE
    // для разных чисел Cr
    for (double Cr = 0.05; Cr < 1.01; Cr += 0.05) {
        calc_quickest_with_cr<quickest_ultimate_fv_solver>(density_initial, density_final, v,
            Cr, experiment_time, path);
    }

    // Задание массива моментов времени для расчета выходного параметра
    double dt_out = 60;
    vector<double> t_out = build_density_out_time_of_interest(dt_out);
    // Расчёт физической диффузии
    calc_physical_diffusion_model(density_initial, density_final, t_out, v, path);

    // Расчёт временных границ области смеси
    calc_physical_diffusion_boundaries(v, path);
}

/// @brief Расчет, демонстрирующий, что QUICKEST-ULTIMATE лучше описывает диффузию 
/// по сравенению с методом QUICKEST
/// Строятся графики QUICKEST-ULTIMATE, QUICKEST для Cr = 0.5 
/// и график Аналитического решения (QUICKEST-ULTIMATE для Cr = 1)
TEST_F(DiffusionOfAdvection, CompareQuickestAndQuickestUltimateDiffusion)
{
    // Создаём папку с результатами и получаем путь к ней
    string path = prepare_research_folder();

    // Получаем значение скорости для эксперимента
    double v = advection_model->getEquationsCoeffs(0, 0);

    // Строим для Cr = 1
    double Cr = 0.5;

    // Расчёт методом QUICKEST-UlTIMATE для Cr = 0.5
    calc_quickest_with_cr<quickest_ultimate_fv_solver>(density_initial, density_final, v,
        Cr, experiment_time, path);

    // Расчёт методом QUICKEST для Cr = 0.5
    calc_quickest_with_cr<quickest_fv_solver>(density_initial, density_final, v,
        Cr, experiment_time, path);

    // Строим для Cr = 1
    Cr = 1;

    // Расчёт методом QUICKEST-UlTIMATE для Cr = 1
    calc_quickest_with_cr<quickest_ultimate_fv_solver>(density_initial, density_final, v,
        Cr, experiment_time, path);
}

/// @brief Расчет, демонстрирующий, что Upstream Differencing и метод характеристик
/// имеют численную диффузию хуже физической
/// Строятся графики Upstream Differencing, метод характеристик для Cr = 0.5 
/// и график Аналитического решения (Upstream Differencing для Cr = 1)
/// Так же строится график физической диффузии
TEST_F(DiffusionOfAdvection, CompareUpstreamDiffAndMOC)
{
    // Создаём папку с результатами и получаем путь к ней
    string path = prepare_research_folder();

    // Получаем значение скорости для эксперимента
    double v = advection_model->getEquationsCoeffs(0, 0);

    // Строим для Cr = 1
    double Cr = 0.5;

    // Расчёт методом Upstream Differencing для Cr = 0.5
    calc_quickest_with_cr<upstream_fv_solver>(density_initial, density_final, v,
        Cr, experiment_time, path);

    // Расчёт методом характеристик для Cr = 0.5
    calc_moc_with_cr<moc_solver<1>>(density_initial, density_final, v,
        Cr, experiment_time, path);

    // Строим для Cr = 1
    Cr = 1;

    // Расчёт методом Upstream Differencing для Cr = 1, аналитическое решение
    calc_quickest_with_cr<upstream_fv_solver>(density_initial, density_final, v,
        Cr, experiment_time, path);

    // Задание массива моментов времени для расчета выходного параметра
    double dt_out = 60;
    vector<double> t_out = build_density_out_time_of_interest(dt_out);
    // Расчёт физической диффузии
    calc_physical_diffusion_model(density_initial, density_final, t_out, v, path);

    // Расчёт временных границ области смеси
    calc_physical_diffusion_boundaries(v, path);
}