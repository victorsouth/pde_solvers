#pragma once

/// @file courant_spatial_profile.h
/// @brief Исследование влияния числа Куранта на пространственный профиль плотности
///
/// Исследуется, как число Куранта влияет на форму профиля плотности по длине трубы
/// при расчёте движения партий методами адвекции: Upstream Differencing, QUICK,
/// QUICKEST, QUICKEST-ULTIMATE и методом характеристик (MOC).
///
/// В отличие от исследования diffusion_of_advection (2023-12), которое анализирует
/// изменение плотности на выходе трубы во времени и сравнивает с физической моделью
/// диффузии, данное исследование фиксирует пространственное распределение плотности
/// по всей длине трубы в выбранный момент времени. Это позволяет визуально оценить
/// форму фронта раздела партий и степень его размытия для каждого солвера.

template <typename Layer, typename StepFn>
void courant_profile_time_loop(pde_solvers::ring_buffer_t<Layer>& buffer, double dt, size_t N,
                                std::ostream& output, StepFn step_fn)
{
    double t = 0;
    for (size_t index = 0; index < N; ++index) {
        if (index == 0) {
            buffer.previous().vars.print(t, output);
        }
        t += dt;
        step_fn(dt);
        buffer.current().vars.print(t, output);
        buffer.advance(+1);
    }
}

/// @brief Выводит пространственный профиль плотности upstream_fv_solver в файл
/// по диапазону чисел Куранта
TEST(CourantSpatialProfile, UpstreamDifferencing)
{
    using namespace upstream_solver_types;

    pipe_properties_t pipe = pipe_properties_t::build_simple_pipe(
        {.length = 700e3, .dx = 100, .diameter = 0.514});
    double rho_initial = 850, rho_in = 860, rho_out = 870, T = 350000, flow = 0.5;

    std::vector<double> Q(pipe.profile.get_point_count(), flow);
    PipeQAdvection advection_model(pipe, Q);
    double dt_ideal = calc_time_step_by_Courant(advection_model, 1.0);
    std::string path = prepare_research_folder();

    for (double Cr = 0.5; Cr < 0.51; Cr += 0.05) {
        ring_buffer_t<layer_t> buffer = build_buffer<layer_t>(pipe, rho_initial);
        double dt = Cr * dt_ideal;

        std::ofstream output(path + "output Cr=" + std::to_string(Cr) + ".csv");
        courant_profile_time_loop(buffer, dt, (size_t)(T / dt), output,
            [&](double dt) {
                upstream_fv_solver(advection_model, buffer).step(dt, rho_in, rho_out);
            });
    }
}

/// @brief Выводит пространственный профиль плотности quick_fv_solver в файл
/// по диапазону чисел Куранта
TEST(CourantSpatialProfile, QuickFvSolver)
{
    using namespace quick_solver_types;

    pipe_properties_t pipe = pipe_properties_t::build_simple_pipe(
        {.length = 700e3, .dx = 100, .diameter = 0.514});
    double rho_initial = 850, rho_in = 860, rho_out = 870, T = 350000, flow = 0.5;

    std::vector<double> Q(pipe.profile.get_point_count(), flow);
    PipeQAdvection advection_model(pipe, Q);
    double dt_ideal = calc_time_step_by_Courant(advection_model, 1.0);
    std::string path = prepare_research_folder();

    for (double Cr = 0.5; Cr < 0.51; Cr += 0.05) {
        ring_buffer_t<layer_t> buffer = build_buffer<layer_t>(pipe, rho_initial);
        double dt = Cr * dt_ideal;

        std::ofstream output(path + "output Cr=" + std::to_string(Cr) + ".csv");
        courant_profile_time_loop(buffer, dt, (size_t)(T / dt), output,
            [&](double dt) {
                quick_fv_solver(advection_model, buffer).step(dt, rho_in, rho_out);
            });
    }
}

/// @brief Выводит пространственный профиль плотности quickest_fv_solver в файл
/// по диапазону чисел Куранта
TEST(CourantSpatialProfile, QuickestFvSolver)
{
    using namespace quickest_solver_types;

    pipe_properties_t pipe = pipe_properties_t::build_simple_pipe(
        {.length = 700e3, .dx = 100, .diameter = 0.514});
    double rho_initial = 850, rho_in = 860, rho_out = 870, T = 350000, flow = 0.5;

    std::vector<double> Q(pipe.profile.get_point_count(), flow);
    PipeQAdvection advection_model(pipe, Q);
    double dt_ideal = calc_time_step_by_Courant(advection_model, 1.0);
    std::string path = prepare_research_folder();

    for (double Cr = 0.5; Cr < 0.51; Cr += 0.05) {
        ring_buffer_t<layer_t> buffer = build_buffer<layer_t>(pipe, rho_initial);
        double dt = Cr * dt_ideal;

        std::ofstream output(path + "output Cr=" + std::to_string(Cr) + ".csv");
        courant_profile_time_loop(buffer, dt, (size_t)(T / dt), output,
            [&](double dt) {
                quickest_fv_solver(advection_model, buffer).step(dt, rho_in, rho_out);
            });
    }
}

/// @brief Выводит пространственный профиль плотности quickest_ultimate_fv_solver в файл
/// по диапазону чисел Куранта
TEST(CourantSpatialProfile, QuickestUltimate)
{
    using namespace quickest_ultimate_solver_types;

    pipe_properties_t pipe = pipe_properties_t::build_simple_pipe(
        {.length = 700e3, .dx = 100, .diameter = 0.514});
    double rho_initial = 850, rho_in = 860, rho_out = 870, T = 350000, flow = 0.5;

    std::vector<double> Q(pipe.profile.get_point_count(), flow);
    PipeQAdvection advection_model(pipe, Q);
    double dt_ideal = calc_time_step_by_Courant(advection_model, 1.0);
    std::string path = prepare_research_folder();

    for (double Cr = 0.05; Cr < 1.01; Cr += 0.05) {
        ring_buffer_t<layer_t> buffer = build_buffer<layer_t>(pipe, rho_initial);
        double dt = Cr * dt_ideal;

        std::ofstream output(path + "output Cr=" + std::to_string(Cr) + ".csv");
        courant_profile_time_loop(buffer, dt, (size_t)(T / dt), output,
            [&](double dt) {
                quickest_sequential_solver(advection_model, buffer).step(dt, rho_in, rho_out);
            });
    }
}

/// @brief Выводит пространственный профиль плотности moc_solver в файл
/// по диапазону чисел Куранта
TEST(CourantSpatialProfile, MocSolver)
{
    typedef pde_solvers::composite_layer_t<profile_collection_t<1>,
        moc_solver<1>::specific_layer> single_var_moc_t;

    pipe_properties_t pipe = pipe_properties_t::build_simple_pipe(
        {.length = 700e3, .dx = 100, .diameter = 0.514});
    double rho_initial = 850, rho_in = 860, rho_out = 870, T = 350000, flow = 0.5;

    std::vector<double> Q(pipe.profile.get_point_count(), flow);
    PipeQAdvection advection_model(pipe, Q);
    double dt_ideal = calc_time_step_by_Courant(advection_model, 1.0);
    std::string path = prepare_research_folder();

    for (double Cr = 0.05; Cr < 0.51; Cr += 0.05) {
        ring_buffer_t<single_var_moc_t> buffer(2, pipe.profile.get_point_count());
        auto& rho_prev = buffer.previous().vars.point_double[0];
        rho_prev = std::vector<double>(rho_prev.size(), rho_initial);
        double dt = Cr * dt_ideal;

        std::ofstream output(path + "output Cr=" + std::to_string(Cr) + ".csv");
        courant_profile_time_loop(buffer, dt, (size_t)(T / dt), output,
            [&](double dt) {
                moc_solver<1>(advection_model, buffer.previous(), buffer.current())
                    .step_optional_boundaries(dt, rho_in, rho_out);
            });
    }
}
