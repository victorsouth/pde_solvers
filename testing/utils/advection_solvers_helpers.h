#pragma once

namespace upstream_solver_types {
;
    typedef pde_solvers::upstream_fv_solver_traits<1>::var_layer_data target_var_t;
    typedef pde_solvers::upstream_fv_solver_traits<1>::specific_layer specific_data_t;
    typedef pde_solvers::composite_layer_t<target_var_t, specific_data_t> layer_t;
}

namespace quick_solver_types {
;
    typedef pde_solvers::quick_fv_solver_traits<1>::var_layer_data target_var_t;
    typedef pde_solvers::quick_fv_solver_traits<1>::specific_layer specific_data_t;
    typedef pde_solvers::composite_layer_t<target_var_t, specific_data_t> layer_t;
}

namespace quickest_solver_types {
;
    typedef pde_solvers::quickest_fv_solver_traits<1>::var_layer_data target_var_t;
    typedef pde_solvers::quickest_fv_solver_traits<1>::specific_layer specific_data_t;
    typedef pde_solvers::composite_layer_t<target_var_t, specific_data_t> layer_t;
}

namespace quickest_ultimate_solver_types {
;
    typedef pde_solvers::quickest_ultimate_fv_solver_traits<1>::var_layer_data target_var_t;
    typedef pde_solvers::quickest_ultimate_fv_solver_traits<1>::specific_layer specific_data_t;
    typedef pde_solvers::composite_layer_t<target_var_t, specific_data_t> layer_t;

    using quickest_sequential_solver = pde_solvers::quickest_ultimate_fv_solver<
        pde_solvers::quickest_cell_compute_mode::sequential>;

    using quickest_parallel_solver = pde_solvers::quickest_ultimate_fv_solver<
        pde_solvers::quickest_cell_compute_mode::parallel>;
}

template <typename layer_t>
inline ring_buffer_t<layer_t> build_buffer(const pipe_properties_t& pipe, double initial_value)
{
    ring_buffer_t<layer_t> result(2, pipe.profile.get_point_count());
    layer_t& prev = result.previous();
    prev.vars.cell_double[0] = std::vector<double>(prev.vars.cell_double[0].size(), initial_value);
    return result;
}

template <typename layer_t>
inline ring_buffer_t<layer_t> build_linear_buffer(
    const pipe_properties_t& pipe, double value_left, double value_right)
{
    ring_buffer_t<layer_t> buffer(2, pipe.profile.get_point_count());
    const size_t cell_count = pipe.profile.get_point_count() - 1;
    const double span = value_right - value_left;
    const double step = span / static_cast<double>(cell_count - 1);
    buffer.previous().vars.cell_double[0] =
        pipe_profile_uniform::generate_uniform_grid(value_left, span, step);
    return buffer;
}

inline double calc_time_step_by_Courant(const PipeQAdvection& advection_model, double courant)
{
    const auto& x = advection_model.get_grid();
    double dx = x[1] - x[0];
    double v = advection_model.getEquationsCoeffs(0, 0);
    double result = courant * dx / abs(v);
    return result;
}

/// @brief Ожидаемое ускорение параллельного расчёта от числа потоков.
/// Консервативная оценка: четверть от степени двойки по числу потоков: 16 потоков → 4x, 8 → 2x.
inline double expected_parallel_speedup(unsigned int threads)
{
    unsigned int log2_threads = static_cast<unsigned int>(std::log2(threads));
    if (log2_threads < 2) {
        return 1.0;
    }
    return static_cast<double>(1u << (log2_threads - 2));
}

/// @brief Один шаг расчета адвекции с замером времени; возвращает профиль и длительность расчета
template <typename Solver, typename layer_t>
inline std::pair<std::vector<double>, double> run_timed_step(
    PipeQAdvection& advection_model, const pipe_properties_t& pipe,
    double dt, double value_left, double value_right, double value_new)
{
    auto buffer = build_linear_buffer<layer_t>(pipe, value_left, value_right);

    const auto t0 = std::chrono::steady_clock::now();
    Solver solver(advection_model, buffer.previous(), buffer.current());
    solver.step(dt, value_new, value_right);
    const auto t1 = std::chrono::steady_clock::now();

    std::pair<std::vector<double>, double> result{
        buffer.current().vars.cell_double[0],
        std::chrono::duration<double>(t1 - t0).count()
    };
    return result;
}
