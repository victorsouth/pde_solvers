#pragma once

#include <iomanip>
#include <tbb/parallel_for.h>
#include <tbb/parallel_invoke.h>
#include <tbb/task_arena.h>
#include <tbb/task_group.h>
#include <tbb/partitioner.h>

/// @brief Альтернативный эндогенный слой с собственным specific на каждый из 4 параметров.
/// В отличие от iso_nonbaro_improver_pipe_endogenious_layer_t (общий specific),
/// здесь гонки данных между параллельными задачами исключены полностью.
namespace separate_specific_layer {

using namespace pde_solvers;

/// @brief Профиль трубопровода: плотность, присадка и их достоверности — у каждого свой буфер потоков
struct endogenous_layer_t {
    confident_layer_t density_std;
    confident_layer_t improver_concentration;
    quickest_ultimate_fv_solver_traits<1>::specific_layer specific_density;
    quickest_ultimate_fv_solver_traits<1>::specific_layer specific_density_conf;
    quickest_ultimate_fv_solver_traits<1>::specific_layer specific_improver;
    quickest_ultimate_fv_solver_traits<1>::specific_layer specific_improver_conf;

    /// @brief Инициализация профиля трубы: равномерная плотность и нулевая присадка
    endogenous_layer_t(size_t point_count)
        : density_std(point_count, 850.0)
        , improver_concentration(point_count, 0.0)
        , specific_density(point_count)
        , specific_density_conf(point_count)
        , specific_improver(point_count)
        , specific_improver_conf(point_count)
    {}

    /// @brief Доступ к профилю плотности для шага адвекции
    static quickest_ultimate_fv_wrapper<1> get_density_wrapper(endogenous_layer_t& layer)
        { return quickest_ultimate_fv_wrapper<1>(layer.density_std.value, layer.specific_density); }
    /// @brief Доступ к достоверности плотности
    static quickest_ultimate_fv_wrapper<1> get_density_confidence_wrapper(endogenous_layer_t& layer)
        { return quickest_ultimate_fv_wrapper<1>(layer.density_std.confidence, layer.specific_density_conf); }
    /// @brief Доступ к концентрации присадки
    static quickest_ultimate_fv_wrapper<1> get_improver_wrapper(endogenous_layer_t& layer)
        { return quickest_ultimate_fv_wrapper<1>(layer.improver_concentration.value, layer.specific_improver); }
    /// @brief Доступ к достоверности концентрации присадки
    static quickest_ultimate_fv_wrapper<1> get_improver_confidence_wrapper(endogenous_layer_t& layer)
        { return quickest_ultimate_fv_wrapper<1>(layer.improver_concentration.confidence, layer.specific_improver_conf); }
};

} // namespace separate_specific_layer

/// @brief Вспомогательные типы и функции для бенчмарков параллельного расчёта QUICKEST
namespace quickest_parallel_bench {

using namespace pde_solvers;

/// @brief Число расчётных ячеек вдоль трубы
constexpr size_t CELL_COUNT  = 500;
/// @brief Число временных шагов в одном прогоне бенчмарка
constexpr size_t STEP_COUNT  = 1'000'000;
/// @brief Длина трубопровода, м
constexpr double PIPE_LENGTH = 100'000.0;
/// @brief Шаг сетки, м
constexpr double PIPE_DX     = PIPE_LENGTH / CELL_COUNT;
/// @brief Внутренний диаметр трубы, м
constexpr double PIPE_DIAM   = 0.514;
/// @brief Объёмный расход по трубе, м³/с
constexpr double FLOW        = 0.5;
/// @brief Начальная плотность нефти в трубе, кг/м³
constexpr double RHO_INIT    = 850.0;
/// @brief Плотность на входе (граничное условие), кг/м³
constexpr double RHO_IN      = 860.0;
/// @brief Достоверность плотности на входе
constexpr double CONF_IN     = 1.0;
/// @brief Концентрация присадки на входе
constexpr double IMPR_IN     = 0.1;
/// @brief Достоверность присадки на входе
constexpr double IMPR_CONF_IN = 0.9;

/// @brief Шаг адвекции одного параметра через обёртку из ring_buffer_t::get_buffer_wrapper
template <quickest_cell_compute_mode Mode, typename BufferType, typename Getter>
void do_step(PipeQAdvection& model, double dt, double bv, BufferType& buffer, Getter getter)
{
    auto wrap = buffer.get_buffer_wrapper(getter);
    quickest_ultimate_fv_solver<Mode> solver(model, wrap);
    solver.step(dt, bv, bv);
}

/// @brief Печать времени параллельного прогона и ускорения относительно последовательного
inline void print_result(const char* label, double seq_time, double par_time)
{
    double speedup = seq_time / par_time;
    std::cout << "  [" << label << "] "
        << std::fixed << std::setprecision(0)
        << par_time * 1000.0 << " ms"
        << "  (vs seq " << seq_time * 1000.0 << " ms"
        << ", speedup " << std::setprecision(2) << speedup << "x)\n";
}

/// @brief Печать эталонного времени последовательного расчёта
inline void print_seq(const char* label, double seq_time)
{
    std::cout << "  [" << label << "] "
        << std::fixed << std::setprecision(0)
        << seq_time * 1000.0 << " ms  (baseline)\n";
}

} // namespace quickest_parallel_bench

/// @brief Бенчмарк 1 — параллелизм по ячейкам одного профиля (500 ячеек);
/// quickest_ultimate_fv_solver не хранит состояния между шагами, стоимость конструктора нулевая

/// @brief Однопоточный расчёт одного параметра (плотность) — эталон для всех вариантов bench1.
TEST(TBB_Bench1_CellParallelism, DISABLED_Sequential_Baseline)
{
    using namespace quickest_parallel_bench;
    using namespace quickest_ultimate_solver_types;

    pipe_properties_t pipe_solo = pipe_properties_t::build_simple_pipe(
        { .length = PIPE_LENGTH, .dx = PIPE_DX, .diameter = PIPE_DIAM });
    std::vector<double> flow_rate_solo(pipe_solo.profile.get_point_count(), FLOW);
    PipeQAdvection model_solo(pipe_solo, flow_rate_solo);
    const double dt = calc_time_step_by_Courant(model_solo, 1.0);
    auto buffer = build_linear_buffer<layer_t>(pipe_solo, RHO_INIT, RHO_IN);

    quickest_sequential_solver solver(model_solo, buffer.previous(), buffer.current());
    const auto solo_start = std::chrono::steady_clock::now();
    for (size_t step = 0; step < STEP_COUNT; ++step) {
        solver.step(dt, RHO_IN, RHO_INIT);
        buffer.advance(+1);
    }
    print_seq("sequential", std::chrono::duration<double>(std::chrono::steady_clock::now() - solo_start).count());
}

/// @brief Параллелизм по ячейкам с auto_partitioner: TBB сам делит диапазон ячеек.
/// Два варианта: солвер живёт вне цикла (pde_solvers) и создаётся заново на каждом шаге (graph_solvers).
TEST(TBB_Bench1_CellParallelism, DISABLED_Parallel_AutoPartitioner)
{
    using namespace quickest_parallel_bench;
    using namespace quickest_ultimate_solver_types;

    pipe_properties_t pipe_seq = pipe_properties_t::build_simple_pipe(
        { .length = PIPE_LENGTH, .dx = PIPE_DX, .diameter = PIPE_DIAM });
    std::vector<double> flow_rate_seq(pipe_seq.profile.get_point_count(), FLOW);
    PipeQAdvection model_seq(pipe_seq, flow_rate_seq);
    const double dt = calc_time_step_by_Courant(model_seq, 1.0);
    auto buffer_seq = build_linear_buffer<layer_t>(pipe_seq, RHO_INIT, RHO_IN);
    quickest_sequential_solver solver_seq(model_seq, buffer_seq.previous(), buffer_seq.current());
    const auto seq_start = std::chrono::steady_clock::now();
    for (size_t step = 0; step < STEP_COUNT; ++step) {
        solver_seq.step(dt, RHO_IN, RHO_INIT);
        buffer_seq.advance(+1);
    }
    const double seq_elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - seq_start).count();

    /// @brief Типовой сценарий pde_solvers: один солвер, множество step()
    pipe_properties_t pipe_par = pipe_properties_t::build_simple_pipe(
        { .length = PIPE_LENGTH, .dx = PIPE_DX, .diameter = PIPE_DIAM });
    std::vector<double> flow_rate_par(pipe_par.profile.get_point_count(), FLOW);
    PipeQAdvection model_par(pipe_par, flow_rate_par);
    auto buffer_par = build_linear_buffer<layer_t>(pipe_par, RHO_INIT, RHO_IN);
    quickest_parallel_solver solver_par(model_par, buffer_par.previous(), buffer_par.current());
    const auto par_start = std::chrono::steady_clock::now();
    for (size_t step = 0; step < STEP_COUNT; ++step) {
        solver_par.step(dt, RHO_IN, RHO_INIT);
        buffer_par.advance(+1);
    }
    print_result("parallel(auto) pde_solvers",
        seq_elapsed, std::chrono::duration<double>(std::chrono::steady_clock::now() - par_start).count());

    /// @brief Типовой сценарий graph_solvers: солвер создаётся на каждом шаге
    pipe_properties_t pipe_graph = pipe_properties_t::build_simple_pipe(
        { .length = PIPE_LENGTH, .dx = PIPE_DX, .diameter = PIPE_DIAM });
    std::vector<double> flow_rate_graph(pipe_graph.profile.get_point_count(), FLOW);
    PipeQAdvection model_graph(pipe_graph, flow_rate_graph);
    auto buffer_graph = build_linear_buffer<layer_t>(pipe_graph, RHO_INIT, RHO_IN);
    const auto graph_start = std::chrono::steady_clock::now();
    for (size_t step = 0; step < STEP_COUNT; ++step) {
        quickest_parallel_solver(model_graph, buffer_graph.previous(), buffer_graph.current()).step(dt, RHO_IN, RHO_INIT);
        buffer_graph.advance(+1);
    }
    print_result("parallel(auto) graph_solvers",
        seq_elapsed, std::chrono::duration<double>(std::chrono::steady_clock::now() - graph_start).count());
}

/// @brief Параллелизм по ячейкам с static_partitioner: детерминированное деление без work stealing.
TEST(TBB_Bench1_CellParallelism, DISABLED_Parallel_StaticPartitioner)
{
    using namespace quickest_parallel_bench;
    using namespace quickest_ultimate_solver_types;

    pipe_properties_t pipe_seq = pipe_properties_t::build_simple_pipe(
        { .length = PIPE_LENGTH, .dx = PIPE_DX, .diameter = PIPE_DIAM });
    std::vector<double> flow_rate_seq(pipe_seq.profile.get_point_count(), FLOW);
    PipeQAdvection model_seq(pipe_seq, flow_rate_seq);
    const double dt = calc_time_step_by_Courant(model_seq, 1.0);
    auto buffer_seq = build_linear_buffer<layer_t>(pipe_seq, RHO_INIT, RHO_IN);
    quickest_sequential_solver solver_seq(model_seq, buffer_seq.previous(), buffer_seq.current());
    const auto seq_start = std::chrono::steady_clock::now();
    for (size_t step = 0; step < STEP_COUNT; ++step) {
        solver_seq.step(dt, RHO_IN, RHO_INIT);
        buffer_seq.advance(+1);
    }
    const double seq_elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - seq_start).count();

    pipe_properties_t pipe_par = pipe_properties_t::build_simple_pipe(
        { .length = PIPE_LENGTH, .dx = PIPE_DX, .diameter = PIPE_DIAM });
    std::vector<double> flow_rate_par(pipe_par.profile.get_point_count(), FLOW);
    PipeQAdvection model_par(pipe_par, flow_rate_par);
    auto buffer_par = build_linear_buffer<layer_t>(pipe_par, RHO_INIT, RHO_IN);
    quickest_parallel_solver solver_par(model_par, buffer_par.previous(), buffer_par.current());
    tbb::static_partitioner static_part;
    const auto par_start = std::chrono::steady_clock::now();
    for (size_t step = 0; step < STEP_COUNT; ++step) {
        solver_par.step(dt, RHO_IN, RHO_INIT, static_part);
        buffer_par.advance(+1);
    }
    print_result("parallel(static)", seq_elapsed, std::chrono::duration<double>(std::chrono::steady_clock::now() - par_start).count());
}

/// @brief Параллелизм по ячейкам с simple_partitioner: агрессивное дробление, каждая ячейка — отдельная задача.
/// Ожидается максимальный overhead планировщика.
TEST(TBB_Bench1_CellParallelism, DISABLED_Parallel_SimplePartitioner)
{
    using namespace quickest_parallel_bench;
    using namespace quickest_ultimate_solver_types;

    pipe_properties_t pipe_seq = pipe_properties_t::build_simple_pipe(
        { .length = PIPE_LENGTH, .dx = PIPE_DX, .diameter = PIPE_DIAM });
    std::vector<double> flow_rate_seq(pipe_seq.profile.get_point_count(), FLOW);
    PipeQAdvection model_seq(pipe_seq, flow_rate_seq);
    const double dt = calc_time_step_by_Courant(model_seq, 1.0);
    auto buffer_seq = build_linear_buffer<layer_t>(pipe_seq, RHO_INIT, RHO_IN);
    quickest_sequential_solver solver_seq(model_seq, buffer_seq.previous(), buffer_seq.current());
    const auto seq_start = std::chrono::steady_clock::now();
    for (size_t step = 0; step < STEP_COUNT; ++step) {
        solver_seq.step(dt, RHO_IN, RHO_INIT);
        buffer_seq.advance(+1);
    }
    const double seq_elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - seq_start).count();

    pipe_properties_t pipe_par = pipe_properties_t::build_simple_pipe(
        { .length = PIPE_LENGTH, .dx = PIPE_DX, .diameter = PIPE_DIAM });
    std::vector<double> flow_rate_par(pipe_par.profile.get_point_count(), FLOW);
    PipeQAdvection model_par(pipe_par, flow_rate_par);
    auto buffer_par = build_linear_buffer<layer_t>(pipe_par, RHO_INIT, RHO_IN);
    quickest_parallel_solver solver_par(model_par, buffer_par.previous(), buffer_par.current());
    tbb::simple_partitioner simple_part;
    const auto par_start = std::chrono::steady_clock::now();
    for (size_t step = 0; step < STEP_COUNT; ++step) {
        solver_par.step(dt, RHO_IN, RHO_INIT, simple_part);
        buffer_par.advance(+1);
    }
    print_result("parallel(simple)", seq_elapsed, std::chrono::duration<double>(std::chrono::steady_clock::now() - par_start).count());
}

/// @brief Параллелизм по ячейкам с affinity_partitioner: накапливает привязку поток→чанк между шагами.
/// Два отдельных партиционера для цикла F и цикла U внутри step().
TEST(TBB_Bench1_CellParallelism, DISABLED_Parallel_AffinityPartitioner)
{
    using namespace quickest_parallel_bench;
    using namespace quickest_ultimate_solver_types;

    pipe_properties_t pipe_seq = pipe_properties_t::build_simple_pipe(
        { .length = PIPE_LENGTH, .dx = PIPE_DX, .diameter = PIPE_DIAM });
    std::vector<double> flow_rate_seq(pipe_seq.profile.get_point_count(), FLOW);
    PipeQAdvection model_seq(pipe_seq, flow_rate_seq);
    const double dt = calc_time_step_by_Courant(model_seq, 1.0);
    auto buffer_seq = build_linear_buffer<layer_t>(pipe_seq, RHO_INIT, RHO_IN);
    quickest_sequential_solver solver_seq(model_seq, buffer_seq.previous(), buffer_seq.current());
    const auto seq_start = std::chrono::steady_clock::now();
    for (size_t step = 0; step < STEP_COUNT; ++step) {
        solver_seq.step(dt, RHO_IN, RHO_INIT);
        buffer_seq.advance(+1);
    }
    const double seq_elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - seq_start).count();

    pipe_properties_t pipe_par = pipe_properties_t::build_simple_pipe(
        { .length = PIPE_LENGTH, .dx = PIPE_DX, .diameter = PIPE_DIAM });
    std::vector<double> flow_rate_par(pipe_par.profile.get_point_count(), FLOW);
    PipeQAdvection model_par(pipe_par, flow_rate_par);
    auto buffer_par = build_linear_buffer<layer_t>(pipe_par, RHO_INIT, RHO_IN);
    quickest_parallel_solver solver_par(model_par, buffer_par.previous(), buffer_par.current());
    tbb::affinity_partitioner ap_flux, ap_update;
    const auto par_start = std::chrono::steady_clock::now();
    for (size_t step = 0; step < STEP_COUNT; ++step) {
        solver_par.step(dt, RHO_IN, RHO_INIT, ap_flux, ap_update);
        buffer_par.advance(+1);
    }
    print_result("parallel(affinity)", seq_elapsed, std::chrono::duration<double>(std::chrono::steady_clock::now() - par_start).count());
}

/// @brief Параллелизм по ячейкам в арене из 2 потоков: проверяет гипотезу о том, что 2 потока
/// держат весь профиль (4 КБ) в L1 без конкуренции за шину памяти.
TEST(TBB_Bench1_CellParallelism, DISABLED_Parallel_Arena2_StaticPartitioner)
{
    using namespace quickest_parallel_bench;
    using namespace quickest_ultimate_solver_types;

    pipe_properties_t pipe_seq = pipe_properties_t::build_simple_pipe(
        { .length = PIPE_LENGTH, .dx = PIPE_DX, .diameter = PIPE_DIAM });
    std::vector<double> flow_rate_seq(pipe_seq.profile.get_point_count(), FLOW);
    PipeQAdvection model_seq(pipe_seq, flow_rate_seq);
    const double dt = calc_time_step_by_Courant(model_seq, 1.0);
    auto buffer_seq = build_linear_buffer<layer_t>(pipe_seq, RHO_INIT, RHO_IN);
    quickest_sequential_solver solver_seq(model_seq, buffer_seq.previous(), buffer_seq.current());
    const auto seq_start = std::chrono::steady_clock::now();
    for (size_t step = 0; step < STEP_COUNT; ++step) {
        solver_seq.step(dt, RHO_IN, RHO_INIT);
        buffer_seq.advance(+1);
    }
    const double seq_elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - seq_start).count();

    pipe_properties_t pipe_par = pipe_properties_t::build_simple_pipe(
        { .length = PIPE_LENGTH, .dx = PIPE_DX, .diameter = PIPE_DIAM });
    std::vector<double> flow_rate_par(pipe_par.profile.get_point_count(), FLOW);
    PipeQAdvection model_par(pipe_par, flow_rate_par);
    auto buffer_par = build_linear_buffer<layer_t>(pipe_par, RHO_INIT, RHO_IN);
    quickest_parallel_solver solver_par(model_par, buffer_par.previous(), buffer_par.current());
    tbb::task_arena arena(2);
    tbb::static_partitioner static_part;
    const auto par_start = std::chrono::steady_clock::now();
    for (size_t step = 0; step < STEP_COUNT; ++step) {
        arena.execute([&] { solver_par.step(dt, RHO_IN, RHO_INIT, static_part); });
        buffer_par.advance(+1);
    }
    print_result("arena(2)+parallel(static)", seq_elapsed, std::chrono::duration<double>(std::chrono::steady_clock::now() - par_start).count());
}

/// @brief Параллелизм по ячейкам в арене из hardware_concurrency() потоков: верхняя граница
/// масштабирования по числу ядер на задаче из 500 ячеек.
TEST(TBB_Bench1_CellParallelism, DISABLED_Parallel_ArenaAll_StaticPartitioner)
{
    using namespace quickest_parallel_bench;
    using namespace quickest_ultimate_solver_types;

    pipe_properties_t pipe_seq = pipe_properties_t::build_simple_pipe(
        { .length = PIPE_LENGTH, .dx = PIPE_DX, .diameter = PIPE_DIAM });
    std::vector<double> flow_rate_seq(pipe_seq.profile.get_point_count(), FLOW);
    PipeQAdvection model_seq(pipe_seq, flow_rate_seq);
    const double dt = calc_time_step_by_Courant(model_seq, 1.0);
    auto buffer_seq = build_linear_buffer<layer_t>(pipe_seq, RHO_INIT, RHO_IN);
    quickest_sequential_solver solver_seq(model_seq, buffer_seq.previous(), buffer_seq.current());
    const auto seq_start = std::chrono::steady_clock::now();
    for (size_t step = 0; step < STEP_COUNT; ++step) {
        solver_seq.step(dt, RHO_IN, RHO_INIT);
        buffer_seq.advance(+1);
    }
    const double seq_elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - seq_start).count();

    pipe_properties_t pipe_par = pipe_properties_t::build_simple_pipe(
        { .length = PIPE_LENGTH, .dx = PIPE_DX, .diameter = PIPE_DIAM });
    std::vector<double> flow_rate_par(pipe_par.profile.get_point_count(), FLOW);
    PipeQAdvection model_par(pipe_par, flow_rate_par);
    auto buffer_par = build_linear_buffer<layer_t>(pipe_par, RHO_INIT, RHO_IN);
    quickest_parallel_solver solver_par(model_par, buffer_par.previous(), buffer_par.current());
    tbb::task_arena arena(static_cast<int>(std::thread::hardware_concurrency()));
    tbb::static_partitioner static_part;
    const auto par_start = std::chrono::steady_clock::now();
    for (size_t step = 0; step < STEP_COUNT; ++step) {
        arena.execute([&] { solver_par.step(dt, RHO_IN, RHO_INIT, static_part); });
        buffer_par.advance(+1);
    }
    print_result("arena(all)+parallel(static)", seq_elapsed, std::chrono::duration<double>(std::chrono::steady_clock::now() - par_start).count());
}

/// @brief Общие утилиты для бенчмарков 3 и 4 (расчёт 4 параметров)

/// @brief Общие сценарии расчёта четырёх параметров трубы (плотность, достоверность, присадка) разными способами параллелизации
namespace bench_4params {

using namespace quickest_parallel_bench;
using namespace pde_solvers;

/// @brief Последовательный шаг по 4 параметрам
template <typename BufferType, typename G0, typename G1, typename G2, typename G3>
void sequential_step(PipeQAdvection& model, double dt, BufferType& buffer,
    double bv0, G0 g0, double bv1, G1 g1, double bv2, G2 g2, double bv3, G3 g3)
{
    do_step<quickest_cell_compute_mode::sequential>(model, dt, bv0, buffer, g0);
    do_step<quickest_cell_compute_mode::sequential>(model, dt, bv1, buffer, g1);
    do_step<quickest_cell_compute_mode::sequential>(model, dt, bv2, buffer, g2);
    do_step<quickest_cell_compute_mode::sequential>(model, dt, bv3, buffer, g3);
}

/// @brief Параллельный шаг по 4 параметрам через parallel_invoke
template <typename BufferType, typename G0, typename G1, typename G2, typename G3>
void parallel_invoke_step(PipeQAdvection& model, double dt, BufferType& buffer,
    double bv0, G0 g0, double bv1, G1 g1, double bv2, G2 g2, double bv3, G3 g3)
{
    tbb::parallel_invoke(
        [&] { do_step<quickest_cell_compute_mode::sequential>(model, dt, bv0, buffer, g0); },
        [&] { do_step<quickest_cell_compute_mode::sequential>(model, dt, bv1, buffer, g1); },
        [&] { do_step<quickest_cell_compute_mode::sequential>(model, dt, bv2, buffer, g2); },
        [&] { do_step<quickest_cell_compute_mode::sequential>(model, dt, bv3, buffer, g3); }
    );
}

/// @brief Параллельный шаг по 4 параметрам через parallel_for с заданным партиционером
template <typename BufferType, typename G0, typename G1, typename G2, typename G3, typename Partitioner>
void parallel_for_step(PipeQAdvection& model, double dt, BufferType& buffer,
    double bv0, G0 g0, double bv1, G1 g1, double bv2, G2 g2, double bv3, G3 g3,
    Partitioner& partitioner)
{
    /// @brief Геттеры и граничные значения упакованы в массивы для индексации в лямбде
    using GType = std::common_type_t<G0, G1, G2, G3>;
    const GType getters[4] = { g0, g1, g2, g3 };
    const double bvs[4] = { bv0, bv1, bv2, bv3 };
    tbb::parallel_for(
        tbb::blocked_range<size_t>(0, 4),
        [&](const tbb::blocked_range<size_t>& r) {
            for (size_t i = r.begin(); i < r.end(); ++i)
                do_step<quickest_cell_compute_mode::sequential>(model, dt, bvs[i], buffer, getters[i]);
        },
        partitioner
    );
}

/// @brief Параллельный шаг по 4 параметрам через task_group
template <typename BufferType, typename G0, typename G1, typename G2, typename G3>
void task_group_step(PipeQAdvection& model, double dt, BufferType& buffer,
    double bv0, G0 g0, double bv1, G1 g1, double bv2, G2 g2, double bv3, G3 g3)
{
    tbb::task_group tg;
    tg.run([&] { do_step<quickest_cell_compute_mode::sequential>(model, dt, bv0, buffer, g0); });
    tg.run([&] { do_step<quickest_cell_compute_mode::sequential>(model, dt, bv1, buffer, g1); });
    tg.run([&] { do_step<quickest_cell_compute_mode::sequential>(model, dt, bv2, buffer, g2); });
    /// @brief Последнюю задачу выполняем текущим потоком через run_and_wait (scheduler bypass)
    tg.run_and_wait([&] { do_step<quickest_cell_compute_mode::sequential>(model, dt, bv3, buffer, g3); });
}

} // namespace bench_4params

/// @brief Бенчмарк 3 — 4 параметра, общий specific (как в текущей библиотеке)

/// @brief Бенчмарк с библиотечной раскладкой: все параметры и общий буфер потоков в одном слое
namespace bench3 {

using namespace quickest_parallel_bench;
using namespace bench_4params;
using namespace pde_solvers;

/// @brief Эндогенный слой как в production-коде
using endogenous_t = iso_nonbaro_improver_pipe_endogenious_layer_t;
/// @brief Слой гидравлики с эндогенными свойствами
using layer_t = hydraulic_pipe_layer<endogenous_t>;
/// @brief Кольцевой буфер prev/curr слоёв
using buffer_t = ring_buffer_t<layer_t>;

using getter_t = quickest_ultimate_fv_wrapper<1>(*)(endogenous_t&);

/// @brief Геттеры четырёх параметров трубы для передачи в общие step-функции
namespace detail {
constexpr getter_t density_getter = &endogenous_t::get_density_wrapper;
constexpr getter_t density_conf_getter = &endogenous_t::get_density_confidence_wrapper;
constexpr getter_t improver_getter = &endogenous_t::get_improver_concentration_wrapper;
constexpr getter_t improver_conf_getter = &endogenous_t::get_improver_concentration_confidence_wrapper;
} // namespace detail

/// @brief Последовательный пересчёт плотности, достоверности и присадки за один шаг по времени
inline void sequential_step(PipeQAdvection& model, double dt, buffer_t& buffer)
{
    bench_4params::sequential_step(model, dt, buffer,
        RHO_IN, detail::density_getter,
        CONF_IN, detail::density_conf_getter,
        IMPR_IN, detail::improver_getter,
        IMPR_CONF_IN, detail::improver_conf_getter);
}

/// @brief Параллельный пересчёт четырёх параметров — по одному потоку на параметр
inline void parallel_invoke_step(PipeQAdvection& model, double dt, buffer_t& buffer)
{
    bench_4params::parallel_invoke_step(model, dt, buffer,
        RHO_IN, detail::density_getter,
        CONF_IN, detail::density_conf_getter,
        IMPR_IN, detail::improver_getter,
        IMPR_CONF_IN, detail::improver_conf_getter);
}

/// @brief Параллельный пересчёт четырёх параметров с выбранной стратегией разбиения задач
template <typename Partitioner>
inline void parallel_for_step(PipeQAdvection& model, double dt, buffer_t& buffer, Partitioner& partitioner)
{
    bench_4params::parallel_for_step(model, dt, buffer,
        RHO_IN, detail::density_getter,
        CONF_IN, detail::density_conf_getter,
        IMPR_IN, detail::improver_getter,
        IMPR_CONF_IN, detail::improver_conf_getter,
        partitioner);
}

/// @brief Параллельный пересчёт четырёх параметров через пул задач TBB
inline void task_group_step(PipeQAdvection& model, double dt, buffer_t& buffer)
{
    bench_4params::task_group_step(model, dt, buffer,
        RHO_IN, detail::density_getter,
        CONF_IN, detail::density_conf_getter,
        IMPR_IN, detail::improver_getter,
        IMPR_CONF_IN, detail::improver_conf_getter);
}

} // namespace bench3

/// @brief Последовательный расчёт 4 параметров с production-раскладкой (общий specific) — baseline для bench3.
TEST(TBB_Bench3_SharedSpecific, DISABLED_Sequential_Baseline)
{
    using namespace quickest_parallel_bench;
    using namespace bench3;
    using namespace bench_4params;

    pipe_properties_t pipe_solo = pipe_properties_t::build_simple_pipe(
        { .length = PIPE_LENGTH, .dx = PIPE_DX, .diameter = PIPE_DIAM });
    std::vector<double> flow_rate_solo(pipe_solo.profile.get_point_count(), FLOW);
    PipeQAdvection model_solo(pipe_solo, flow_rate_solo);
    const double dt = calc_time_step_by_Courant(model_solo, 1.0);
    buffer_t buffer(2, pipe_solo.profile.get_point_count());

    /// @brief Типовой сценарий pde_solvers: солвер не создаётся (только step), буфер переиспользуется
    const auto solo_start = std::chrono::steady_clock::now();
    for (size_t step = 0; step < STEP_COUNT; ++step) {
        sequential_step(model_solo, dt, buffer);
        buffer.advance(+1);
    }
    print_seq("sequential-4p(shared_F)", std::chrono::duration<double>(std::chrono::steady_clock::now() - solo_start).count());
}

/// @brief 4 параметра параллельно через parallel_invoke без арены: потоки берутся из глобального пула
/// и пробуждаются заново на каждом шаге. Общий specific — возможны гонки данных на буфере F.
TEST(TBB_Bench3_SharedSpecific, DISABLED_ParallelInvoke)
{
    using namespace quickest_parallel_bench;
    using namespace bench3;
    using namespace bench_4params;

    pipe_properties_t pipe_seq = pipe_properties_t::build_simple_pipe(
        { .length = PIPE_LENGTH, .dx = PIPE_DX, .diameter = PIPE_DIAM });
    std::vector<double> flow_rate_seq(pipe_seq.profile.get_point_count(), FLOW);
    PipeQAdvection model_seq(pipe_seq, flow_rate_seq);
    const double dt = calc_time_step_by_Courant(model_seq, 1.0);
    buffer_t buffer_seq(2, pipe_seq.profile.get_point_count());
    const auto seq_start = std::chrono::steady_clock::now();
    for (size_t step = 0; step < STEP_COUNT; ++step) {
        sequential_step(model_seq, dt, buffer_seq);
        buffer_seq.advance(+1);
    }
    const double seq_elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - seq_start).count();

    /// @brief Типовой сценарий pde_solvers
    pipe_properties_t pipe_par = pipe_properties_t::build_simple_pipe(
        { .length = PIPE_LENGTH, .dx = PIPE_DX, .diameter = PIPE_DIAM });
    std::vector<double> flow_rate_par(pipe_par.profile.get_point_count(), FLOW);
    PipeQAdvection model_par(pipe_par, flow_rate_par);
    buffer_t buffer_par(2, pipe_par.profile.get_point_count());
    const auto par_start = std::chrono::steady_clock::now();
    for (size_t step = 0; step < STEP_COUNT; ++step) {
        parallel_invoke_step(model_par, dt, buffer_par);
        buffer_par.advance(+1);
    }
    print_result("parallel_invoke(shared_F) pde_solvers",
        seq_elapsed, std::chrono::duration<double>(std::chrono::steady_clock::now() - par_start).count());

    /// @brief Типовой сценарий graph_solvers: do_step создаёт solver внутри;
    /// оба сценария идентичны, вывод одного достаточен
}

/// @brief 4 параметра параллельно через parallel_for с auto_partitioner без арены.
TEST(TBB_Bench3_SharedSpecific, DISABLED_ParallelFor_Auto)
{
    using namespace quickest_parallel_bench;
    using namespace bench3;
    using namespace bench_4params;

    pipe_properties_t pipe_seq = pipe_properties_t::build_simple_pipe(
        { .length = PIPE_LENGTH, .dx = PIPE_DX, .diameter = PIPE_DIAM });
    std::vector<double> flow_rate_seq(pipe_seq.profile.get_point_count(), FLOW);
    PipeQAdvection model_seq(pipe_seq, flow_rate_seq);
    const double dt = calc_time_step_by_Courant(model_seq, 1.0);
    buffer_t buffer_seq(2, pipe_seq.profile.get_point_count());
    const auto seq_start = std::chrono::steady_clock::now();
    for (size_t step = 0; step < STEP_COUNT; ++step) {
        sequential_step(model_seq, dt, buffer_seq);
        buffer_seq.advance(+1);
    }
    const double seq_elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - seq_start).count();

    pipe_properties_t pipe_par = pipe_properties_t::build_simple_pipe(
        { .length = PIPE_LENGTH, .dx = PIPE_DX, .diameter = PIPE_DIAM });
    std::vector<double> flow_rate_par(pipe_par.profile.get_point_count(), FLOW);
    PipeQAdvection model_par(pipe_par, flow_rate_par);
    buffer_t buffer_par(2, pipe_par.profile.get_point_count());
    tbb::auto_partitioner auto_part;
    const auto par_start = std::chrono::steady_clock::now();
    for (size_t step = 0; step < STEP_COUNT; ++step) {
        parallel_for_step(model_par, dt, buffer_par, auto_part);
        buffer_par.advance(+1);
    }
    print_result("parallel_for(auto,shared_F)",
        seq_elapsed, std::chrono::duration<double>(std::chrono::steady_clock::now() - par_start).count());
}

/// @brief 4 параметра параллельно через parallel_for с static_partitioner без арены.
TEST(TBB_Bench3_SharedSpecific, DISABLED_ParallelFor_Static)
{
    using namespace quickest_parallel_bench;
    using namespace bench3;
    using namespace bench_4params;

    pipe_properties_t pipe_seq = pipe_properties_t::build_simple_pipe(
        { .length = PIPE_LENGTH, .dx = PIPE_DX, .diameter = PIPE_DIAM });
    std::vector<double> flow_rate_seq(pipe_seq.profile.get_point_count(), FLOW);
    PipeQAdvection model_seq(pipe_seq, flow_rate_seq);
    const double dt = calc_time_step_by_Courant(model_seq, 1.0);
    buffer_t buffer_seq(2, pipe_seq.profile.get_point_count());
    const auto seq_start = std::chrono::steady_clock::now();
    for (size_t step = 0; step < STEP_COUNT; ++step) {
        sequential_step(model_seq, dt, buffer_seq);
        buffer_seq.advance(+1);
    }
    const double seq_elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - seq_start).count();

    pipe_properties_t pipe_par = pipe_properties_t::build_simple_pipe(
        { .length = PIPE_LENGTH, .dx = PIPE_DX, .diameter = PIPE_DIAM });
    std::vector<double> flow_rate_par(pipe_par.profile.get_point_count(), FLOW);
    PipeQAdvection model_par(pipe_par, flow_rate_par);
    buffer_t buffer_par(2, pipe_par.profile.get_point_count());
    tbb::static_partitioner static_part;
    const auto par_start = std::chrono::steady_clock::now();
    for (size_t step = 0; step < STEP_COUNT; ++step) {
        parallel_for_step(model_par, dt, buffer_par, static_part);
        buffer_par.advance(+1);
    }
    print_result("parallel_for(static,shared_F)",
        seq_elapsed, std::chrono::duration<double>(std::chrono::steady_clock::now() - par_start).count());
}

/// @brief 4 параметра параллельно через parallel_for с affinity_partitioner без арены.
/// Два варианта: партиционер живёт вне цикла (pde_solvers) и создаётся заново на каждом шаге (graph_solvers).
TEST(TBB_Bench3_SharedSpecific, DISABLED_ParallelFor_Affinity)
{
    using namespace quickest_parallel_bench;
    using namespace bench3;
    using namespace bench_4params;

    pipe_properties_t pipe_seq = pipe_properties_t::build_simple_pipe(
        { .length = PIPE_LENGTH, .dx = PIPE_DX, .diameter = PIPE_DIAM });
    std::vector<double> flow_rate_seq(pipe_seq.profile.get_point_count(), FLOW);
    PipeQAdvection model_seq(pipe_seq, flow_rate_seq);
    const double dt = calc_time_step_by_Courant(model_seq, 1.0);
    buffer_t buffer_seq(2, pipe_seq.profile.get_point_count());
    const auto seq_start = std::chrono::steady_clock::now();
    for (size_t step = 0; step < STEP_COUNT; ++step) {
        sequential_step(model_seq, dt, buffer_seq);
        buffer_seq.advance(+1);
    }
    const double seq_elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - seq_start).count();

    pipe_properties_t pipe_par = pipe_properties_t::build_simple_pipe(
        { .length = PIPE_LENGTH, .dx = PIPE_DX, .diameter = PIPE_DIAM });
    std::vector<double> flow_rate_par(pipe_par.profile.get_point_count(), FLOW);
    PipeQAdvection model_par(pipe_par, flow_rate_par);
    buffer_t buffer_par(2, pipe_par.profile.get_point_count());

    /// @brief Типовой сценарий pde_solvers: affinity_partitioner живёт вне цикла — запоминает привязки
    tbb::affinity_partitioner affinity_part;
    const auto par_start = std::chrono::steady_clock::now();
    for (size_t step = 0; step < STEP_COUNT; ++step) {
        parallel_for_step(model_par, dt, buffer_par, affinity_part);
        buffer_par.advance(+1);
    }
    print_result("parallel_for(affinity,shared_F) pde_solvers",
        seq_elapsed, std::chrono::duration<double>(std::chrono::steady_clock::now() - par_start).count());

    /// @brief Типовой сценарий graph_solvers: affinity_partitioner создаётся заново на каждом шаге — привязки не накапливаются
    pipe_properties_t pipe_graph = pipe_properties_t::build_simple_pipe(
        { .length = PIPE_LENGTH, .dx = PIPE_DX, .diameter = PIPE_DIAM });
    std::vector<double> flow_rate_graph(pipe_graph.profile.get_point_count(), FLOW);
    PipeQAdvection model_graph(pipe_graph, flow_rate_graph);
    buffer_t buffer_graph(2, pipe_graph.profile.get_point_count());
    const auto graph_start = std::chrono::steady_clock::now();
    for (size_t step = 0; step < STEP_COUNT; ++step) {
        tbb::affinity_partitioner affinity_part_local;
        parallel_for_step(model_graph, dt, buffer_graph, affinity_part_local);
        buffer_graph.advance(+1);
    }
    print_result("parallel_for(affinity,shared_F) graph_solvers",
        seq_elapsed, std::chrono::duration<double>(std::chrono::steady_clock::now() - graph_start).count());
}

/// @brief 4 параметра параллельно через task_group: последняя задача выполняется главным потоком
/// через run_and_wait() (scheduler bypass) без арены.
TEST(TBB_Bench3_SharedSpecific, DISABLED_TaskGroup)
{
    using namespace quickest_parallel_bench;
    using namespace bench3;
    using namespace bench_4params;

    pipe_properties_t pipe_seq = pipe_properties_t::build_simple_pipe(
        { .length = PIPE_LENGTH, .dx = PIPE_DX, .diameter = PIPE_DIAM });
    std::vector<double> flow_rate_seq(pipe_seq.profile.get_point_count(), FLOW);
    PipeQAdvection model_seq(pipe_seq, flow_rate_seq);
    const double dt = calc_time_step_by_Courant(model_seq, 1.0);
    buffer_t buffer_seq(2, pipe_seq.profile.get_point_count());
    const auto seq_start = std::chrono::steady_clock::now();
    for (size_t step = 0; step < STEP_COUNT; ++step) {
        sequential_step(model_seq, dt, buffer_seq);
        buffer_seq.advance(+1);
    }
    const double seq_elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - seq_start).count();

    /// @brief Типовой сценарий pde_solvers: task_group создаётся один раз вне цикла
    pipe_properties_t pipe_par = pipe_properties_t::build_simple_pipe(
        { .length = PIPE_LENGTH, .dx = PIPE_DX, .diameter = PIPE_DIAM });
    std::vector<double> flow_rate_par(pipe_par.profile.get_point_count(), FLOW);
    PipeQAdvection model_par(pipe_par, flow_rate_par);
    buffer_t buffer_par(2, pipe_par.profile.get_point_count());
    const auto par_start = std::chrono::steady_clock::now();
    for (size_t step = 0; step < STEP_COUNT; ++step) {
        task_group_step(model_par, dt, buffer_par);
        buffer_par.advance(+1);
    }
    print_result("task_group(shared_F)",
        seq_elapsed, std::chrono::duration<double>(std::chrono::steady_clock::now() - par_start).count());
}

/// @brief 4 параметра параллельно: arena(4) вне цикла устраняет overhead пробуждения потоков,
/// static_partitioner фиксирует деление. Два варианта: arena вне цикла (pde_solvers) и внутри (graph_solvers).
TEST(TBB_Bench3_SharedSpecific, DISABLED_Arena4_ParallelFor_Static)
{
    using namespace quickest_parallel_bench;
    using namespace bench3;
    using namespace bench_4params;

    pipe_properties_t pipe_seq = pipe_properties_t::build_simple_pipe(
        { .length = PIPE_LENGTH, .dx = PIPE_DX, .diameter = PIPE_DIAM });
    std::vector<double> flow_rate_seq(pipe_seq.profile.get_point_count(), FLOW);
    PipeQAdvection model_seq(pipe_seq, flow_rate_seq);
    const double dt = calc_time_step_by_Courant(model_seq, 1.0);
    buffer_t buffer_seq(2, pipe_seq.profile.get_point_count());
    const auto seq_start = std::chrono::steady_clock::now();
    for (size_t step = 0; step < STEP_COUNT; ++step) {
        sequential_step(model_seq, dt, buffer_seq);
        buffer_seq.advance(+1);
    }
    const double seq_elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - seq_start).count();

    /// @brief Типовой сценарий pde_solvers: arena живёт вне цикла — потоки не уходят в сон
    tbb::task_arena arena(4);
    pipe_properties_t pipe_par = pipe_properties_t::build_simple_pipe(
        { .length = PIPE_LENGTH, .dx = PIPE_DX, .diameter = PIPE_DIAM });
    std::vector<double> flow_rate_par(pipe_par.profile.get_point_count(), FLOW);
    PipeQAdvection model_par(pipe_par, flow_rate_par);
    buffer_t buffer_par(2, pipe_par.profile.get_point_count());
    tbb::static_partitioner static_part;
    const auto par_start = std::chrono::steady_clock::now();
    for (size_t step = 0; step < STEP_COUNT; ++step) {
        arena.execute([&] { parallel_for_step(model_par, dt, buffer_par, static_part); });
        buffer_par.advance(+1);
    }
    print_result("arena(4)+parallel_for(static,shared_F) pde_solvers",
        seq_elapsed, std::chrono::duration<double>(std::chrono::steady_clock::now() - par_start).count());

    /// @brief Типовой сценарий graph_solvers: arena создаётся на каждом шаге
    pipe_properties_t pipe_graph = pipe_properties_t::build_simple_pipe(
        { .length = PIPE_LENGTH, .dx = PIPE_DX, .diameter = PIPE_DIAM });
    std::vector<double> flow_rate_graph(pipe_graph.profile.get_point_count(), FLOW);
    PipeQAdvection model_graph(pipe_graph, flow_rate_graph);
    buffer_t buffer_graph(2, pipe_graph.profile.get_point_count());
    const auto graph_start = std::chrono::steady_clock::now();
    for (size_t step = 0; step < STEP_COUNT; ++step) {
        tbb::task_arena arena_local(4);
        tbb::static_partitioner static_part_local;
        arena_local.execute([&] { parallel_for_step(model_graph, dt, buffer_graph, static_part_local); });
        buffer_graph.advance(+1);
    }
    print_result("arena(4)+parallel_for(static,shared_F) graph_solvers",
        seq_elapsed, std::chrono::duration<double>(std::chrono::steady_clock::now() - graph_start).count());
}

/// @brief 4 параметра параллельно: arena(4) вне цикла + parallel_invoke.
/// Два варианта: arena вне цикла (pde_solvers) и внутри (graph_solvers).
TEST(TBB_Bench3_SharedSpecific, DISABLED_Arena4_ParallelInvoke)
{
    using namespace quickest_parallel_bench;
    using namespace bench3;
    using namespace bench_4params;

    pipe_properties_t pipe_seq = pipe_properties_t::build_simple_pipe(
        { .length = PIPE_LENGTH, .dx = PIPE_DX, .diameter = PIPE_DIAM });
    std::vector<double> flow_rate_seq(pipe_seq.profile.get_point_count(), FLOW);
    PipeQAdvection model_seq(pipe_seq, flow_rate_seq);
    const double dt = calc_time_step_by_Courant(model_seq, 1.0);
    buffer_t buffer_seq(2, pipe_seq.profile.get_point_count());
    const auto seq_start = std::chrono::steady_clock::now();
    for (size_t step = 0; step < STEP_COUNT; ++step) {
        sequential_step(model_seq, dt, buffer_seq);
        buffer_seq.advance(+1);
    }
    const double seq_elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - seq_start).count();

    /// @brief Типовой сценарий pde_solvers: arena живёт вне цикла
    tbb::task_arena arena(4);
    pipe_properties_t pipe_par = pipe_properties_t::build_simple_pipe(
        { .length = PIPE_LENGTH, .dx = PIPE_DX, .diameter = PIPE_DIAM });
    std::vector<double> flow_rate_par(pipe_par.profile.get_point_count(), FLOW);
    PipeQAdvection model_par(pipe_par, flow_rate_par);
    buffer_t buffer_par(2, pipe_par.profile.get_point_count());
    const auto par_start = std::chrono::steady_clock::now();
    for (size_t step = 0; step < STEP_COUNT; ++step) {
        arena.execute([&] { parallel_invoke_step(model_par, dt, buffer_par); });
        buffer_par.advance(+1);
    }
    print_result("arena(4)+parallel_invoke(shared_F) pde_solvers",
        seq_elapsed, std::chrono::duration<double>(std::chrono::steady_clock::now() - par_start).count());

    /// @brief Типовой сценарий graph_solvers: arena создаётся на каждом шаге
    pipe_properties_t pipe_graph = pipe_properties_t::build_simple_pipe(
        { .length = PIPE_LENGTH, .dx = PIPE_DX, .diameter = PIPE_DIAM });
    std::vector<double> flow_rate_graph(pipe_graph.profile.get_point_count(), FLOW);
    PipeQAdvection model_graph(pipe_graph, flow_rate_graph);
    buffer_t buffer_graph(2, pipe_graph.profile.get_point_count());
    const auto graph_start = std::chrono::steady_clock::now();
    for (size_t step = 0; step < STEP_COUNT; ++step) {
        tbb::task_arena arena_local(4);
        arena_local.execute([&] { parallel_invoke_step(model_graph, dt, buffer_graph); });
        buffer_graph.advance(+1);
    }
    print_result("arena(4)+parallel_invoke(shared_F) graph_solvers",
        seq_elapsed, std::chrono::duration<double>(std::chrono::steady_clock::now() - graph_start).count());
}

/// @brief 4 параметра параллельно: arena(2) — проверка гипотезы, что 2 потока дадут лучший результат
/// при общем specific, т.к. два параметра вместе помещаются в L1.
/// Два варианта: arena вне цикла (pde_solvers) и внутри (graph_solvers).
TEST(TBB_Bench3_SharedSpecific, DISABLED_Arena2_ParallelFor_Static)
{
    using namespace quickest_parallel_bench;
    using namespace bench3;
    using namespace bench_4params;

    pipe_properties_t pipe_seq = pipe_properties_t::build_simple_pipe(
        { .length = PIPE_LENGTH, .dx = PIPE_DX, .diameter = PIPE_DIAM });
    std::vector<double> flow_rate_seq(pipe_seq.profile.get_point_count(), FLOW);
    PipeQAdvection model_seq(pipe_seq, flow_rate_seq);
    const double dt = calc_time_step_by_Courant(model_seq, 1.0);
    buffer_t buffer_seq(2, pipe_seq.profile.get_point_count());
    const auto seq_start = std::chrono::steady_clock::now();
    for (size_t step = 0; step < STEP_COUNT; ++step) {
        sequential_step(model_seq, dt, buffer_seq);
        buffer_seq.advance(+1);
    }
    const double seq_elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - seq_start).count();

    /// @brief Типовой сценарий pde_solvers: arena живёт вне цикла
    tbb::task_arena arena(2);
    pipe_properties_t pipe_par = pipe_properties_t::build_simple_pipe(
        { .length = PIPE_LENGTH, .dx = PIPE_DX, .diameter = PIPE_DIAM });
    std::vector<double> flow_rate_par(pipe_par.profile.get_point_count(), FLOW);
    PipeQAdvection model_par(pipe_par, flow_rate_par);
    buffer_t buffer_par(2, pipe_par.profile.get_point_count());
    tbb::static_partitioner static_part;
    const auto par_start = std::chrono::steady_clock::now();
    for (size_t step = 0; step < STEP_COUNT; ++step) {
        arena.execute([&] { parallel_for_step(model_par, dt, buffer_par, static_part); });
        buffer_par.advance(+1);
    }
    print_result("arena(2)+parallel_for(static,shared_F) pde_solvers",
        seq_elapsed, std::chrono::duration<double>(std::chrono::steady_clock::now() - par_start).count());

    /// @brief Типовой сценарий graph_solvers: arena создаётся на каждом шаге
    pipe_properties_t pipe_graph = pipe_properties_t::build_simple_pipe(
        { .length = PIPE_LENGTH, .dx = PIPE_DX, .diameter = PIPE_DIAM });
    std::vector<double> flow_rate_graph(pipe_graph.profile.get_point_count(), FLOW);
    PipeQAdvection model_graph(pipe_graph, flow_rate_graph);
    buffer_t buffer_graph(2, pipe_graph.profile.get_point_count());
    const auto graph_start = std::chrono::steady_clock::now();
    for (size_t step = 0; step < STEP_COUNT; ++step) {
        tbb::task_arena arena_local(2);
        tbb::static_partitioner static_part_local;
        arena_local.execute([&] { parallel_for_step(model_graph, dt, buffer_graph, static_part_local); });
        buffer_graph.advance(+1);
    }
    print_result("arena(2)+parallel_for(static,shared_F) graph_solvers",
        seq_elapsed, std::chrono::duration<double>(std::chrono::steady_clock::now() - graph_start).count());
}

/// @brief Бенчмарк 4 — 4 параметра, собственный specific на каждый параметр

/// @brief Бенчмарк с раздельными буферами потоков: у каждого параметра свой specific, без гонок при параллелизме
namespace bench4 {

using namespace quickest_parallel_bench;
using namespace bench_4params;
using namespace pde_solvers;

/// @brief Эндогенный слой с изолированными буферами потоков
using endogenous_t = separate_specific_layer::endogenous_layer_t;
using layer_t = hydraulic_pipe_layer<endogenous_t>;
using buffer_t = ring_buffer_t<layer_t>;

/// @brief Геттеры четырёх параметров трубы с раздельными буферами потоков
namespace detail {
constexpr auto density_getter = &endogenous_t::get_density_wrapper;
constexpr auto density_conf_getter = &endogenous_t::get_density_confidence_wrapper;
constexpr auto improver_getter = &endogenous_t::get_improver_wrapper;
constexpr auto improver_conf_getter = &endogenous_t::get_improver_confidence_wrapper;
} // namespace detail

/// @brief Последовательный пересчёт плотности, достоверности и присадки за один шаг по времени
inline void sequential_step(PipeQAdvection& model, double dt, buffer_t& buffer)
{
    bench_4params::sequential_step(model, dt, buffer,
        RHO_IN, detail::density_getter,
        CONF_IN, detail::density_conf_getter,
        IMPR_IN, detail::improver_getter,
        IMPR_CONF_IN, detail::improver_conf_getter);
}

/// @brief Параллельный пересчёт четырёх параметров — по одному потоку на параметр
inline void parallel_invoke_step(PipeQAdvection& model, double dt, buffer_t& buffer)
{
    bench_4params::parallel_invoke_step(model, dt, buffer,
        RHO_IN, detail::density_getter,
        CONF_IN, detail::density_conf_getter,
        IMPR_IN, detail::improver_getter,
        IMPR_CONF_IN, detail::improver_conf_getter);
}

/// @brief Параллельный пересчёт четырёх параметров с выбранной стратегией разбиения задач
template <typename Partitioner>
inline void parallel_for_step(PipeQAdvection& model, double dt, buffer_t& buffer, Partitioner& partitioner)
{
    bench_4params::parallel_for_step(model, dt, buffer,
        RHO_IN, detail::density_getter,
        CONF_IN, detail::density_conf_getter,
        IMPR_IN, detail::improver_getter,
        IMPR_CONF_IN, detail::improver_conf_getter,
        partitioner);
}

/// @brief Параллельный пересчёт четырёх параметров через пул задач TBB
inline void task_group_step(PipeQAdvection& model, double dt, buffer_t& buffer)
{
    bench_4params::task_group_step(model, dt, buffer,
        RHO_IN, detail::density_getter,
        CONF_IN, detail::density_conf_getter,
        IMPR_IN, detail::improver_getter,
        IMPR_CONF_IN, detail::improver_conf_getter);
}

} // namespace bench4

/// @brief Последовательный расчёт 4 параметров с раздельным specific на каждый — baseline для bench4.
TEST(TBB_Bench4_SplitSpecific, DISABLED_Sequential_Baseline)
{
    using namespace quickest_parallel_bench;
    using namespace bench4;
    using namespace bench_4params;

    pipe_properties_t pipe_solo = pipe_properties_t::build_simple_pipe(
        { .length = PIPE_LENGTH, .dx = PIPE_DX, .diameter = PIPE_DIAM });
    std::vector<double> flow_rate_solo(pipe_solo.profile.get_point_count(), FLOW);
    PipeQAdvection model_solo(pipe_solo, flow_rate_solo);
    const double dt = calc_time_step_by_Courant(model_solo, 1.0);
    buffer_t buffer(2, pipe_solo.profile.get_point_count());

    /// @brief Типовой сценарий pde_solvers
    const auto solo_start = std::chrono::steady_clock::now();
    for (size_t step = 0; step < STEP_COUNT; ++step) {
        sequential_step(model_solo, dt, buffer);
        buffer.advance(+1);
    }
    print_seq("sequential-4p(split_F)", std::chrono::duration<double>(std::chrono::steady_clock::now() - solo_start).count());
}

/// @brief 4 параметра параллельно через parallel_invoke без арены. Раздельный specific исключает
/// гонки данных на буфере F, в отличие от bench3.
TEST(TBB_Bench4_SplitSpecific, DISABLED_ParallelInvoke)
{
    using namespace quickest_parallel_bench;
    using namespace bench4;
    using namespace bench_4params;

    pipe_properties_t pipe_seq = pipe_properties_t::build_simple_pipe(
        { .length = PIPE_LENGTH, .dx = PIPE_DX, .diameter = PIPE_DIAM });
    std::vector<double> flow_rate_seq(pipe_seq.profile.get_point_count(), FLOW);
    PipeQAdvection model_seq(pipe_seq, flow_rate_seq);
    const double dt = calc_time_step_by_Courant(model_seq, 1.0);
    buffer_t buffer_seq(2, pipe_seq.profile.get_point_count());
    const auto seq_start = std::chrono::steady_clock::now();
    for (size_t step = 0; step < STEP_COUNT; ++step) {
        sequential_step(model_seq, dt, buffer_seq);
        buffer_seq.advance(+1);
    }
    const double seq_elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - seq_start).count();

    /// @brief Типовой сценарий pde_solvers
    pipe_properties_t pipe_par = pipe_properties_t::build_simple_pipe(
        { .length = PIPE_LENGTH, .dx = PIPE_DX, .diameter = PIPE_DIAM });
    std::vector<double> flow_rate_par(pipe_par.profile.get_point_count(), FLOW);
    PipeQAdvection model_par(pipe_par, flow_rate_par);
    buffer_t buffer_par(2, pipe_par.profile.get_point_count());
    const auto par_start = std::chrono::steady_clock::now();
    for (size_t step = 0; step < STEP_COUNT; ++step) {
        parallel_invoke_step(model_par, dt, buffer_par);
        buffer_par.advance(+1);
    }
    print_result("parallel_invoke(split_F) pde_solvers",
        seq_elapsed, std::chrono::duration<double>(std::chrono::steady_clock::now() - par_start).count());
}

/// @brief 4 параметра параллельно через parallel_for с auto_partitioner без арены, раздельный specific.
TEST(TBB_Bench4_SplitSpecific, DISABLED_ParallelFor_Auto)
{
    using namespace quickest_parallel_bench;
    using namespace bench4;
    using namespace bench_4params;

    pipe_properties_t pipe_seq = pipe_properties_t::build_simple_pipe(
        { .length = PIPE_LENGTH, .dx = PIPE_DX, .diameter = PIPE_DIAM });
    std::vector<double> flow_rate_seq(pipe_seq.profile.get_point_count(), FLOW);
    PipeQAdvection model_seq(pipe_seq, flow_rate_seq);
    const double dt = calc_time_step_by_Courant(model_seq, 1.0);
    buffer_t buffer_seq(2, pipe_seq.profile.get_point_count());
    const auto seq_start = std::chrono::steady_clock::now();
    for (size_t step = 0; step < STEP_COUNT; ++step) {
        sequential_step(model_seq, dt, buffer_seq);
        buffer_seq.advance(+1);
    }
    const double seq_elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - seq_start).count();

    pipe_properties_t pipe_par = pipe_properties_t::build_simple_pipe(
        { .length = PIPE_LENGTH, .dx = PIPE_DX, .diameter = PIPE_DIAM });
    std::vector<double> flow_rate_par(pipe_par.profile.get_point_count(), FLOW);
    PipeQAdvection model_par(pipe_par, flow_rate_par);
    buffer_t buffer_par(2, pipe_par.profile.get_point_count());
    tbb::auto_partitioner auto_part;
    const auto par_start = std::chrono::steady_clock::now();
    for (size_t step = 0; step < STEP_COUNT; ++step) {
        parallel_for_step(model_par, dt, buffer_par, auto_part);
        buffer_par.advance(+1);
    }
    print_result("parallel_for(auto,split_F)",
        seq_elapsed, std::chrono::duration<double>(std::chrono::steady_clock::now() - par_start).count());
}

/// @brief 4 параметра параллельно через parallel_for с static_partitioner без арены, раздельный specific.
TEST(TBB_Bench4_SplitSpecific, DISABLED_ParallelFor_Static)
{
    using namespace quickest_parallel_bench;
    using namespace bench4;
    using namespace bench_4params;

    pipe_properties_t pipe_seq = pipe_properties_t::build_simple_pipe(
        { .length = PIPE_LENGTH, .dx = PIPE_DX, .diameter = PIPE_DIAM });
    std::vector<double> flow_rate_seq(pipe_seq.profile.get_point_count(), FLOW);
    PipeQAdvection model_seq(pipe_seq, flow_rate_seq);
    const double dt = calc_time_step_by_Courant(model_seq, 1.0);
    buffer_t buffer_seq(2, pipe_seq.profile.get_point_count());
    const auto seq_start = std::chrono::steady_clock::now();
    for (size_t step = 0; step < STEP_COUNT; ++step) {
        sequential_step(model_seq, dt, buffer_seq);
        buffer_seq.advance(+1);
    }
    const double seq_elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - seq_start).count();

    pipe_properties_t pipe_par = pipe_properties_t::build_simple_pipe(
        { .length = PIPE_LENGTH, .dx = PIPE_DX, .diameter = PIPE_DIAM });
    std::vector<double> flow_rate_par(pipe_par.profile.get_point_count(), FLOW);
    PipeQAdvection model_par(pipe_par, flow_rate_par);
    buffer_t buffer_par(2, pipe_par.profile.get_point_count());
    tbb::static_partitioner static_part;
    const auto par_start = std::chrono::steady_clock::now();
    for (size_t step = 0; step < STEP_COUNT; ++step) {
        parallel_for_step(model_par, dt, buffer_par, static_part);
        buffer_par.advance(+1);
    }
    print_result("parallel_for(static,split_F)",
        seq_elapsed, std::chrono::duration<double>(std::chrono::steady_clock::now() - par_start).count());
}

/// @brief 4 параметра параллельно через parallel_for с affinity_partitioner без арены, раздельный specific.
/// Два варианта: партиционер вне цикла (pde_solvers) и внутри (graph_solvers).
TEST(TBB_Bench4_SplitSpecific, DISABLED_ParallelFor_Affinity)
{
    using namespace quickest_parallel_bench;
    using namespace bench4;
    using namespace bench_4params;

    pipe_properties_t pipe_seq = pipe_properties_t::build_simple_pipe(
        { .length = PIPE_LENGTH, .dx = PIPE_DX, .diameter = PIPE_DIAM });
    std::vector<double> flow_rate_seq(pipe_seq.profile.get_point_count(), FLOW);
    PipeQAdvection model_seq(pipe_seq, flow_rate_seq);
    const double dt = calc_time_step_by_Courant(model_seq, 1.0);
    buffer_t buffer_seq(2, pipe_seq.profile.get_point_count());
    const auto seq_start = std::chrono::steady_clock::now();
    for (size_t step = 0; step < STEP_COUNT; ++step) {
        sequential_step(model_seq, dt, buffer_seq);
        buffer_seq.advance(+1);
    }
    const double seq_elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - seq_start).count();

    /// @brief Типовой сценарий pde_solvers: affinity_partitioner живёт вне цикла
    pipe_properties_t pipe_par = pipe_properties_t::build_simple_pipe(
        { .length = PIPE_LENGTH, .dx = PIPE_DX, .diameter = PIPE_DIAM });
    std::vector<double> flow_rate_par(pipe_par.profile.get_point_count(), FLOW);
    PipeQAdvection model_par(pipe_par, flow_rate_par);
    buffer_t buffer_par(2, pipe_par.profile.get_point_count());
    tbb::affinity_partitioner affinity_part;
    const auto par_start = std::chrono::steady_clock::now();
    for (size_t step = 0; step < STEP_COUNT; ++step) {
        parallel_for_step(model_par, dt, buffer_par, affinity_part);
        buffer_par.advance(+1);
    }
    print_result("parallel_for(affinity,split_F) pde_solvers",
        seq_elapsed, std::chrono::duration<double>(std::chrono::steady_clock::now() - par_start).count());

    /// @brief Типовой сценарий graph_solvers: affinity_partitioner создаётся заново на каждом шаге
    pipe_properties_t pipe_graph = pipe_properties_t::build_simple_pipe(
        { .length = PIPE_LENGTH, .dx = PIPE_DX, .diameter = PIPE_DIAM });
    std::vector<double> flow_rate_graph(pipe_graph.profile.get_point_count(), FLOW);
    PipeQAdvection model_graph(pipe_graph, flow_rate_graph);
    buffer_t buffer_graph(2, pipe_graph.profile.get_point_count());
    const auto graph_start = std::chrono::steady_clock::now();
    for (size_t step = 0; step < STEP_COUNT; ++step) {
        tbb::affinity_partitioner affinity_part_local;
        parallel_for_step(model_graph, dt, buffer_graph, affinity_part_local);
        buffer_graph.advance(+1);
    }
    print_result("parallel_for(affinity,split_F) graph_solvers",
        seq_elapsed, std::chrono::duration<double>(std::chrono::steady_clock::now() - graph_start).count());
}

/// @brief 4 параметра параллельно через task_group + run_and_wait без арены, раздельный specific.
TEST(TBB_Bench4_SplitSpecific, DISABLED_TaskGroup)
{
    using namespace quickest_parallel_bench;
    using namespace bench4;
    using namespace bench_4params;

    pipe_properties_t pipe_seq = pipe_properties_t::build_simple_pipe(
        { .length = PIPE_LENGTH, .dx = PIPE_DX, .diameter = PIPE_DIAM });
    std::vector<double> flow_rate_seq(pipe_seq.profile.get_point_count(), FLOW);
    PipeQAdvection model_seq(pipe_seq, flow_rate_seq);
    const double dt = calc_time_step_by_Courant(model_seq, 1.0);
    buffer_t buffer_seq(2, pipe_seq.profile.get_point_count());
    const auto seq_start = std::chrono::steady_clock::now();
    for (size_t step = 0; step < STEP_COUNT; ++step) {
        sequential_step(model_seq, dt, buffer_seq);
        buffer_seq.advance(+1);
    }
    const double seq_elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - seq_start).count();

    pipe_properties_t pipe_par = pipe_properties_t::build_simple_pipe(
        { .length = PIPE_LENGTH, .dx = PIPE_DX, .diameter = PIPE_DIAM });
    std::vector<double> flow_rate_par(pipe_par.profile.get_point_count(), FLOW);
    PipeQAdvection model_par(pipe_par, flow_rate_par);
    buffer_t buffer_par(2, pipe_par.profile.get_point_count());
    const auto par_start = std::chrono::steady_clock::now();
    for (size_t step = 0; step < STEP_COUNT; ++step) {
        task_group_step(model_par, dt, buffer_par);
        buffer_par.advance(+1);
    }
    print_result("task_group(split_F)",
        seq_elapsed, std::chrono::duration<double>(std::chrono::steady_clock::now() - par_start).count());
}

/// @brief 4 параметра параллельно: arena(4) + static_partitioner + раздельный specific.
/// Два варианта: arena вне цикла (pde_solvers) и внутри (graph_solvers).
TEST(TBB_Bench4_SplitSpecific, DISABLED_Arena4_ParallelFor_Static)
{
    using namespace quickest_parallel_bench;
    using namespace bench4;
    using namespace bench_4params;

    pipe_properties_t pipe_seq = pipe_properties_t::build_simple_pipe(
        { .length = PIPE_LENGTH, .dx = PIPE_DX, .diameter = PIPE_DIAM });
    std::vector<double> flow_rate_seq(pipe_seq.profile.get_point_count(), FLOW);
    PipeQAdvection model_seq(pipe_seq, flow_rate_seq);
    const double dt = calc_time_step_by_Courant(model_seq, 1.0);
    buffer_t buffer_seq(2, pipe_seq.profile.get_point_count());
    const auto seq_start = std::chrono::steady_clock::now();
    for (size_t step = 0; step < STEP_COUNT; ++step) {
        sequential_step(model_seq, dt, buffer_seq);
        buffer_seq.advance(+1);
    }
    const double seq_elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - seq_start).count();

    /// @brief Типовой сценарий pde_solvers: arena живёт вне цикла
    tbb::task_arena arena(4);
    pipe_properties_t pipe_par = pipe_properties_t::build_simple_pipe(
        { .length = PIPE_LENGTH, .dx = PIPE_DX, .diameter = PIPE_DIAM });
    std::vector<double> flow_rate_par(pipe_par.profile.get_point_count(), FLOW);
    PipeQAdvection model_par(pipe_par, flow_rate_par);
    buffer_t buffer_par(2, pipe_par.profile.get_point_count());
    tbb::static_partitioner static_part;
    const auto par_start = std::chrono::steady_clock::now();
    for (size_t step = 0; step < STEP_COUNT; ++step) {
        arena.execute([&] { parallel_for_step(model_par, dt, buffer_par, static_part); });
        buffer_par.advance(+1);
    }
    print_result("arena(4)+parallel_for(static,split_F) pde_solvers",
        seq_elapsed, std::chrono::duration<double>(std::chrono::steady_clock::now() - par_start).count());

    /// @brief Типовой сценарий graph_solvers: arena создаётся на каждом шаге
    pipe_properties_t pipe_graph = pipe_properties_t::build_simple_pipe(
        { .length = PIPE_LENGTH, .dx = PIPE_DX, .diameter = PIPE_DIAM });
    std::vector<double> flow_rate_graph(pipe_graph.profile.get_point_count(), FLOW);
    PipeQAdvection model_graph(pipe_graph, flow_rate_graph);
    buffer_t buffer_graph(2, pipe_graph.profile.get_point_count());
    const auto graph_start = std::chrono::steady_clock::now();
    for (size_t step = 0; step < STEP_COUNT; ++step) {
        tbb::task_arena arena_local(4);
        tbb::static_partitioner static_part_local;
        arena_local.execute([&] { parallel_for_step(model_graph, dt, buffer_graph, static_part_local); });
        buffer_graph.advance(+1);
    }
    print_result("arena(4)+parallel_for(static,split_F) graph_solvers",
        seq_elapsed, std::chrono::duration<double>(std::chrono::steady_clock::now() - graph_start).count());
}

/// @brief 4 параметра параллельно: arena(4) + affinity_partitioner + раздельный specific — лучший
/// результат среди всех бенчмарков (6715 ms, 2.12×).
/// Два варианта: arena и партиционер вне цикла (pde_solvers) и внутри (graph_solvers).
TEST(TBB_Bench4_SplitSpecific, DISABLED_Arena4_ParallelFor_Affinity)
{
    using namespace quickest_parallel_bench;
    using namespace bench4;
    using namespace bench_4params;

    pipe_properties_t pipe_seq = pipe_properties_t::build_simple_pipe(
        { .length = PIPE_LENGTH, .dx = PIPE_DX, .diameter = PIPE_DIAM });
    std::vector<double> flow_rate_seq(pipe_seq.profile.get_point_count(), FLOW);
    PipeQAdvection model_seq(pipe_seq, flow_rate_seq);
    const double dt = calc_time_step_by_Courant(model_seq, 1.0);
    buffer_t buffer_seq(2, pipe_seq.profile.get_point_count());
    const auto seq_start = std::chrono::steady_clock::now();
    for (size_t step = 0; step < STEP_COUNT; ++step) {
        sequential_step(model_seq, dt, buffer_seq);
        buffer_seq.advance(+1);
    }
    const double seq_elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - seq_start).count();

    /// @brief Типовой сценарий pde_solvers: arena + affinity_partitioner живут вне цикла
    tbb::task_arena arena(4);
    pipe_properties_t pipe_par = pipe_properties_t::build_simple_pipe(
        { .length = PIPE_LENGTH, .dx = PIPE_DX, .diameter = PIPE_DIAM });
    std::vector<double> flow_rate_par(pipe_par.profile.get_point_count(), FLOW);
    PipeQAdvection model_par(pipe_par, flow_rate_par);
    buffer_t buffer_par(2, pipe_par.profile.get_point_count());
    tbb::affinity_partitioner affinity_part;
    const auto par_start = std::chrono::steady_clock::now();
    for (size_t step = 0; step < STEP_COUNT; ++step) {
        arena.execute([&] { parallel_for_step(model_par, dt, buffer_par, affinity_part); });
        buffer_par.advance(+1);
    }
    print_result("arena(4)+parallel_for(affinity,split_F) pde_solvers",
        seq_elapsed, std::chrono::duration<double>(std::chrono::steady_clock::now() - par_start).count());

    /// @brief Типовой сценарий graph_solvers: arena и affinity_partitioner создаются на каждом шаге
    pipe_properties_t pipe_graph = pipe_properties_t::build_simple_pipe(
        { .length = PIPE_LENGTH, .dx = PIPE_DX, .diameter = PIPE_DIAM });
    std::vector<double> flow_rate_graph(pipe_graph.profile.get_point_count(), FLOW);
    PipeQAdvection model_graph(pipe_graph, flow_rate_graph);
    buffer_t buffer_graph(2, pipe_graph.profile.get_point_count());
    const auto graph_start = std::chrono::steady_clock::now();
    for (size_t step = 0; step < STEP_COUNT; ++step) {
        tbb::task_arena arena_local(4);
        tbb::affinity_partitioner affinity_part_local;
        arena_local.execute([&] { parallel_for_step(model_graph, dt, buffer_graph, affinity_part_local); });
        buffer_graph.advance(+1);
    }
    print_result("arena(4)+parallel_for(affinity,split_F) graph_solvers",
        seq_elapsed, std::chrono::duration<double>(std::chrono::steady_clock::now() - graph_start).count());
}

/// @brief 4 параметра параллельно: arena(4) + parallel_invoke + раздельный specific.
/// Два варианта: arena вне цикла (pde_solvers) и внутри (graph_solvers).
TEST(TBB_Bench4_SplitSpecific, DISABLED_Arena4_ParallelInvoke)
{
    using namespace quickest_parallel_bench;
    using namespace bench4;
    using namespace bench_4params;

    pipe_properties_t pipe_seq = pipe_properties_t::build_simple_pipe(
        { .length = PIPE_LENGTH, .dx = PIPE_DX, .diameter = PIPE_DIAM });
    std::vector<double> flow_rate_seq(pipe_seq.profile.get_point_count(), FLOW);
    PipeQAdvection model_seq(pipe_seq, flow_rate_seq);
    const double dt = calc_time_step_by_Courant(model_seq, 1.0);
    buffer_t buffer_seq(2, pipe_seq.profile.get_point_count());
    const auto seq_start = std::chrono::steady_clock::now();
    for (size_t step = 0; step < STEP_COUNT; ++step) {
        sequential_step(model_seq, dt, buffer_seq);
        buffer_seq.advance(+1);
    }
    const double seq_elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - seq_start).count();

    /// @brief Типовой сценарий pde_solvers: arena живёт вне цикла
    tbb::task_arena arena(4);
    pipe_properties_t pipe_par = pipe_properties_t::build_simple_pipe(
        { .length = PIPE_LENGTH, .dx = PIPE_DX, .diameter = PIPE_DIAM });
    std::vector<double> flow_rate_par(pipe_par.profile.get_point_count(), FLOW);
    PipeQAdvection model_par(pipe_par, flow_rate_par);
    buffer_t buffer_par(2, pipe_par.profile.get_point_count());
    const auto par_start = std::chrono::steady_clock::now();
    for (size_t step = 0; step < STEP_COUNT; ++step) {
        arena.execute([&] { parallel_invoke_step(model_par, dt, buffer_par); });
        buffer_par.advance(+1);
    }
    print_result("arena(4)+parallel_invoke(split_F) pde_solvers",
        seq_elapsed, std::chrono::duration<double>(std::chrono::steady_clock::now() - par_start).count());

    /// @brief Типовой сценарий graph_solvers: arena создаётся на каждом шаге
    pipe_properties_t pipe_graph = pipe_properties_t::build_simple_pipe(
        { .length = PIPE_LENGTH, .dx = PIPE_DX, .diameter = PIPE_DIAM });
    std::vector<double> flow_rate_graph(pipe_graph.profile.get_point_count(), FLOW);
    PipeQAdvection model_graph(pipe_graph, flow_rate_graph);
    buffer_t buffer_graph(2, pipe_graph.profile.get_point_count());
    const auto graph_start = std::chrono::steady_clock::now();
    for (size_t step = 0; step < STEP_COUNT; ++step) {
        tbb::task_arena arena_local(4);
        arena_local.execute([&] { parallel_invoke_step(model_graph, dt, buffer_graph); });
        buffer_graph.advance(+1);
    }
    print_result("arena(4)+parallel_invoke(split_F) graph_solvers",
        seq_elapsed, std::chrono::duration<double>(std::chrono::steady_clock::now() - graph_start).count());
}

/// @brief 4 параметра параллельно: arena(2) + static_partitioner + раздельный specific.
/// Два варианта: arena вне цикла (pde_solvers) и внутри (graph_solvers).
TEST(TBB_Bench4_SplitSpecific, DISABLED_Arena2_ParallelFor_Static)
{
    using namespace quickest_parallel_bench;
    using namespace bench4;
    using namespace bench_4params;

    pipe_properties_t pipe_seq = pipe_properties_t::build_simple_pipe(
        { .length = PIPE_LENGTH, .dx = PIPE_DX, .diameter = PIPE_DIAM });
    std::vector<double> flow_rate_seq(pipe_seq.profile.get_point_count(), FLOW);
    PipeQAdvection model_seq(pipe_seq, flow_rate_seq);
    const double dt = calc_time_step_by_Courant(model_seq, 1.0);
    buffer_t buffer_seq(2, pipe_seq.profile.get_point_count());
    const auto seq_start = std::chrono::steady_clock::now();
    for (size_t step = 0; step < STEP_COUNT; ++step) {
        sequential_step(model_seq, dt, buffer_seq);
        buffer_seq.advance(+1);
    }
    const double seq_elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - seq_start).count();

    /// @brief Типовой сценарий pde_solvers: arena живёт вне цикла
    tbb::task_arena arena(2);
    pipe_properties_t pipe_par = pipe_properties_t::build_simple_pipe(
        { .length = PIPE_LENGTH, .dx = PIPE_DX, .diameter = PIPE_DIAM });
    std::vector<double> flow_rate_par(pipe_par.profile.get_point_count(), FLOW);
    PipeQAdvection model_par(pipe_par, flow_rate_par);
    buffer_t buffer_par(2, pipe_par.profile.get_point_count());
    tbb::static_partitioner static_part;
    const auto par_start = std::chrono::steady_clock::now();
    for (size_t step = 0; step < STEP_COUNT; ++step) {
        arena.execute([&] { parallel_for_step(model_par, dt, buffer_par, static_part); });
        buffer_par.advance(+1);
    }
    print_result("arena(2)+parallel_for(static,split_F) pde_solvers",
        seq_elapsed, std::chrono::duration<double>(std::chrono::steady_clock::now() - par_start).count());

    /// @brief Типовой сценарий graph_solvers: arena создаётся на каждом шаге
    pipe_properties_t pipe_graph = pipe_properties_t::build_simple_pipe(
        { .length = PIPE_LENGTH, .dx = PIPE_DX, .diameter = PIPE_DIAM });
    std::vector<double> flow_rate_graph(pipe_graph.profile.get_point_count(), FLOW);
    PipeQAdvection model_graph(pipe_graph, flow_rate_graph);
    buffer_t buffer_graph(2, pipe_graph.profile.get_point_count());
    const auto graph_start = std::chrono::steady_clock::now();
    for (size_t step = 0; step < STEP_COUNT; ++step) {
        tbb::task_arena arena_local(2);
        tbb::static_partitioner static_part_local;
        arena_local.execute([&] { parallel_for_step(model_graph, dt, buffer_graph, static_part_local); });
        buffer_graph.advance(+1);
    }
    print_result("arena(2)+parallel_for(static,split_F) graph_solvers",
        seq_elapsed, std::chrono::duration<double>(std::chrono::steady_clock::now() - graph_start).count());
}

/// @brief Бенчмарк 5 — кэш-дружелюбная раскладка: U_prev и U_curr одного параметра в одном блоке памяти;
/// каждый параметр изолирован в single_param_data_t с alignas(64), отдельный specific на параметр

/// @brief Бенчмарк с разнесением параметров по памяти: каждый поток работает со своим непрерывным блоком данных
namespace bench5 {

using namespace quickest_parallel_bench;
using namespace bench_4params;
using namespace pde_solvers;

/// @brief Буфер потоков на границах ячеек
using specific_layer = quickest_ultimate_fv_solver_traits<1>::specific_layer;

/// @brief Один адвектируемый параметр трубы: предыдущий и текущий профиль изолированы от остальных параметров
struct alignas(64) single_param_data_t {
    std::vector<double> U_prev;
    std::vector<double> U_curr;
    specific_layer      spec;

    explicit single_param_data_t(size_t point_count)
        : U_prev(point_count - 1, 850.0)
        , U_curr(point_count - 1, 850.0)
        , spec(point_count)
    {}

    single_param_data_t(const single_param_data_t&) = delete;
    single_param_data_t& operator=(const single_param_data_t&) = delete;
    single_param_data_t(single_param_data_t&&) = default;

    /// @brief Сдвиг по времени: предыдущий профиль становится текущим без копирования
    void advance() { std::swap(U_prev, U_curr); }

    /// @brief Профиль параметра на новом временном слое для солвера
    quickest_ultimate_fv_wrapper<1> get_curr_wrapper()
        { return quickest_ultimate_fv_wrapper<1>(U_curr, spec); }

    /// @brief Профиль параметра на предыдущем временном слое для солвера
    quickest_ultimate_fv_wrapper<1> get_prev_wrapper()
        { return quickest_ultimate_fv_wrapper<1>(U_prev, spec); }
};

/// @brief Полный набор адвектируемых свойств нефти: плотность, достоверность, присадка и её достоверность
struct quad_params_t {
    single_param_data_t density;
    single_param_data_t density_conf;
    single_param_data_t improver;
    single_param_data_t improver_conf;

    explicit quad_params_t(size_t point_count)
        : density(point_count)
        , density_conf(point_count)
        , improver(point_count)
        , improver_conf(point_count)
    {}

    /// @brief Синхронный сдвиг всех четырёх параметров на следующий шаг по времени
    void advance() {
        density.advance();
        density_conf.advance();
        improver.advance();
        improver_conf.advance();
    }
};

/// @brief Один шаг адвекции выбранного параметра (плотность, присадка и т.д.) с заданным граничным значением
inline void do_param_step(PipeQAdvection& model, double dt, double bv, single_param_data_t& p)
{
    auto prev_w = p.get_prev_wrapper();
    auto curr_w = p.get_curr_wrapper();
    /// @brief Ring_buffer_t из двух обёрток prev/curr для вызова библиотечного солвера
    std::vector<quickest_ultimate_fv_wrapper<1>> wrappers{ prev_w, curr_w };
    ring_buffer_t<quickest_ultimate_fv_wrapper<1>> wrap_buf(wrappers, 1);
    quickest_ultimate_fv_solver<quickest_cell_compute_mode::sequential> solver(model, wrap_buf);
    solver.step(dt, bv, bv);
}

/// @brief Последовательный пересчёт всех четырёх параметров с заданными граничными условиями
inline void sequential_step(PipeQAdvection& model, double dt, quad_params_t& qp,
    double bv0, double bv1, double bv2, double bv3)
{
    do_param_step(model, dt, bv0, qp.density);
    do_param_step(model, dt, bv1, qp.density_conf);
    do_param_step(model, dt, bv2, qp.improver);
    do_param_step(model, dt, bv3, qp.improver_conf);
}

/// @brief Параллельный пересчёт четырёх параметров — по одному потоку на параметр
inline void parallel_invoke_step(PipeQAdvection& model, double dt, quad_params_t& qp,
    double bv0, double bv1, double bv2, double bv3)
{
    tbb::parallel_invoke(
        [&] { do_param_step(model, dt, bv0, qp.density); },
        [&] { do_param_step(model, dt, bv1, qp.density_conf); },
        [&] { do_param_step(model, dt, bv2, qp.improver); },
        [&] { do_param_step(model, dt, bv3, qp.improver_conf); }
    );
}

/// @brief Параллельный пересчёт четырёх параметров с выбранной стратегией разбиения задач
template <typename Partitioner>
void parallel_for_step(PipeQAdvection& model, double dt, quad_params_t& qp,
    double bv0, double bv1, double bv2, double bv3, Partitioner& partitioner)
{
    single_param_data_t* params[4] = { &qp.density, &qp.density_conf, &qp.improver, &qp.improver_conf };
    const double         bvs[4]    = { bv0, bv1, bv2, bv3 };
    tbb::parallel_for(
        tbb::blocked_range<size_t>(0, 4),
        [&](const tbb::blocked_range<size_t>& r) {
            for (size_t i = r.begin(); i < r.end(); ++i)
                do_param_step(model, dt, bvs[i], *params[i]);
        },
        partitioner
    );
}

/// @brief Параллельный пересчёт четырёх параметров через пул задач TBB
inline void task_group_step(PipeQAdvection& model, double dt, quad_params_t& qp,
    double bv0, double bv1, double bv2, double bv3)
{
    tbb::task_group tg;
    tg.run([&] { do_param_step(model, dt, bv0, qp.density); });
    tg.run([&] { do_param_step(model, dt, bv1, qp.density_conf); });
    tg.run([&] { do_param_step(model, dt, bv2, qp.improver); });
    tg.run_and_wait([&] { do_param_step(model, dt, bv3, qp.improver_conf); });
}

/// @brief Последовательный шаг с граничными условиями по умолчанию
inline void sequential_step(PipeQAdvection& model, double dt, quad_params_t& qp)
{
    sequential_step(model, dt, qp, RHO_IN, CONF_IN, IMPR_IN, IMPR_CONF_IN);
}

/// @brief Параллельный шаг с граничными условиями по умолчанию
inline void parallel_invoke_step(PipeQAdvection& model, double dt, quad_params_t& qp)
{
    parallel_invoke_step(model, dt, qp, RHO_IN, CONF_IN, IMPR_IN, IMPR_CONF_IN);
}

/// @brief Параллельный шаг с граничными условиями по умолчанию и выбранным партиционером
template <typename Partitioner>
inline void parallel_for_step(PipeQAdvection& model, double dt, quad_params_t& qp, Partitioner& partitioner)
{
    parallel_for_step(model, dt, qp, RHO_IN, CONF_IN, IMPR_IN, IMPR_CONF_IN, partitioner);
}

/// @brief Параллельный шаг через task_group с граничными условиями по умолчанию
inline void task_group_step(PipeQAdvection& model, double dt, quad_params_t& qp)
{
    task_group_step(model, dt, qp, RHO_IN, CONF_IN, IMPR_IN, IMPR_CONF_IN);
}

} // namespace bench5

/// @brief Тесты Bench5: кэш-дружелюбная раскладка (отдельный буфер на параметр);
/// те же сценарии, что Bench4, но данные сгруппированы по параметрам

/// @brief Последовательный расчёт 4 параметров с параметро-ориентированной раскладкой (alignas(64)) — baseline для bench5.
TEST(TBB_Bench5_CacheLayout, DISABLED_Sequential_Baseline)
{
    using namespace quickest_parallel_bench;
    using namespace bench5;

    pipe_properties_t pipe_solo = pipe_properties_t::build_simple_pipe(
        { .length = PIPE_LENGTH, .dx = PIPE_DX, .diameter = PIPE_DIAM });
    std::vector<double> flow_rate_solo(pipe_solo.profile.get_point_count(), FLOW);
    PipeQAdvection model_solo(pipe_solo, flow_rate_solo);
    const double dt = calc_time_step_by_Courant(model_solo, 1.0);
    quad_params_t qb(pipe_solo.profile.get_point_count());

    /// @brief Типовой сценарий pde_solvers
    const auto solo_start = std::chrono::steady_clock::now();
    for (size_t step = 0; step < STEP_COUNT; ++step) {
        sequential_step(model_solo, dt, qb);
        qb.advance();
    }
    print_seq("sequential-4p(cache_layout)", std::chrono::duration<double>(std::chrono::steady_clock::now() - solo_start).count());
}

/// @brief 4 параметра параллельно через parallel_invoke без арены, alignas(64)-раскладка.
TEST(TBB_Bench5_CacheLayout, DISABLED_ParallelInvoke)
{
    using namespace quickest_parallel_bench;
    using namespace bench5;

    pipe_properties_t pipe_seq = pipe_properties_t::build_simple_pipe(
        { .length = PIPE_LENGTH, .dx = PIPE_DX, .diameter = PIPE_DIAM });
    std::vector<double> flow_rate_seq(pipe_seq.profile.get_point_count(), FLOW);
    PipeQAdvection model_seq(pipe_seq, flow_rate_seq);
    const double dt = calc_time_step_by_Courant(model_seq, 1.0);
    quad_params_t quad_seq(pipe_seq.profile.get_point_count());
    const auto seq_start = std::chrono::steady_clock::now();
    for (size_t step = 0; step < STEP_COUNT; ++step) {
        sequential_step(model_seq, dt, quad_seq);
        quad_seq.advance();
    }
    const double seq_elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - seq_start).count();

    /// @brief Типовой сценарий pde_solvers
    pipe_properties_t pipe_par = pipe_properties_t::build_simple_pipe(
        { .length = PIPE_LENGTH, .dx = PIPE_DX, .diameter = PIPE_DIAM });
    std::vector<double> flow_rate_par(pipe_par.profile.get_point_count(), FLOW);
    PipeQAdvection model_par(pipe_par, flow_rate_par);
    quad_params_t quad_par(pipe_par.profile.get_point_count());
    const auto par_start = std::chrono::steady_clock::now();
    for (size_t step = 0; step < STEP_COUNT; ++step) {
        parallel_invoke_step(model_par, dt, quad_par);
        quad_par.advance();
    }
    print_result("parallel_invoke(cache_layout) pde_solvers",
        seq_elapsed, std::chrono::duration<double>(std::chrono::steady_clock::now() - par_start).count());
}

/// @brief 4 параметра параллельно через parallel_for с auto_partitioner без арены, alignas(64)-раскладка.
TEST(TBB_Bench5_CacheLayout, DISABLED_ParallelFor_Auto)
{
    using namespace quickest_parallel_bench;
    using namespace bench5;

    pipe_properties_t pipe_seq = pipe_properties_t::build_simple_pipe(
        { .length = PIPE_LENGTH, .dx = PIPE_DX, .diameter = PIPE_DIAM });
    std::vector<double> flow_rate_seq(pipe_seq.profile.get_point_count(), FLOW);
    PipeQAdvection model_seq(pipe_seq, flow_rate_seq);
    const double dt = calc_time_step_by_Courant(model_seq, 1.0);
    quad_params_t quad_seq(pipe_seq.profile.get_point_count());
    const auto seq_start = std::chrono::steady_clock::now();
    for (size_t step = 0; step < STEP_COUNT; ++step) {
        sequential_step(model_seq, dt, quad_seq);
        quad_seq.advance();
    }
    const double seq_elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - seq_start).count();

    pipe_properties_t pipe_par = pipe_properties_t::build_simple_pipe(
        { .length = PIPE_LENGTH, .dx = PIPE_DX, .diameter = PIPE_DIAM });
    std::vector<double> flow_rate_par(pipe_par.profile.get_point_count(), FLOW);
    PipeQAdvection model_par(pipe_par, flow_rate_par);
    quad_params_t quad_par(pipe_par.profile.get_point_count());
    tbb::auto_partitioner auto_part;
    const auto par_start = std::chrono::steady_clock::now();
    for (size_t step = 0; step < STEP_COUNT; ++step) {
        parallel_for_step(model_par, dt, quad_par, auto_part);
        quad_par.advance();
    }
    print_result("parallel_for(auto,cache_layout)",
        seq_elapsed, std::chrono::duration<double>(std::chrono::steady_clock::now() - par_start).count());
}

/// @brief 4 параметра параллельно через parallel_for с static_partitioner без арены, alignas(64)-раскладка.
TEST(TBB_Bench5_CacheLayout, DISABLED_ParallelFor_Static)
{
    using namespace quickest_parallel_bench;
    using namespace bench5;

    pipe_properties_t pipe_seq = pipe_properties_t::build_simple_pipe(
        { .length = PIPE_LENGTH, .dx = PIPE_DX, .diameter = PIPE_DIAM });
    std::vector<double> flow_rate_seq(pipe_seq.profile.get_point_count(), FLOW);
    PipeQAdvection model_seq(pipe_seq, flow_rate_seq);
    const double dt = calc_time_step_by_Courant(model_seq, 1.0);
    quad_params_t quad_seq(pipe_seq.profile.get_point_count());
    const auto seq_start = std::chrono::steady_clock::now();
    for (size_t step = 0; step < STEP_COUNT; ++step) {
        sequential_step(model_seq, dt, quad_seq);
        quad_seq.advance();
    }
    const double seq_elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - seq_start).count();

    pipe_properties_t pipe_par = pipe_properties_t::build_simple_pipe(
        { .length = PIPE_LENGTH, .dx = PIPE_DX, .diameter = PIPE_DIAM });
    std::vector<double> flow_rate_par(pipe_par.profile.get_point_count(), FLOW);
    PipeQAdvection model_par(pipe_par, flow_rate_par);
    quad_params_t quad_par(pipe_par.profile.get_point_count());
    tbb::static_partitioner static_part;
    const auto par_start = std::chrono::steady_clock::now();
    for (size_t step = 0; step < STEP_COUNT; ++step) {
        parallel_for_step(model_par, dt, quad_par, static_part);
        quad_par.advance();
    }
    print_result("parallel_for(static,cache_layout)",
        seq_elapsed, std::chrono::duration<double>(std::chrono::steady_clock::now() - par_start).count());
}

/// @brief 4 параметра параллельно через parallel_for с affinity_partitioner без арены, alignas(64)-раскладка.
/// Два варианта: партиционер вне цикла (pde_solvers) и внутри (graph_solvers).
TEST(TBB_Bench5_CacheLayout, DISABLED_ParallelFor_Affinity)
{
    using namespace quickest_parallel_bench;
    using namespace bench5;

    pipe_properties_t pipe_seq = pipe_properties_t::build_simple_pipe(
        { .length = PIPE_LENGTH, .dx = PIPE_DX, .diameter = PIPE_DIAM });
    std::vector<double> flow_rate_seq(pipe_seq.profile.get_point_count(), FLOW);
    PipeQAdvection model_seq(pipe_seq, flow_rate_seq);
    const double dt = calc_time_step_by_Courant(model_seq, 1.0);
    quad_params_t quad_seq(pipe_seq.profile.get_point_count());
    const auto seq_start = std::chrono::steady_clock::now();
    for (size_t step = 0; step < STEP_COUNT; ++step) {
        sequential_step(model_seq, dt, quad_seq);
        quad_seq.advance();
    }
    const double seq_elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - seq_start).count();

    /// @brief Типовой сценарий pde_solvers: affinity_partitioner живёт вне цикла
    pipe_properties_t pipe_par = pipe_properties_t::build_simple_pipe(
        { .length = PIPE_LENGTH, .dx = PIPE_DX, .diameter = PIPE_DIAM });
    std::vector<double> flow_rate_par(pipe_par.profile.get_point_count(), FLOW);
    PipeQAdvection model_par(pipe_par, flow_rate_par);
    quad_params_t quad_par(pipe_par.profile.get_point_count());
    tbb::affinity_partitioner affinity_part;
    const auto par_start = std::chrono::steady_clock::now();
    for (size_t step = 0; step < STEP_COUNT; ++step) {
        parallel_for_step(model_par, dt, quad_par, affinity_part);
        quad_par.advance();
    }
    print_result("parallel_for(affinity,cache_layout) pde_solvers",
        seq_elapsed, std::chrono::duration<double>(std::chrono::steady_clock::now() - par_start).count());

    /// @brief Типовой сценарий graph_solvers: affinity_partitioner создаётся заново на каждом шаге
    pipe_properties_t pipe_graph = pipe_properties_t::build_simple_pipe(
        { .length = PIPE_LENGTH, .dx = PIPE_DX, .diameter = PIPE_DIAM });
    std::vector<double> flow_rate_graph(pipe_graph.profile.get_point_count(), FLOW);
    PipeQAdvection model_graph(pipe_graph, flow_rate_graph);
    quad_params_t quad_graph(pipe_graph.profile.get_point_count());
    const auto graph_start = std::chrono::steady_clock::now();
    for (size_t step = 0; step < STEP_COUNT; ++step) {
        tbb::affinity_partitioner affinity_part_local;
        parallel_for_step(model_graph, dt, quad_graph, affinity_part_local);
        quad_graph.advance();
    }
    print_result("parallel_for(affinity,cache_layout) graph_solvers",
        seq_elapsed, std::chrono::duration<double>(std::chrono::steady_clock::now() - graph_start).count());
}

/// @brief 4 параметра параллельно через task_group + run_and_wait без арены, alignas(64)-раскладка.
TEST(TBB_Bench5_CacheLayout, DISABLED_TaskGroup)
{
    using namespace quickest_parallel_bench;
    using namespace bench5;

    pipe_properties_t pipe_seq = pipe_properties_t::build_simple_pipe(
        { .length = PIPE_LENGTH, .dx = PIPE_DX, .diameter = PIPE_DIAM });
    std::vector<double> flow_rate_seq(pipe_seq.profile.get_point_count(), FLOW);
    PipeQAdvection model_seq(pipe_seq, flow_rate_seq);
    const double dt = calc_time_step_by_Courant(model_seq, 1.0);
    quad_params_t quad_seq(pipe_seq.profile.get_point_count());
    const auto seq_start = std::chrono::steady_clock::now();
    for (size_t step = 0; step < STEP_COUNT; ++step) {
        sequential_step(model_seq, dt, quad_seq);
        quad_seq.advance();
    }
    const double seq_elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - seq_start).count();

    pipe_properties_t pipe_par = pipe_properties_t::build_simple_pipe(
        { .length = PIPE_LENGTH, .dx = PIPE_DX, .diameter = PIPE_DIAM });
    std::vector<double> flow_rate_par(pipe_par.profile.get_point_count(), FLOW);
    PipeQAdvection model_par(pipe_par, flow_rate_par);
    quad_params_t quad_par(pipe_par.profile.get_point_count());
    const auto par_start = std::chrono::steady_clock::now();
    for (size_t step = 0; step < STEP_COUNT; ++step) {
        task_group_step(model_par, dt, quad_par);
        quad_par.advance();
    }
    print_result("task_group(cache_layout)",
        seq_elapsed, std::chrono::duration<double>(std::chrono::steady_clock::now() - par_start).count());
}

/// @brief 4 параметра параллельно: arena(4) + static_partitioner + alignas(64)-раскладка.
/// Два варианта: arena вне цикла (pde_solvers) и внутри (graph_solvers).
TEST(TBB_Bench5_CacheLayout, DISABLED_Arena4_ParallelFor_Static)
{
    using namespace quickest_parallel_bench;
    using namespace bench5;

    pipe_properties_t pipe_seq = pipe_properties_t::build_simple_pipe(
        { .length = PIPE_LENGTH, .dx = PIPE_DX, .diameter = PIPE_DIAM });
    std::vector<double> flow_rate_seq(pipe_seq.profile.get_point_count(), FLOW);
    PipeQAdvection model_seq(pipe_seq, flow_rate_seq);
    const double dt = calc_time_step_by_Courant(model_seq, 1.0);
    quad_params_t quad_seq(pipe_seq.profile.get_point_count());
    const auto seq_start = std::chrono::steady_clock::now();
    for (size_t step = 0; step < STEP_COUNT; ++step) {
        sequential_step(model_seq, dt, quad_seq);
        quad_seq.advance();
    }
    const double seq_elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - seq_start).count();

    /// @brief Типовой сценарий pde_solvers: arena живёт вне цикла
    tbb::task_arena arena(4);
    pipe_properties_t pipe_par = pipe_properties_t::build_simple_pipe(
        { .length = PIPE_LENGTH, .dx = PIPE_DX, .diameter = PIPE_DIAM });
    std::vector<double> flow_rate_par(pipe_par.profile.get_point_count(), FLOW);
    PipeQAdvection model_par(pipe_par, flow_rate_par);
    quad_params_t quad_par(pipe_par.profile.get_point_count());
    tbb::static_partitioner static_part;
    const auto par_start = std::chrono::steady_clock::now();
    for (size_t step = 0; step < STEP_COUNT; ++step) {
        arena.execute([&] { parallel_for_step(model_par, dt, quad_par, static_part); });
        quad_par.advance();
    }
    print_result("arena(4)+parallel_for(static,cache_layout) pde_solvers",
        seq_elapsed, std::chrono::duration<double>(std::chrono::steady_clock::now() - par_start).count());

    /// @brief Типовой сценарий graph_solvers: arena создаётся на каждом шаге
    pipe_properties_t pipe_graph = pipe_properties_t::build_simple_pipe(
        { .length = PIPE_LENGTH, .dx = PIPE_DX, .diameter = PIPE_DIAM });
    std::vector<double> flow_rate_graph(pipe_graph.profile.get_point_count(), FLOW);
    PipeQAdvection model_graph(pipe_graph, flow_rate_graph);
    quad_params_t quad_graph(pipe_graph.profile.get_point_count());
    const auto graph_start = std::chrono::steady_clock::now();
    for (size_t step = 0; step < STEP_COUNT; ++step) {
        tbb::task_arena arena_local(4);
        tbb::static_partitioner static_part_local;
        arena_local.execute([&] { parallel_for_step(model_graph, dt, quad_graph, static_part_local); });
        quad_graph.advance();
    }
    print_result("arena(4)+parallel_for(static,cache_layout) graph_solvers",
        seq_elapsed, std::chrono::duration<double>(std::chrono::steady_clock::now() - graph_start).count());
}

/// @brief 4 параметра параллельно: arena(4) + affinity_partitioner + alignas(64)-раскладка —
/// лучший результат bench5 (7047 ms, 2.21×). Аффинность устойчива благодаря изоляции кэш-линий.
/// Два варианта: arena и партиционер вне цикла (pde_solvers) и внутри (graph_solvers).
TEST(TBB_Bench5_CacheLayout, DISABLED_Arena4_ParallelFor_Affinity)
{
    using namespace quickest_parallel_bench;
    using namespace bench5;

    pipe_properties_t pipe_seq = pipe_properties_t::build_simple_pipe(
        { .length = PIPE_LENGTH, .dx = PIPE_DX, .diameter = PIPE_DIAM });
    std::vector<double> flow_rate_seq(pipe_seq.profile.get_point_count(), FLOW);
    PipeQAdvection model_seq(pipe_seq, flow_rate_seq);
    const double dt = calc_time_step_by_Courant(model_seq, 1.0);
    quad_params_t quad_seq(pipe_seq.profile.get_point_count());
    const auto seq_start = std::chrono::steady_clock::now();
    for (size_t step = 0; step < STEP_COUNT; ++step) {
        sequential_step(model_seq, dt, quad_seq);
        quad_seq.advance();
    }
    const double seq_elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - seq_start).count();

    /// @brief Типовой сценарий pde_solvers: arena + affinity_partitioner живут вне цикла
    tbb::task_arena arena(4);
    pipe_properties_t pipe_par = pipe_properties_t::build_simple_pipe(
        { .length = PIPE_LENGTH, .dx = PIPE_DX, .diameter = PIPE_DIAM });
    std::vector<double> flow_rate_par(pipe_par.profile.get_point_count(), FLOW);
    PipeQAdvection model_par(pipe_par, flow_rate_par);
    quad_params_t quad_par(pipe_par.profile.get_point_count());
    tbb::affinity_partitioner affinity_part;
    const auto par_start = std::chrono::steady_clock::now();
    for (size_t step = 0; step < STEP_COUNT; ++step) {
        arena.execute([&] { parallel_for_step(model_par, dt, quad_par, affinity_part); });
        quad_par.advance();
    }
    print_result("arena(4)+parallel_for(affinity,cache_layout) pde_solvers",
        seq_elapsed, std::chrono::duration<double>(std::chrono::steady_clock::now() - par_start).count());

    /// @brief Типовой сценарий graph_solvers: arena и affinity_partitioner создаются на каждом шаге
    pipe_properties_t pipe_graph = pipe_properties_t::build_simple_pipe(
        { .length = PIPE_LENGTH, .dx = PIPE_DX, .diameter = PIPE_DIAM });
    std::vector<double> flow_rate_graph(pipe_graph.profile.get_point_count(), FLOW);
    PipeQAdvection model_graph(pipe_graph, flow_rate_graph);
    quad_params_t quad_graph(pipe_graph.profile.get_point_count());
    const auto graph_start = std::chrono::steady_clock::now();
    for (size_t step = 0; step < STEP_COUNT; ++step) {
        tbb::task_arena arena_local(4);
        tbb::affinity_partitioner affinity_part_local;
        arena_local.execute([&] { parallel_for_step(model_graph, dt, quad_graph, affinity_part_local); });
        quad_graph.advance();
    }
    print_result("arena(4)+parallel_for(affinity,cache_layout) graph_solvers",
        seq_elapsed, std::chrono::duration<double>(std::chrono::steady_clock::now() - graph_start).count());
}

/// @brief 4 параметра параллельно: arena(4) + parallel_invoke + alignas(64)-раскладка.
/// Два варианта: arena вне цикла (pde_solvers) и внутри (graph_solvers).
TEST(TBB_Bench5_CacheLayout, DISABLED_Arena4_ParallelInvoke)
{
    using namespace quickest_parallel_bench;
    using namespace bench5;

    pipe_properties_t pipe_seq = pipe_properties_t::build_simple_pipe(
        { .length = PIPE_LENGTH, .dx = PIPE_DX, .diameter = PIPE_DIAM });
    std::vector<double> flow_rate_seq(pipe_seq.profile.get_point_count(), FLOW);
    PipeQAdvection model_seq(pipe_seq, flow_rate_seq);
    const double dt = calc_time_step_by_Courant(model_seq, 1.0);
    quad_params_t quad_seq(pipe_seq.profile.get_point_count());
    const auto seq_start = std::chrono::steady_clock::now();
    for (size_t step = 0; step < STEP_COUNT; ++step) {
        sequential_step(model_seq, dt, quad_seq);
        quad_seq.advance();
    }
    const double seq_elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - seq_start).count();

    /// @brief Типовой сценарий pde_solvers: arena живёт вне цикла
    tbb::task_arena arena(4);
    pipe_properties_t pipe_par = pipe_properties_t::build_simple_pipe(
        { .length = PIPE_LENGTH, .dx = PIPE_DX, .diameter = PIPE_DIAM });
    std::vector<double> flow_rate_par(pipe_par.profile.get_point_count(), FLOW);
    PipeQAdvection model_par(pipe_par, flow_rate_par);
    quad_params_t quad_par(pipe_par.profile.get_point_count());
    const auto par_start = std::chrono::steady_clock::now();
    for (size_t step = 0; step < STEP_COUNT; ++step) {
        arena.execute([&] { parallel_invoke_step(model_par, dt, quad_par); });
        quad_par.advance();
    }
    print_result("arena(4)+parallel_invoke(cache_layout) pde_solvers",
        seq_elapsed, std::chrono::duration<double>(std::chrono::steady_clock::now() - par_start).count());

    /// @brief Типовой сценарий graph_solvers: arena создаётся на каждом шаге
    pipe_properties_t pipe_graph = pipe_properties_t::build_simple_pipe(
        { .length = PIPE_LENGTH, .dx = PIPE_DX, .diameter = PIPE_DIAM });
    std::vector<double> flow_rate_graph(pipe_graph.profile.get_point_count(), FLOW);
    PipeQAdvection model_graph(pipe_graph, flow_rate_graph);
    quad_params_t quad_graph(pipe_graph.profile.get_point_count());
    const auto graph_start = std::chrono::steady_clock::now();
    for (size_t step = 0; step < STEP_COUNT; ++step) {
        tbb::task_arena arena_local(4);
        arena_local.execute([&] { parallel_invoke_step(model_graph, dt, quad_graph); });
        quad_graph.advance();
    }
    print_result("arena(4)+parallel_invoke(cache_layout) graph_solvers",
        seq_elapsed, std::chrono::duration<double>(std::chrono::steady_clock::now() - graph_start).count());
}

/// @brief 4 параметра параллельно: arena(2) + static_partitioner + alignas(64)-раскладка.
/// Два варианта: arena вне цикла (pde_solvers) и внутри (graph_solvers).
TEST(TBB_Bench5_CacheLayout, DISABLED_Arena2_ParallelFor_Static)
{
    using namespace quickest_parallel_bench;
    using namespace bench5;

    pipe_properties_t pipe_seq = pipe_properties_t::build_simple_pipe(
        { .length = PIPE_LENGTH, .dx = PIPE_DX, .diameter = PIPE_DIAM });
    std::vector<double> flow_rate_seq(pipe_seq.profile.get_point_count(), FLOW);
    PipeQAdvection model_seq(pipe_seq, flow_rate_seq);
    const double dt = calc_time_step_by_Courant(model_seq, 1.0);
    quad_params_t quad_seq(pipe_seq.profile.get_point_count());
    const auto seq_start = std::chrono::steady_clock::now();
    for (size_t step = 0; step < STEP_COUNT; ++step) {
        sequential_step(model_seq, dt, quad_seq);
        quad_seq.advance();
    }
    const double seq_elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - seq_start).count();

    /// @brief Типовой сценарий pde_solvers: arena живёт вне цикла
    tbb::task_arena arena(2);
    pipe_properties_t pipe_par = pipe_properties_t::build_simple_pipe(
        { .length = PIPE_LENGTH, .dx = PIPE_DX, .diameter = PIPE_DIAM });
    std::vector<double> flow_rate_par(pipe_par.profile.get_point_count(), FLOW);
    PipeQAdvection model_par(pipe_par, flow_rate_par);
    quad_params_t quad_par(pipe_par.profile.get_point_count());
    tbb::static_partitioner static_part;
    const auto par_start = std::chrono::steady_clock::now();
    for (size_t step = 0; step < STEP_COUNT; ++step) {
        arena.execute([&] { parallel_for_step(model_par, dt, quad_par, static_part); });
        quad_par.advance();
    }
    print_result("arena(2)+parallel_for(static,cache_layout) pde_solvers",
        seq_elapsed, std::chrono::duration<double>(std::chrono::steady_clock::now() - par_start).count());

    /// @brief Типовой сценарий graph_solvers: arena создаётся на каждом шаге
    pipe_properties_t pipe_graph = pipe_properties_t::build_simple_pipe(
        { .length = PIPE_LENGTH, .dx = PIPE_DX, .diameter = PIPE_DIAM });
    std::vector<double> flow_rate_graph(pipe_graph.profile.get_point_count(), FLOW);
    PipeQAdvection model_graph(pipe_graph, flow_rate_graph);
    quad_params_t quad_graph(pipe_graph.profile.get_point_count());
    const auto graph_start = std::chrono::steady_clock::now();
    for (size_t step = 0; step < STEP_COUNT; ++step) {
        tbb::task_arena arena_local(2);
        tbb::static_partitioner static_part_local;
        arena_local.execute([&] { parallel_for_step(model_graph, dt, quad_graph, static_part_local); });
        quad_graph.advance();
    }
    print_result("arena(2)+parallel_for(static,cache_layout) graph_solvers",
        seq_elapsed, std::chrono::duration<double>(std::chrono::steady_clock::now() - graph_start).count());
}

/// @brief Альтернативный wrapper и solver для непрерывных буферов в памяти
/// (quickest_ultimate_span_wrapper/solver вместо std::vector<double>&)

namespace pde_solvers {

/// @brief Wrapper над одним параметром QUICKEST-ULTIMATE через raw-указатели.
/// Совместим с любым непрерывным буфером — в отличие от quickest_ultimate_fv_wrapper<1>,
/// который требует std::vector<double>&.
struct quickest_ultimate_span_wrapper {
    /// @brief Профиль параметра по ячейкам трубы (плотность, присадка и т.д.)
    double* U;
    /// @brief Потоки параметра на границах ячеек
    double* F;
    /// @brief Число расчётных ячеек
    size_t  cell_count;
    /// @brief Число граничных точек (на одну больше числа ячеек)
    size_t  point_count;
};

/// @brief Солвер QUICKEST-ULTIMATE, работающий с quickest_ultimate_span_wrapper;
/// логика step() идентична quickest_ultimate_fv_solver, данные — через raw-указатели
/// @tparam QuickestMode Режим последовательной или параллельной обработки ячеек
template <quickest_cell_compute_mode QuickestMode>
class quickest_ultimate_span_solver {
    pde_t<1>&                          pde;
    const std::vector<double>&         grid;
    const quickest_ultimate_span_wrapper prev;
    quickest_ultimate_span_wrapper       curr;
public:
    /// @brief Солвер привязывается к модели трубы и двум временным срезам одного параметра
    quickest_ultimate_span_solver(
        pde_t<1>& pde,
        const quickest_ultimate_span_wrapper& prev,
        const quickest_ultimate_span_wrapper& curr)
        : pde(pde)
        , grid(pde.get_grid())
        , prev(prev)
        , curr(curr)
    {}

    /// @brief Один шаг адвекции параметра с граничными значениями на входе и выходе трубы
    void step(double dt, double u_in, double u_out)
    {
        const double* U     = prev.U;
        double*       U_new = curr.U;
        double*       F     = curr.F;
        const size_t  cells = prev.cell_count;
        const size_t  points = prev.point_count;

        double v_in  = pde.getEquationsCoeffs(0, U[0]);
        double v_out = pde.getEquationsCoeffs(points - 1, U[cells - 1]);
        if (v_in  >= 0) F[0]        = v_in  * u_in;
        if (v_out <= 0) F[points-1] = v_out * u_out;

        double v_pipe = pde.getEquationsCoeffs(0, U[0]);
        if (v_pipe >= 0) {
            for_each_cell<QuickestMode>(cells, [&](size_t cell) {
                double Ub = (cell == 0 || cell == cells - 1)
                    ? U[cell]
                    : quickest_ultimate_border_approximation(
                        U[cell-1], U[cell], U[cell+1],
                        0, grid[cell+1] - grid[cell], dt, v_pipe);
                F[cell + 1] = Ub * v_pipe;
            });
        }
        else {
            for_each_cell<QuickestMode>(cells, [&](size_t cell) {
                double Ub = (cell == 0 || cell == cells - 1)
                    ? U[cell]
                    : quickest_ultimate_border_approximation(
                        U[cell+1], U[cell], U[cell-1],
                        0, grid[cell+1] - grid[cell], dt, std::abs(v_pipe));
                F[cell] = Ub * v_pipe;
            });
        }

        for_each_cell<QuickestMode>(cells, [&](size_t cell) {
            double dx = grid[cell+1] - grid[cell];
            double v_cell = pde.getEquationsCoeffs(cell, U[cell]);
            double Cr = std::abs(v_cell * dt / dx);
            if (Cr > 1.0 + std::numeric_limits<double>::epsilon())
                throw std::runtime_error("Quickest-ultimate is called with Cr > 1");
            U_new[cell] = U[cell] + dt / dx * (F[cell] - F[cell+1])
                        + dt * pde.getSourceTerm(cell, U[cell]);
        });
    }
};

} // namespace pde_solvers

/// @brief Бенчмарк 6 — плоский буфер [U_prev | U_curr | F] в одном аллоке на параметр

/// @brief Бенчмарк с плоским буфером: профиль, новый слой и потоки одного параметра в одном блоке памяти
namespace bench6 {

using namespace quickest_parallel_bench;
using namespace pde_solvers;

/// @brief Данные одного адвектируемого параметра: предыдущий и текущий профиль + потоки в одном векторе
struct alignas(64) single_param_data_t {
    /// @brief Непрерывный буфер [U_prev | U_curr | F]
    std::vector<double> storage;
    /// @brief Профиль на предыдущем шаге
    double* U_prev;
    /// @brief Профиль на текущем шаге
    double* U_curr;
    /// @brief Потоки на границах ячеек
    double* F;
    size_t  cell_count;
    size_t  point_count;

    /// @brief Инициализация профиля: равномерная плотность 850 кг/м³ по всей трубе
    explicit single_param_data_t(size_t point_count_)
        : storage(2 * (point_count_ - 1) + point_count_, 0.0)
        , U_prev(storage.data())
        , U_curr(storage.data() + (point_count_ - 1))
        , F(storage.data() + 2 * (point_count_ - 1))
        , cell_count(point_count_ - 1)
        , point_count(point_count_)
    {
        std::fill(U_prev, U_prev + cell_count, 850.0);
        std::fill(U_curr, U_curr + cell_count, 850.0);
    }
    single_param_data_t(const single_param_data_t&) = delete;
    single_param_data_t& operator=(const single_param_data_t&) = delete;
    single_param_data_t(single_param_data_t&&) = default;

    /// @brief Сдвиг по времени без перекладки данных в памяти
    void advance() { std::swap(U_prev, U_curr); }

    /// @brief Представление предыдущего временного слоя для солвера
    quickest_ultimate_span_wrapper prev_wrapper() const
        { return { U_prev, F, cell_count, point_count }; }
    /// @brief Представление текущего временного слоя для солвера
    quickest_ultimate_span_wrapper curr_wrapper() const
        { return { U_curr, F, cell_count, point_count }; }
};

/// @brief Полный набор адвектируемых свойств нефти в плоской раскладке памяти
struct quad_params_t {
    single_param_data_t density;
    single_param_data_t density_conf;
    single_param_data_t improver;
    single_param_data_t improver_conf;

    explicit quad_params_t(size_t point_count)
        : density(point_count)
        , density_conf(point_count)
        , improver(point_count)
        , improver_conf(point_count)
    {}

    /// @brief Синхронный сдвиг всех четырёх параметров на следующий шаг по времени
    void advance() {
        density.advance();
        density_conf.advance();
        improver.advance();
        improver_conf.advance();
    }
};

/// @brief Один шаг адвекции параметра через span-солвер без промежуточных обёрток
inline void do_param_step(PipeQAdvection& model, double dt, double bv, single_param_data_t& p)
{
    quickest_ultimate_span_solver<quickest_cell_compute_mode::sequential> solver(
        model, p.prev_wrapper(), p.curr_wrapper());
    solver.step(dt, bv, bv);
}

/// @brief Последовательный пересчёт всех четырёх параметров за шаг
inline void sequential_step(PipeQAdvection& model, double dt, quad_params_t& qp,
    double bv0, double bv1, double bv2, double bv3)
{
    do_param_step(model, dt, bv0, qp.density);
    do_param_step(model, dt, bv1, qp.density_conf);
    do_param_step(model, dt, bv2, qp.improver);
    do_param_step(model, dt, bv3, qp.improver_conf);
}

/// @brief Параллельный пересчёт четырёх параметров — по одному потоку на параметр
inline void parallel_invoke_step(PipeQAdvection& model, double dt, quad_params_t& qp,
    double bv0, double bv1, double bv2, double bv3)
{
    tbb::parallel_invoke(
        [&] { do_param_step(model, dt, bv0, qp.density); },
        [&] { do_param_step(model, dt, bv1, qp.density_conf); },
        [&] { do_param_step(model, dt, bv2, qp.improver); },
        [&] { do_param_step(model, dt, bv3, qp.improver_conf); }
    );
}

/// @brief Параллельный пересчёт четырёх параметров с выбранной стратегией разбиения задач
template <typename Partitioner>
void parallel_for_step(PipeQAdvection& model, double dt, quad_params_t& qp,
    double bv0, double bv1, double bv2, double bv3, Partitioner& partitioner)
{
    single_param_data_t* params[4] = { &qp.density, &qp.density_conf, &qp.improver, &qp.improver_conf };
    const double         bvs[4]    = { bv0, bv1, bv2, bv3 };
    tbb::parallel_for(
        tbb::blocked_range<size_t>(0, 4),
        [&](const tbb::blocked_range<size_t>& r) {
            for (size_t i = r.begin(); i < r.end(); ++i)
                do_param_step(model, dt, bvs[i], *params[i]);
        },
        partitioner
    );
}

/// @brief Параллельный пересчёт четырёх параметров через пул задач TBB
inline void task_group_step(PipeQAdvection& model, double dt, quad_params_t& qp,
    double bv0, double bv1, double bv2, double bv3)
{
    tbb::task_group tg;
    tg.run([&] { do_param_step(model, dt, bv0, qp.density); });
    tg.run([&] { do_param_step(model, dt, bv1, qp.density_conf); });
    tg.run([&] { do_param_step(model, dt, bv2, qp.improver); });
    tg.run_and_wait([&] { do_param_step(model, dt, bv3, qp.improver_conf); });
}

/// @brief Последовательный шаг с граничными условиями по умолчанию
inline void sequential_step(PipeQAdvection& model, double dt, quad_params_t& qp)
{
    sequential_step(model, dt, qp, RHO_IN, CONF_IN, IMPR_IN, IMPR_CONF_IN);
}

/// @brief Параллельный шаг с граничными условиями по умолчанию
inline void parallel_invoke_step(PipeQAdvection& model, double dt, quad_params_t& qp)
{
    parallel_invoke_step(model, dt, qp, RHO_IN, CONF_IN, IMPR_IN, IMPR_CONF_IN);
}

/// @brief Параллельный шаг с граничными условиями по умолчанию и выбранным партиционером
template <typename Partitioner>
inline void parallel_for_step(PipeQAdvection& model, double dt, quad_params_t& qp, Partitioner& partitioner)
{
    parallel_for_step(model, dt, qp, RHO_IN, CONF_IN, IMPR_IN, IMPR_CONF_IN, partitioner);
}

/// @brief Параллельный шаг через task_group с граничными условиями по умолчанию
inline void task_group_step(PipeQAdvection& model, double dt, quad_params_t& qp)
{
    task_group_step(model, dt, qp, RHO_IN, CONF_IN, IMPR_IN, IMPR_CONF_IN);
}

} // namespace bench6

/// @brief Тесты Bench6: плоский буфер [U_prev | U_curr | F] в одном аллоке

/// @brief Последовательный расчёт 4 параметров с плоским размещением [U_prev|U_curr|F] в одном
/// векторе на параметр — baseline для bench6.
TEST(TBB_Bench6_FlatStorage, DISABLED_Sequential_Baseline)
{
    using namespace quickest_parallel_bench;
    using namespace bench6;

    pipe_properties_t pipe_solo = pipe_properties_t::build_simple_pipe(
        { .length = PIPE_LENGTH, .dx = PIPE_DX, .diameter = PIPE_DIAM });
    std::vector<double> flow_rate_solo(pipe_solo.profile.get_point_count(), FLOW);
    PipeQAdvection model_solo(pipe_solo, flow_rate_solo);
    const double dt = calc_time_step_by_Courant(model_solo, 1.0);
    quad_params_t qp(pipe_solo.profile.get_point_count());

    /// @brief Типовой сценарий pde_solvers
    const auto solo_start = std::chrono::steady_clock::now();
    for (size_t step = 0; step < STEP_COUNT; ++step) {
        sequential_step(model_solo, dt, qp);
        qp.advance();
    }
    print_seq("sequential-4p(flat_storage)", std::chrono::duration<double>(std::chrono::steady_clock::now() - solo_start).count());
}

/// @brief 4 параметра параллельно через parallel_invoke без арены, плоское хранилище.
TEST(TBB_Bench6_FlatStorage, DISABLED_ParallelInvoke)
{
    using namespace quickest_parallel_bench;
    using namespace bench6;

    pipe_properties_t pipe_seq = pipe_properties_t::build_simple_pipe(
        { .length = PIPE_LENGTH, .dx = PIPE_DX, .diameter = PIPE_DIAM });
    std::vector<double> flow_rate_seq(pipe_seq.profile.get_point_count(), FLOW);
    PipeQAdvection model_seq(pipe_seq, flow_rate_seq);
    const double dt = calc_time_step_by_Courant(model_seq, 1.0);
    quad_params_t quad_seq(pipe_seq.profile.get_point_count());
    const auto seq_start = std::chrono::steady_clock::now();
    for (size_t step = 0; step < STEP_COUNT; ++step) {
        sequential_step(model_seq, dt, quad_seq);
        quad_seq.advance();
    }
    const double seq_elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - seq_start).count();

    /// @brief Типовой сценарий pde_solvers
    pipe_properties_t pipe_par = pipe_properties_t::build_simple_pipe(
        { .length = PIPE_LENGTH, .dx = PIPE_DX, .diameter = PIPE_DIAM });
    std::vector<double> flow_rate_par(pipe_par.profile.get_point_count(), FLOW);
    PipeQAdvection model_par(pipe_par, flow_rate_par);
    quad_params_t quad_par(pipe_par.profile.get_point_count());
    const auto par_start = std::chrono::steady_clock::now();
    for (size_t step = 0; step < STEP_COUNT; ++step) {
        parallel_invoke_step(model_par, dt, quad_par);
        quad_par.advance();
    }
    print_result("parallel_invoke(flat_storage) pde_solvers",
        seq_elapsed, std::chrono::duration<double>(std::chrono::steady_clock::now() - par_start).count());
}

/// @brief 4 параметра параллельно через parallel_for с auto_partitioner без арены, плоское хранилище.
TEST(TBB_Bench6_FlatStorage, DISABLED_ParallelFor_Auto)
{
    using namespace quickest_parallel_bench;
    using namespace bench6;

    pipe_properties_t pipe_seq = pipe_properties_t::build_simple_pipe(
        { .length = PIPE_LENGTH, .dx = PIPE_DX, .diameter = PIPE_DIAM });
    std::vector<double> flow_rate_seq(pipe_seq.profile.get_point_count(), FLOW);
    PipeQAdvection model_seq(pipe_seq, flow_rate_seq);
    const double dt = calc_time_step_by_Courant(model_seq, 1.0);
    quad_params_t quad_seq(pipe_seq.profile.get_point_count());
    const auto seq_start = std::chrono::steady_clock::now();
    for (size_t step = 0; step < STEP_COUNT; ++step) {
        sequential_step(model_seq, dt, quad_seq);
        quad_seq.advance();
    }
    const double seq_elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - seq_start).count();

    pipe_properties_t pipe_par = pipe_properties_t::build_simple_pipe(
        { .length = PIPE_LENGTH, .dx = PIPE_DX, .diameter = PIPE_DIAM });
    std::vector<double> flow_rate_par(pipe_par.profile.get_point_count(), FLOW);
    PipeQAdvection model_par(pipe_par, flow_rate_par);
    quad_params_t quad_par(pipe_par.profile.get_point_count());
    tbb::auto_partitioner auto_part;
    const auto par_start = std::chrono::steady_clock::now();
    for (size_t step = 0; step < STEP_COUNT; ++step) {
        parallel_for_step(model_par, dt, quad_par, auto_part);
        quad_par.advance();
    }
    print_result("parallel_for(auto,flat_storage)",
        seq_elapsed, std::chrono::duration<double>(std::chrono::steady_clock::now() - par_start).count());
}

/// @brief 4 параметра параллельно через parallel_for с static_partitioner без арены, плоское хранилище.
TEST(TBB_Bench6_FlatStorage, DISABLED_ParallelFor_Static)
{
    using namespace quickest_parallel_bench;
    using namespace bench6;

    pipe_properties_t pipe_seq = pipe_properties_t::build_simple_pipe(
        { .length = PIPE_LENGTH, .dx = PIPE_DX, .diameter = PIPE_DIAM });
    std::vector<double> flow_rate_seq(pipe_seq.profile.get_point_count(), FLOW);
    PipeQAdvection model_seq(pipe_seq, flow_rate_seq);
    const double dt = calc_time_step_by_Courant(model_seq, 1.0);
    quad_params_t quad_seq(pipe_seq.profile.get_point_count());
    const auto seq_start = std::chrono::steady_clock::now();
    for (size_t step = 0; step < STEP_COUNT; ++step) {
        sequential_step(model_seq, dt, quad_seq);
        quad_seq.advance();
    }
    const double seq_elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - seq_start).count();

    pipe_properties_t pipe_par = pipe_properties_t::build_simple_pipe(
        { .length = PIPE_LENGTH, .dx = PIPE_DX, .diameter = PIPE_DIAM });
    std::vector<double> flow_rate_par(pipe_par.profile.get_point_count(), FLOW);
    PipeQAdvection model_par(pipe_par, flow_rate_par);
    quad_params_t quad_par(pipe_par.profile.get_point_count());
    tbb::static_partitioner static_part;
    const auto par_start = std::chrono::steady_clock::now();
    for (size_t step = 0; step < STEP_COUNT; ++step) {
        parallel_for_step(model_par, dt, quad_par, static_part);
        quad_par.advance();
    }
    print_result("parallel_for(static,flat_storage)",
        seq_elapsed, std::chrono::duration<double>(std::chrono::steady_clock::now() - par_start).count());
}

/// @brief 4 параметра параллельно через parallel_for с affinity_partitioner без арены, плоское хранилище.
/// Два варианта: партиционер вне цикла (pde_solvers) и внутри (graph_solvers).
TEST(TBB_Bench6_FlatStorage, DISABLED_ParallelFor_Affinity)
{
    using namespace quickest_parallel_bench;
    using namespace bench6;

    pipe_properties_t pipe_seq = pipe_properties_t::build_simple_pipe(
        { .length = PIPE_LENGTH, .dx = PIPE_DX, .diameter = PIPE_DIAM });
    std::vector<double> flow_rate_seq(pipe_seq.profile.get_point_count(), FLOW);
    PipeQAdvection model_seq(pipe_seq, flow_rate_seq);
    const double dt = calc_time_step_by_Courant(model_seq, 1.0);
    quad_params_t quad_seq(pipe_seq.profile.get_point_count());
    const auto seq_start = std::chrono::steady_clock::now();
    for (size_t step = 0; step < STEP_COUNT; ++step) {
        sequential_step(model_seq, dt, quad_seq);
        quad_seq.advance();
    }
    const double seq_elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - seq_start).count();

    /// @brief Типовой сценарий pde_solvers: affinity_partitioner живёт вне цикла
    pipe_properties_t pipe_par = pipe_properties_t::build_simple_pipe(
        { .length = PIPE_LENGTH, .dx = PIPE_DX, .diameter = PIPE_DIAM });
    std::vector<double> flow_rate_par(pipe_par.profile.get_point_count(), FLOW);
    PipeQAdvection model_par(pipe_par, flow_rate_par);
    quad_params_t quad_par(pipe_par.profile.get_point_count());
    tbb::affinity_partitioner affinity_part;
    const auto par_start = std::chrono::steady_clock::now();
    for (size_t step = 0; step < STEP_COUNT; ++step) {
        parallel_for_step(model_par, dt, quad_par, affinity_part);
        quad_par.advance();
    }
    print_result("parallel_for(affinity,flat_storage) pde_solvers",
        seq_elapsed, std::chrono::duration<double>(std::chrono::steady_clock::now() - par_start).count());

    /// @brief Типовой сценарий graph_solvers: affinity_partitioner создаётся заново на каждом шаге
    pipe_properties_t pipe_graph = pipe_properties_t::build_simple_pipe(
        { .length = PIPE_LENGTH, .dx = PIPE_DX, .diameter = PIPE_DIAM });
    std::vector<double> flow_rate_graph(pipe_graph.profile.get_point_count(), FLOW);
    PipeQAdvection model_graph(pipe_graph, flow_rate_graph);
    quad_params_t quad_graph(pipe_graph.profile.get_point_count());
    const auto graph_start = std::chrono::steady_clock::now();
    for (size_t step = 0; step < STEP_COUNT; ++step) {
        tbb::affinity_partitioner affinity_part_local;
        parallel_for_step(model_graph, dt, quad_graph, affinity_part_local);
        quad_graph.advance();
    }
    print_result("parallel_for(affinity,flat_storage) graph_solvers",
        seq_elapsed, std::chrono::duration<double>(std::chrono::steady_clock::now() - graph_start).count());
}

/// @brief 4 параметра параллельно через task_group + run_and_wait без арены, плоское хранилище.
TEST(TBB_Bench6_FlatStorage, DISABLED_TaskGroup)
{
    using namespace quickest_parallel_bench;
    using namespace bench6;

    pipe_properties_t pipe_seq = pipe_properties_t::build_simple_pipe(
        { .length = PIPE_LENGTH, .dx = PIPE_DX, .diameter = PIPE_DIAM });
    std::vector<double> flow_rate_seq(pipe_seq.profile.get_point_count(), FLOW);
    PipeQAdvection model_seq(pipe_seq, flow_rate_seq);
    const double dt = calc_time_step_by_Courant(model_seq, 1.0);
    quad_params_t quad_seq(pipe_seq.profile.get_point_count());
    const auto seq_start = std::chrono::steady_clock::now();
    for (size_t step = 0; step < STEP_COUNT; ++step) {
        sequential_step(model_seq, dt, quad_seq);
        quad_seq.advance();
    }
    const double seq_elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - seq_start).count();

    pipe_properties_t pipe_par = pipe_properties_t::build_simple_pipe(
        { .length = PIPE_LENGTH, .dx = PIPE_DX, .diameter = PIPE_DIAM });
    std::vector<double> flow_rate_par(pipe_par.profile.get_point_count(), FLOW);
    PipeQAdvection model_par(pipe_par, flow_rate_par);
    quad_params_t quad_par(pipe_par.profile.get_point_count());
    const auto par_start = std::chrono::steady_clock::now();
    for (size_t step = 0; step < STEP_COUNT; ++step) {
        task_group_step(model_par, dt, quad_par);
        quad_par.advance();
    }
    print_result("task_group(flat_storage)",
        seq_elapsed, std::chrono::duration<double>(std::chrono::steady_clock::now() - par_start).count());
}

/// @brief 4 параметра параллельно: arena(4) + static_partitioner + плоское хранилище —
/// лучший результат bench6 (6993 ms, 2.04×).
/// Два варианта: arena вне цикла (pde_solvers) и внутри (graph_solvers).
TEST(TBB_Bench6_FlatStorage, DISABLED_Arena4_ParallelFor_Static)
{
    using namespace quickest_parallel_bench;
    using namespace bench6;

    pipe_properties_t pipe_seq = pipe_properties_t::build_simple_pipe(
        { .length = PIPE_LENGTH, .dx = PIPE_DX, .diameter = PIPE_DIAM });
    std::vector<double> flow_rate_seq(pipe_seq.profile.get_point_count(), FLOW);
    PipeQAdvection model_seq(pipe_seq, flow_rate_seq);
    const double dt = calc_time_step_by_Courant(model_seq, 1.0);
    quad_params_t quad_seq(pipe_seq.profile.get_point_count());
    const auto seq_start = std::chrono::steady_clock::now();
    for (size_t step = 0; step < STEP_COUNT; ++step) {
        sequential_step(model_seq, dt, quad_seq);
        quad_seq.advance();
    }
    const double seq_elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - seq_start).count();

    /// @brief Типовой сценарий pde_solvers: arena живёт вне цикла
    tbb::task_arena arena(4);
    pipe_properties_t pipe_par = pipe_properties_t::build_simple_pipe(
        { .length = PIPE_LENGTH, .dx = PIPE_DX, .diameter = PIPE_DIAM });
    std::vector<double> flow_rate_par(pipe_par.profile.get_point_count(), FLOW);
    PipeQAdvection model_par(pipe_par, flow_rate_par);
    quad_params_t quad_par(pipe_par.profile.get_point_count());
    tbb::static_partitioner static_part;
    const auto par_start = std::chrono::steady_clock::now();
    for (size_t step = 0; step < STEP_COUNT; ++step) {
        arena.execute([&] { parallel_for_step(model_par, dt, quad_par, static_part); });
        quad_par.advance();
    }
    print_result("arena(4)+parallel_for(static,flat_storage) pde_solvers",
        seq_elapsed, std::chrono::duration<double>(std::chrono::steady_clock::now() - par_start).count());

    /// @brief Типовой сценарий graph_solvers: arena создаётся на каждом шаге
    pipe_properties_t pipe_graph = pipe_properties_t::build_simple_pipe(
        { .length = PIPE_LENGTH, .dx = PIPE_DX, .diameter = PIPE_DIAM });
    std::vector<double> flow_rate_graph(pipe_graph.profile.get_point_count(), FLOW);
    PipeQAdvection model_graph(pipe_graph, flow_rate_graph);
    quad_params_t quad_graph(pipe_graph.profile.get_point_count());
    const auto graph_start = std::chrono::steady_clock::now();
    for (size_t step = 0; step < STEP_COUNT; ++step) {
        tbb::task_arena arena_local(4);
        tbb::static_partitioner static_part_local;
        arena_local.execute([&] { parallel_for_step(model_graph, dt, quad_graph, static_part_local); });
        quad_graph.advance();
    }
    print_result("arena(4)+parallel_for(static,flat_storage) graph_solvers",
        seq_elapsed, std::chrono::duration<double>(std::chrono::steady_clock::now() - graph_start).count());
}

/// @brief 4 параметра параллельно: arena(4) + affinity_partitioner + плоское хранилище.
/// Два варианта: arena и партиционер вне цикла (pde_solvers) и внутри (graph_solvers).
TEST(TBB_Bench6_FlatStorage, DISABLED_Arena4_ParallelFor_Affinity)
{
    using namespace quickest_parallel_bench;
    using namespace bench6;

    pipe_properties_t pipe_seq = pipe_properties_t::build_simple_pipe(
        { .length = PIPE_LENGTH, .dx = PIPE_DX, .diameter = PIPE_DIAM });
    std::vector<double> flow_rate_seq(pipe_seq.profile.get_point_count(), FLOW);
    PipeQAdvection model_seq(pipe_seq, flow_rate_seq);
    const double dt = calc_time_step_by_Courant(model_seq, 1.0);
    quad_params_t quad_seq(pipe_seq.profile.get_point_count());
    const auto seq_start = std::chrono::steady_clock::now();
    for (size_t step = 0; step < STEP_COUNT; ++step) {
        sequential_step(model_seq, dt, quad_seq);
        quad_seq.advance();
    }
    const double seq_elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - seq_start).count();

    /// @brief Типовой сценарий pde_solvers: arena + affinity_partitioner живут вне цикла
    tbb::task_arena arena(4);
    pipe_properties_t pipe_par = pipe_properties_t::build_simple_pipe(
        { .length = PIPE_LENGTH, .dx = PIPE_DX, .diameter = PIPE_DIAM });
    std::vector<double> flow_rate_par(pipe_par.profile.get_point_count(), FLOW);
    PipeQAdvection model_par(pipe_par, flow_rate_par);
    quad_params_t quad_par(pipe_par.profile.get_point_count());
    tbb::affinity_partitioner affinity_part;
    const auto par_start = std::chrono::steady_clock::now();
    for (size_t step = 0; step < STEP_COUNT; ++step) {
        arena.execute([&] { parallel_for_step(model_par, dt, quad_par, affinity_part); });
        quad_par.advance();
    }
    print_result("arena(4)+parallel_for(affinity,flat_storage) pde_solvers",
        seq_elapsed, std::chrono::duration<double>(std::chrono::steady_clock::now() - par_start).count());

    /// @brief Типовой сценарий graph_solvers: arena и affinity_partitioner создаются на каждом шаге
    pipe_properties_t pipe_graph = pipe_properties_t::build_simple_pipe(
        { .length = PIPE_LENGTH, .dx = PIPE_DX, .diameter = PIPE_DIAM });
    std::vector<double> flow_rate_graph(pipe_graph.profile.get_point_count(), FLOW);
    PipeQAdvection model_graph(pipe_graph, flow_rate_graph);
    quad_params_t quad_graph(pipe_graph.profile.get_point_count());
    const auto graph_start = std::chrono::steady_clock::now();
    for (size_t step = 0; step < STEP_COUNT; ++step) {
        tbb::task_arena arena_local(4);
        tbb::affinity_partitioner affinity_part_local;
        arena_local.execute([&] { parallel_for_step(model_graph, dt, quad_graph, affinity_part_local); });
        quad_graph.advance();
    }
    print_result("arena(4)+parallel_for(affinity,flat_storage) graph_solvers",
        seq_elapsed, std::chrono::duration<double>(std::chrono::steady_clock::now() - graph_start).count());
}

/// @brief 4 параметра параллельно: arena(4) + parallel_invoke + плоское хранилище.
/// Два варианта: arena вне цикла (pde_solvers) и внутри (graph_solvers).
TEST(TBB_Bench6_FlatStorage, DISABLED_Arena4_ParallelInvoke)
{
    using namespace quickest_parallel_bench;
    using namespace bench6;

    pipe_properties_t pipe_seq = pipe_properties_t::build_simple_pipe(
        { .length = PIPE_LENGTH, .dx = PIPE_DX, .diameter = PIPE_DIAM });
    std::vector<double> flow_rate_seq(pipe_seq.profile.get_point_count(), FLOW);
    PipeQAdvection model_seq(pipe_seq, flow_rate_seq);
    const double dt = calc_time_step_by_Courant(model_seq, 1.0);
    quad_params_t quad_seq(pipe_seq.profile.get_point_count());
    const auto seq_start = std::chrono::steady_clock::now();
    for (size_t step = 0; step < STEP_COUNT; ++step) {
        sequential_step(model_seq, dt, quad_seq);
        quad_seq.advance();
    }
    const double seq_elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - seq_start).count();

    /// @brief Типовой сценарий pde_solvers: arena живёт вне цикла
    tbb::task_arena arena(4);
    pipe_properties_t pipe_par = pipe_properties_t::build_simple_pipe(
        { .length = PIPE_LENGTH, .dx = PIPE_DX, .diameter = PIPE_DIAM });
    std::vector<double> flow_rate_par(pipe_par.profile.get_point_count(), FLOW);
    PipeQAdvection model_par(pipe_par, flow_rate_par);
    quad_params_t quad_par(pipe_par.profile.get_point_count());
    const auto par_start = std::chrono::steady_clock::now();
    for (size_t step = 0; step < STEP_COUNT; ++step) {
        arena.execute([&] { parallel_invoke_step(model_par, dt, quad_par); });
        quad_par.advance();
    }
    print_result("arena(4)+parallel_invoke(flat_storage) pde_solvers",
        seq_elapsed, std::chrono::duration<double>(std::chrono::steady_clock::now() - par_start).count());

    /// @brief Типовой сценарий graph_solvers: arena создаётся на каждом шаге
    pipe_properties_t pipe_graph = pipe_properties_t::build_simple_pipe(
        { .length = PIPE_LENGTH, .dx = PIPE_DX, .diameter = PIPE_DIAM });
    std::vector<double> flow_rate_graph(pipe_graph.profile.get_point_count(), FLOW);
    PipeQAdvection model_graph(pipe_graph, flow_rate_graph);
    quad_params_t quad_graph(pipe_graph.profile.get_point_count());
    const auto graph_start = std::chrono::steady_clock::now();
    for (size_t step = 0; step < STEP_COUNT; ++step) {
        tbb::task_arena arena_local(4);
        arena_local.execute([&] { parallel_invoke_step(model_graph, dt, quad_graph); });
        quad_graph.advance();
    }
    print_result("arena(4)+parallel_invoke(flat_storage) graph_solvers",
        seq_elapsed, std::chrono::duration<double>(std::chrono::steady_clock::now() - graph_start).count());
}

/// @brief 4 параметра параллельно: arena(2) + static_partitioner + плоское хранилище.
/// Два варианта: arena вне цикла (pde_solvers) и внутри (graph_solvers).
TEST(TBB_Bench6_FlatStorage, DISABLED_Arena2_ParallelFor_Static)
{
    using namespace quickest_parallel_bench;
    using namespace bench6;

    pipe_properties_t pipe_seq = pipe_properties_t::build_simple_pipe(
        { .length = PIPE_LENGTH, .dx = PIPE_DX, .diameter = PIPE_DIAM });
    std::vector<double> flow_rate_seq(pipe_seq.profile.get_point_count(), FLOW);
    PipeQAdvection model_seq(pipe_seq, flow_rate_seq);
    const double dt = calc_time_step_by_Courant(model_seq, 1.0);
    quad_params_t quad_seq(pipe_seq.profile.get_point_count());
    const auto seq_start = std::chrono::steady_clock::now();
    for (size_t step = 0; step < STEP_COUNT; ++step) {
        sequential_step(model_seq, dt, quad_seq);
        quad_seq.advance();
    }
    const double seq_elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - seq_start).count();

    /// @brief Типовой сценарий pde_solvers: arena живёт вне цикла
    tbb::task_arena arena(2);
    pipe_properties_t pipe_par = pipe_properties_t::build_simple_pipe(
        { .length = PIPE_LENGTH, .dx = PIPE_DX, .diameter = PIPE_DIAM });
    std::vector<double> flow_rate_par(pipe_par.profile.get_point_count(), FLOW);
    PipeQAdvection model_par(pipe_par, flow_rate_par);
    quad_params_t quad_par(pipe_par.profile.get_point_count());
    tbb::static_partitioner static_part;
    const auto par_start = std::chrono::steady_clock::now();
    for (size_t step = 0; step < STEP_COUNT; ++step) {
        arena.execute([&] { parallel_for_step(model_par, dt, quad_par, static_part); });
        quad_par.advance();
    }
    print_result("arena(2)+parallel_for(static,flat_storage) pde_solvers",
        seq_elapsed, std::chrono::duration<double>(std::chrono::steady_clock::now() - par_start).count());

    /// @brief Типовой сценарий graph_solvers: arena создаётся на каждом шаге
    pipe_properties_t pipe_graph = pipe_properties_t::build_simple_pipe(
        { .length = PIPE_LENGTH, .dx = PIPE_DX, .diameter = PIPE_DIAM });
    std::vector<double> flow_rate_graph(pipe_graph.profile.get_point_count(), FLOW);
    PipeQAdvection model_graph(pipe_graph, flow_rate_graph);
    quad_params_t quad_graph(pipe_graph.profile.get_point_count());
    const auto graph_start = std::chrono::steady_clock::now();
    for (size_t step = 0; step < STEP_COUNT; ++step) {
        tbb::task_arena arena_local(2);
        tbb::static_partitioner static_part_local;
        arena_local.execute([&] { parallel_for_step(model_graph, dt, quad_graph, static_part_local); });
        quad_graph.advance();
    }
    print_result("arena(2)+parallel_for(static,flat_storage) graph_solvers",
        seq_elapsed, std::chrono::duration<double>(std::chrono::steady_clock::now() - graph_start).count());
}
