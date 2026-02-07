#pragma once


namespace pde_solvers {

/// @brief Интерфейс для транспортных солверов трубы (движение партий)
class pipe_solver_transport_interface_t {
public:
    /// @brief Выполнение шага эндогенного расчета
    /// @param dt Временной шаг
    /// @param volumetric_flow Объемный расход
    /// @param boundaries Эндогенные свойства на входе трубы
    virtual void transport_step(double dt, double volumetric_flow, const pde_solvers::endogenous_values_t& boundaries) = 0;
    /// @brief Получение эндогенных значений на выходе трубы
    /// @param volumetric_flow Объемный расход для определения направления движения партий
    virtual pde_solvers::endogenous_values_t get_endogenous_output(double volumetric_flow) const = 0;
};

/// @brief Интерфейс для гидравлического расчета трубы
class pipe_solver_hydro_interface_t {
public:
    /// @brief Решение гидравлической задачи PP
    virtual double hydro_solve_PP(double pressure_input, double pressure_output) = 0;
    /// @brief Решение гидравлической задачи QP/PQ (заданы расход и давление на одной границе)
    /// @param volumetric_flow Объемный расход
    /// @param bound_pressure Граничное давление (граница задаётся solve_direction)
    /// @param solve_direction Если >0: bound_pressure — давление на входе (PQ), возвращается давление на выходе.
    ///        Если <0: bound_pressure — давление на выходе (QP), возвращается давление на входе.
    virtual double hydro_solve_QP(double volumetric_flow, double bound_pressure, int solve_direction) = 0;
    /// @brief Вычисление якобиана для задачи PP (реализация по умолчанию через численное дифференцирование)
    /// TODO: Не очень понятно, как это дружит с самотеками (вероятно, никак - там вообще нужен МГГ вместо МД)
    /// @return Массив из двух элементов: [dQ/dP_in, dQ/dP_out]
    virtual std::array<double, 2> hydro_solve_PP_jacobian(double pressure_input, double pressure_output)  {
        // Вычисляем базовое решение - расход при заданных давлениях
        double Q_base = hydro_solve_PP(pressure_input, pressure_output);

        // Малое приращение для численного дифференцирования (0.1% от расхода)
        const double eps = std::max(1e-6, std::abs(Q_base) * 1e-3);

        // Вычисляем производную перепада давления по расходу dP/dQ
        // Вычисляем давление на выходе при базовом и увеличенном расходе (PQ: solve_direction = +1)
        double P_out_base = hydro_solve_QP(Q_base, pressure_input, +1);
        double P_out_plus = hydro_solve_QP(Q_base + eps, pressure_input, +1);

        // Производная давления на выходе по расходу: dP_out/dQ
        double dP_out_dQ = (P_out_plus - P_out_base) / eps;

        // Производная перепада давления по расходу: dP/dQ = -dP_out/dQ
        double dP_dQ = -dP_out_dQ;

        // Используем формулу переворота производной:
        double dQ_dP_in = 1.0 / dP_dQ;
        double dQ_dP_out = -1.0 / dP_dQ;

        return { dQ_dP_in, dQ_dP_out };
    }
};

/// @brief Объединенный интерфейс для совмещенных солвером (умеют выполнять транспортный и гидравлический расчеты)
class pipe_solver_hydrotransport_interface_t
    : public pipe_solver_hydro_interface_t
    , public pipe_solver_transport_interface_t
{

};

/// @brief Traits для солвера, реализующего гидравлический интерфейс. Используется для выявления типа солвера в тасках
template<typename T>
constexpr bool has_hydro_interface_v = std::is_base_of_v<pde_solvers::pipe_solver_hydro_interface_t, T>;

/// @brief Traits для солвера, реализующего только транспортный интерфейс (без гидравлического). 
/// Используется для выявления типа солвера в тасках, которые могут работать только с транспортным солвером
template<typename T>
constexpr bool is_transport_only_solver_v = 
    std::is_base_of_v<pde_solvers::pipe_solver_transport_interface_t, T> &&
    !std::is_base_of_v<pde_solvers::pipe_solver_hydro_interface_t, T>;

}
