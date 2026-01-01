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
    /// @brief Решение гидравлической задачи QP
    virtual double hydro_solve_QP(double volumetric_flow, double pressure_output) = 0;
    /// @brief Решение гидравлической задачи PQ
    virtual double hydro_solve_PQ(double volumetric_flow, double pressure_in) = 0;
    /// @brief Вычисление якобиана для задачи PP
    virtual std::array<double, 2> hydro_solve_PP_jacobian(double pressure_input, double pressure_output) = 0;
};

/// @brief Объединенный интерфейс для совмещенных солвером (умеют выполнять транспортный и гидравлический расчеты)
class pipe_solver_interface_t
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
