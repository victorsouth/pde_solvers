#pragma once

namespace pde_solvers {
;

/// @brief Структура, содержащая в себе краевые условия задачи PP
struct iso_nonbaro_pipe_PP_task_boundaries_t {
    /// @brief Изначальное давление на входе
    double pressure_in;
    /// @brief Изначальное давление на выходе
    double pressure_out;
    /// @brief Изначальная плотность на входе
    double density;
    /// @brief Создание структуры со значениями по умолчанию
    static iso_nonbaro_pipe_PP_task_boundaries_t default_values() {
        iso_nonbaro_pipe_PP_task_boundaries_t result;
        result.pressure_out = 0.6e6;
        result.pressure_in = 6e6;
        result.density = 850;
        return result;
    }
};




/// @brief Солвер квазистационарного гидравлического расчета для конденсатопровода
class iso_nonbaro_pipe_solver_t : public pipe_solver_hydrotransport_interface_t {
public:
    /// @brief Тип слоя
    using layer_type = iso_nonbaro_pipe_layer_t;
    /// @brief Тип буфера
    using buffer_type = ring_buffer_t<layer_type>;
    /// @brief Тип параметров трубы
    using pipe_parameters_type = iso_nonbaro_pipe_properties_t;

private:
    /// @brief Ссылка на свойства конденсатопровода
    const iso_nonbaro_pipe_properties_t& pipe;
    /// @brief Ссылка на буфер слоев
    buffer_type& buffer;

public:
    /// @brief Фиктивный констуктор для совместмости с селектором рассчитываемых свойств
    iso_nonbaro_pipe_solver_t(
        const iso_nonbaro_pipe_properties_t& pipe,
        buffer_type& buffer,
        const pde_solvers::endogenous_selector_t& endogenous_selector)
        : iso_nonbaro_pipe_solver_t(pipe, buffer)
    {
    }

    /// @brief Конструктор
    /// @param pipe Ссылка на свойства конденсатопровода
    /// @param buffer Ссылка на буфер слоев
    iso_nonbaro_pipe_solver_t(
        const iso_nonbaro_pipe_properties_t& pipe,
        buffer_type& buffer)
        : pipe(pipe)
        , buffer(buffer)
    {
    }

    /// @brief Решение гидравлической задачи PP (заданы давления на входе и выходе, найти расход)
    virtual double hydro_solve_PP(double pressure_input, double pressure_output) override {
        auto& current = buffer.current();
        double volumetric_flow_initial = 0.0;
        if (std::isfinite(current.std_volumetric_flow)) {
            volumetric_flow_initial = current.std_volumetric_flow;
        }

        rigorous_impulse_solver_PP<iso_nonbaro_impulse_equation_t>
            solver_pp(pipe, current, pressure_input, pressure_output);

        double std_vol_flow = solver_pp.solve(volumetric_flow_initial);
        return std_vol_flow;
    }

    /// @brief Решение гидравлической задачи QP/PQ (заданы расход и давление на одной границе)
    virtual double hydro_solve_QP(double volumetric_flow, double bound_pressure, int solve_direction) override {
        auto& current = buffer.current();
        return rigorous_impulse_solve_QP<iso_nonbaro_impulse_equation_t>(
            pipe, current, volumetric_flow, bound_pressure, solve_direction);
    }

    /// @brief Выполнение транспортного шага (расчет движения партий)
    virtual void transport_step(double dt, double volumetric_flow, const pde_solvers::endogenous_values_t& boundaries) override {
        if (!std::isfinite(dt)) {
            throw std::runtime_error("transport_step: dt must be finite; use transport_solve for steady-state (infinite dt)");
        }

        buffer.current().std_volumetric_flow = volumetric_flow;

        step_advection(dt, volumetric_flow, boundaries.density_std,
            pipe.wall.diameter, pipe.profile.coordinates, buffer,
            &iso_nonbaro_pipe_endogenious_layer_t::get_density_wrapper,
            &iso_nonbaro_pipe_endogenious_layer_t::get_density_confidence_wrapper);
    }

    /// @brief Транспортное решение при бесконечном dt (заполнение трубы граничными значениями)
    virtual void transport_solve(double volumetric_flow, const pde_solvers::endogenous_values_t& boundaries) override {
        buffer.current().std_volumetric_flow = volumetric_flow;

        solve_advection(boundaries.density_std, buffer,
            &iso_nonbaro_pipe_endogenious_layer_t::get_density_wrapper,
            &iso_nonbaro_pipe_endogenious_layer_t::get_density_confidence_wrapper);
    }

    ///// @brief Получает эндогенные значения на выходе трубы
    ///// @param volumetric_flow Объемный расход
    ///// @return Эндогенные значения на выходе
    pde_solvers::endogenous_values_t get_endogenous_output(double volumetric_flow) const override
    {
        pde_solvers::endogenous_values_t result;
        layer_type& current_layer = buffer.current();
        result.density_std = current_layer.density_std.get_boundary_value(volumetric_flow);
        return result;
    }

    /// @brief Геттер для текущего слоя
    iso_nonbaro_pipe_layer_t& get_current_layer() {
        return buffer.current();
    }
};


/// @brief Структура, содержащая в себе краевые условия задачи PQ
struct iso_nonbarotropic_pipe_PQ_task_boundaries_t {
    /// @brief Изначальный объемный расход
    double volumetric_flow;
    /// @brief Изначальное давление на входе
    double pressure_in;
    /// @brief Изначальная плотность на входе
    double density;
    /// @brief Создание структуры со значениями по умолчанию
    static iso_nonbarotropic_pipe_PQ_task_boundaries_t default_values() {
        iso_nonbarotropic_pipe_PQ_task_boundaries_t result;
        result.volumetric_flow = 0.2;
        result.pressure_in = 6e6;
        result.density = 850;
        return result;
    }
};

/// @brief Квазистационарная задача PQ в условиях переменной плотности
class iso_nonbarotropic_pipe_PQ_task_t {
public:
    /// @brief Тип слоя
    using layer_type = iso_nonbaro_pipe_layer_t;
    /// @brief Тип буфера
    using buffer_type = ring_buffer_t<layer_type>;
    /// @brief Тип граничных условий
    using boundaries_type = iso_nonbarotropic_pipe_PQ_task_boundaries_t;
private:
    // Модель трубы
    iso_nonbaro_pipe_properties_t pipe;
    // Создаётся буфер, тип слоя которого определяется в зависимости от типа солвера
    buffer_type buffer;
public:
    /// @brief Конструктор
    /// @param pipe Модель трубопровода
    iso_nonbarotropic_pipe_PQ_task_t(const iso_nonbaro_pipe_properties_t& pipe)
        : pipe(pipe)
        , buffer(2, pipe.profile.get_point_count())
    {
    }


    /// @brief Начальный стационарный расчёт.
    /// Ставим по всей трубе реологию из initial_conditions, делаем гидравлический расчет
    /// @param initial_conditions Начальные условия
    void solve(const boundaries_type& initial_conditions)
    {
        pde_solvers::endogenous_values_t endogenous_boundaries;
        endogenous_boundaries.density_std.value = initial_conditions.density;
        endogenous_boundaries.density_std.confidence = true;

        iso_nonbaro_pipe_solver_t solver(pipe, buffer);
        solver.transport_solve(initial_conditions.volumetric_flow, endogenous_boundaries);
        solver.hydro_solve_QP(initial_conditions.volumetric_flow, initial_conditions.pressure_in, +1);
    }

    /// @brief Рассчёт шага моделирования, включающий в себя расчёт шага движения партии и гидравлический расчёт
    /// buffer.current после вызова содержит свежерасчитанный слой
    /// @param dt временной шаг моделирования
    /// @param boundaries Краевые условие
    void step(double dt, const boundaries_type& boundaries) {
        advance();

        pde_solvers::endogenous_values_t endogenous_boundaries;
        endogenous_boundaries.density_std.value = boundaries.density;
        endogenous_boundaries.density_std.confidence = true;

        iso_nonbaro_pipe_solver_t solver(pipe, buffer);
        if (std::isfinite(dt)) {
            solver.transport_step(dt, boundaries.volumetric_flow, endogenous_boundaries);
        } else {
            solver.transport_solve(boundaries.volumetric_flow, endogenous_boundaries);
        }
        solver.hydro_solve_QP(boundaries.volumetric_flow, boundaries.pressure_in, +1);
    }

    /// @brief Сдвиг текущего слоя в буфере
    void advance()
    {
        buffer.advance(+1);
    }
    /// @brief Возвращает ссылку на буфер
    auto& get_buffer()
    {
        return buffer;
    }
    /// @brief Геттер для текущего слоя  
    iso_nonbaro_pipe_layer_t& get_current_layer() {
        return buffer.current();
    }
};

/// @brief Квазистационарная задача PP в условиях переменной плотности
class iso_nonbarotropic_pipe_PP_task_t {
public:
    /// @brief Тип слоя
    using layer_type = iso_nonbaro_pipe_layer_t;
    /// @brief Тип буфера
    using buffer_type = ring_buffer_t<layer_type>;
    /// @brief Тип граничных условий
    using boundaries_type = iso_nonbaro_pipe_PP_task_boundaries_t;
private:
    // Модель трубы
    iso_nonbaro_pipe_properties_t pipe;
    // Создаётся буфер, тип слоя которого определяется в зависимости от типа солвера
    buffer_type buffer;

public:
    /// @brief Конструктор
    /// @param pipe Модель трубопровода
    iso_nonbarotropic_pipe_PP_task_t(const iso_nonbaro_pipe_properties_t& pipe)
        : pipe(pipe)
        , buffer(2, pipe.profile.get_point_count())
    {
    }


    /// @brief Начальный стационарный расчёт.
    /// Ставим по всей трубе реологию из initial_conditions, делаем гидравлический расчет
    /// @param initial_conditions Начальные условия
    /// @param volumetric_flow_initial Начальное значение расхода для заполнения трубы
    void solve(const boundaries_type& initial_conditions, double volumetric_flow_initial = 0.2)
    {
        pde_solvers::endogenous_values_t endogenous_boundaries;
        endogenous_boundaries.density_std.value = initial_conditions.density;
        endogenous_boundaries.density_std.confidence = true;

        iso_nonbaro_pipe_solver_t solver(pipe, buffer);
        solver.transport_solve(volumetric_flow_initial, endogenous_boundaries);
        solver.hydro_solve_PP(initial_conditions.pressure_in, initial_conditions.pressure_out);
    }

    /// @brief Рассчёт шага моделирования, включающий в себя расчёт шага движения партии и гидравлический расчёт
    /// buffer.current после вызова содержит свежерасчитанный слой
    /// @param dt временной шаг моделирования
    /// @param boundaries Краевые условие
    void step(double dt, const boundaries_type& boundaries) {
        double volumetric_flow = buffer.current().std_volumetric_flow;
        advance();

        pde_solvers::endogenous_values_t endogenous_boundaries;
        endogenous_boundaries.density_std.value = boundaries.density;
        endogenous_boundaries.density_std.confidence = true;

        iso_nonbaro_pipe_solver_t solver(pipe, buffer);
        if (std::isfinite(dt)) {
            solver.transport_step(dt, volumetric_flow, endogenous_boundaries);
        } else {
            solver.transport_solve(volumetric_flow, endogenous_boundaries);
        }
        solver.hydro_solve_PP(boundaries.pressure_in, boundaries.pressure_out);
    }

    /// @brief Сдвиг текущего слоя в буфере
    void advance()
    {
        buffer.advance(+1);
    }
    /// @brief Возвращает ссылку на буфер
    auto& get_buffer()
    {
        return buffer;
    }
    /// @brief Геттер для текущего слоя  
    iso_nonbaro_pipe_layer_t& get_current_layer() {
        return buffer.current();
    }
};

/// @brief Солвер квазистационарного гидравлического расчета для конденсатопровода
class iso_nonbaro_improver_pipe_solver_t : public pipe_solver_hydrotransport_interface_t {
public:
    /// @brief Тип слоя
    using layer_type = iso_nonbaro_improver_pipe_layer_t;
    /// @brief Тип буфера
    using buffer_type = ring_buffer_t<layer_type>;
    /// @brief Тип параметров трубы
    using pipe_parameters_type = iso_nonbaro_improver_pipe_properties_t;

private:
    /// @brief Ссылка на свойства конденсатопровода
    const pipe_parameters_type& pipe;
    /// @brief Ссылка на буфер слоев
    buffer_type& buffer;

public:
    /// @brief Фиктивный констуктор для совместмости с селектором рассчитываемых свойств
    iso_nonbaro_improver_pipe_solver_t(
        const iso_nonbaro_improver_pipe_properties_t& pipe,
        buffer_type& buffer,
        const pde_solvers::endogenous_selector_t& endogenous_selector)
        : iso_nonbaro_improver_pipe_solver_t(pipe, buffer)
    {
    }

    /// @brief Конструктор
    /// @param pipe Ссылка на свойства конденсатопровода
    /// @param buffer Ссылка на буфер слоев
    iso_nonbaro_improver_pipe_solver_t(
        const iso_nonbaro_improver_pipe_properties_t& pipe,
        buffer_type& buffer)
        : pipe(pipe)
        , buffer(buffer)
    {
    }

    /// @brief Решение гидравлической задачи PP (заданы давления на входе и выходе, найти расход)
    virtual double hydro_solve_PP(double pressure_input, double pressure_output) override {
        auto& current = buffer.current();
        double volumetric_flow_initial = 0.0;
        if (std::isfinite(current.std_volumetric_flow)) {
            volumetric_flow_initial = current.std_volumetric_flow;
        }

        rigorous_impulse_solver_PP<iso_nonbaro_improver_impulse_equation_t>
            solver_pp(pipe, current, pressure_input, pressure_output);
        return solver_pp.solve(volumetric_flow_initial);
    }

    /// @brief Решение гидравлической задачи QP/PQ (заданы расход и давление на одной границе)
    virtual double hydro_solve_QP(double volumetric_flow, double bound_pressure, int solve_direction) override {
        auto& current = buffer.current();

        return rigorous_impulse_solve_QP<iso_nonbaro_improver_impulse_equation_t>(
            pipe, current, volumetric_flow, bound_pressure, solve_direction);
    }

    /// @brief Выполнение транспортного шага (расчет движения партий)
    virtual void transport_step(double dt, double volumetric_flow, const pde_solvers::endogenous_values_t& boundaries) override {
        if (!std::isfinite(dt)) {
            throw std::runtime_error("transport_step: dt must be finite; use transport_solve for steady-state (infinite dt)");
        }

        buffer.current().std_volumetric_flow = volumetric_flow;

        step_advection(dt, volumetric_flow, boundaries.density_std,
            pipe.wall.diameter, pipe.profile.coordinates, buffer,
            &iso_nonbaro_improver_pipe_endogenious_layer_t::get_density_wrapper,
            &iso_nonbaro_improver_pipe_endogenious_layer_t::get_density_confidence_wrapper);
        step_advection(dt, volumetric_flow, boundaries.improver,
            pipe.wall.diameter, pipe.profile.coordinates, buffer,
            &iso_nonbaro_improver_pipe_endogenious_layer_t::get_improver_concentration_wrapper,
            &iso_nonbaro_improver_pipe_endogenious_layer_t::get_improver_concentration_confidence_wrapper);
    }

    /// @brief Транспортное решение при бесконечном dt (заполнение трубы граничными значениями)
    virtual void transport_solve(double volumetric_flow, const pde_solvers::endogenous_values_t& boundaries) override {
        buffer.current().std_volumetric_flow = volumetric_flow;

        solve_advection(boundaries.density_std, buffer,
            &iso_nonbaro_improver_pipe_endogenious_layer_t::get_density_wrapper,
            &iso_nonbaro_improver_pipe_endogenious_layer_t::get_density_confidence_wrapper);
        solve_advection(boundaries.improver, buffer,
            &iso_nonbaro_improver_pipe_endogenious_layer_t::get_improver_concentration_wrapper,
            &iso_nonbaro_improver_pipe_endogenious_layer_t::get_improver_concentration_confidence_wrapper);
    }

    ///// @brief Получает эндогенные значения на выходе трубы
    ///// @param volumetric_flow Объемный расход
    ///// @return Эндогенные значения на выходе
    pde_solvers::endogenous_values_t get_endogenous_output(double volumetric_flow) const override
    {
        pde_solvers::endogenous_values_t result;
        layer_type& current_layer = buffer.current();
        result.density_std = current_layer.density_std.get_boundary_value(volumetric_flow);
        result.improver = current_layer.improver_concentration.get_boundary_value(volumetric_flow);
        return result;
    }

    /// @brief Геттер для текущего слоя
    iso_nonbaro_improver_pipe_layer_t& get_current_layer() {
        return buffer.current();
    }
};


}

