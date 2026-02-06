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
class iso_nonbarotropic_pipe_solver_t : public pipe_solver_hydrotransport_interface_t {
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
    iso_nonbarotropic_pipe_solver_t(
        const iso_nonbaro_pipe_properties_t& pipe,
        buffer_type& buffer,
        const pde_solvers::endogenous_selector_t& endogenous_selector)
        : iso_nonbarotropic_pipe_solver_t(pipe, buffer)
    {
    }

    /// @brief Конструктор
    /// @param pipe Ссылка на свойства конденсатопровода
    /// @param buffer Ссылка на буфер слоев
    iso_nonbarotropic_pipe_solver_t(
        const iso_nonbaro_pipe_properties_t& pipe,
        buffer_type& buffer)
        : pipe(pipe)
        , buffer(buffer)
    {
    }

    /// @brief Решение гидравлической задачи PP (заданы давления на входе и выходе, найти расход)
    virtual double hydro_solve_PP(double pressure_input, double pressure_output) override {
        auto& current = buffer.current();
        double volumetric_flow_initial = current.std_volumetric_flow;

        rigorous_impulse_solver_PP<iso_nonbaro_impulse_equation_t>
            solver_pp(pipe, current, pressure_input, pressure_output);

        double std_vol_flow = solver_pp.solve(volumetric_flow_initial);
        return std_vol_flow;
    }

    /// @brief Решение гидравлической задачи QP
    virtual double hydro_solve_QP(double volumetric_flow, double pressure_output) override {
        auto& current = buffer.current();

        // Рассчитываем профиль давления методом Эйлера в обратном направлении (от выхода ко входу)
        std::vector<double>& p_profile = current.pressure;
        int euler_direction = -1;
        iso_nonbaro_impulse_equation_t pipeModel(pipe, current, volumetric_flow, euler_direction);
        solve_euler<1>(pipeModel, euler_direction, pressure_output, &p_profile);

        // Обновляем расход в текущем слое
        current.std_volumetric_flow = volumetric_flow;

        // Возвращаем давление на входе (первый элемент профиля)
        return p_profile.front();
    }

    /// @brief Решение гидравлической задачи PQ
    /// @return Давление на выходе, Па
    virtual double hydro_solve_PQ(double volumetric_flow, double pressure_in) override {
        auto& current = buffer.current();

        // Проверяем наличие данных о плотности
        if (current.density_std.value.empty()) {
            throw std::runtime_error("iso_nonbarotropic_pipe_solver_t::hydro_solve_PQ: density profile is empty");
        }

        // Рассчитываем профиль давления методом Эйлера
        std::vector<double>& p_profile = current.pressure;
        int euler_direction = +1; // Задаем направление для Эйлера
        iso_nonbaro_impulse_equation_t pipeModel(pipe, current, volumetric_flow, euler_direction);
        solve_euler<1>(pipeModel, euler_direction, pressure_in, &p_profile);

        // Обновляем расход в текущем слое
        current.std_volumetric_flow = volumetric_flow;

        // Возвращаем давление на выходе (последний элемент профиля)
        return p_profile.back();
    }

    /// @brief Выполнение транспортного шага (расчет движения партий)
    virtual void transport_step(double dt, double volumetric_flow, const pde_solvers::endogenous_values_t& boundaries) override {
        buffer.current().std_volumetric_flow = volumetric_flow;

        if (std::isinf(dt)) {
            auto value_buffer_wrapper = buffer.get_buffer_wrapper(
                &iso_nonbaro_pipe_endogenious_layer_t::get_density_wrapper);
            std::vector<double>& val = value_buffer_wrapper.current().vars;
            std::fill(val.begin(), val.end(), boundaries.density_std.value);

            auto confidence_buffer_wrapper = buffer.get_buffer_wrapper(
                &iso_nonbaro_pipe_endogenious_layer_t::get_density_confidence_wrapper);
            std::vector<double>& conf = confidence_buffer_wrapper.current().vars;
            std::fill(conf.begin(), conf.end(), boundaries.density_std.confidence);
        }
        else {
            // считаем партии с помощью QUICKEST-ULTIMATE
            pde_solvers::pipe_advection_pde_t advection_model(pde_solvers::circle_area(pipe.wall.diameter),
                volumetric_flow, pipe.profile.coordinates);
            {
                auto buffer_wrapper = buffer.get_buffer_wrapper(
                    &iso_nonbaro_pipe_endogenious_layer_t::get_density_wrapper);
                pde_solvers::quickest_ultimate_fv_solver solver(advection_model, buffer_wrapper);
                solver.step(dt, boundaries.density_std.value, boundaries.density_std.value);
            }
            {
                auto buffer_wrapper = buffer.get_buffer_wrapper(
                    &iso_nonbaro_pipe_endogenious_layer_t::get_density_confidence_wrapper);
                pde_solvers::quickest_ultimate_fv_solver solver(advection_model, buffer_wrapper);
                solver.step(dt, boundaries.density_std.confidence, boundaries.density_std.confidence);
            }
        }
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
        // Количество точек
        size_t n = pipe.profile.get_point_count();

        // Инициализация реологии
        auto& current = buffer.current();

        // Инициализация начального профиля плотности (не важно, ячейки или точки)
        for (double& density : current.density_std.value) {
            density = initial_conditions.density;
        }

        current.std_volumetric_flow = initial_conditions.volumetric_flow;

        //// Начальный гидравлический расчет
        calc_pressure_layer(initial_conditions);
    }
private:
    /// @brief Проводится расчёт шага движения партии
    /// @param dt Временной шаг моделирования
    /// @param boundaries Краевые условия
    void make_rheology_step(double dt, const boundaries_type& boundaries) {
        advance(); // Сдвигаем текущий и предыдущий слои

        // Преобразуем граничные условия в формат endogenous_values_t
        pde_solvers::endogenous_values_t endogenous_boundaries;
        endogenous_boundaries.density_std.value = boundaries.density;
        endogenous_boundaries.density_std.confidence = true;

        iso_nonbarotropic_pipe_solver_t solver(pipe, buffer);
        solver.transport_step(dt, boundaries.volumetric_flow, endogenous_boundaries);
    }

    /// @brief Рассчёт профиля давления методом Эйлера (задача PQ)
    /// @param boundaries Краевые условия
    void calc_pressure_layer(const boundaries_type& boundaries) {
        iso_nonbarotropic_pipe_solver_t solver(pipe, buffer);
        solver.hydro_solve_PQ(boundaries.volumetric_flow, boundaries.pressure_in);
    }
public:
    /// @brief Рассчёт шага моделирования, включающий в себя расчёт шага движения партии и гидравлический расчёт
    /// buffer.current после вызова содержит свежерасчитанный слой
    /// @param dt временной шаг моделирования
    /// @param boundaries Краевые условие
    void step(double dt, const boundaries_type& boundaries) {
        make_rheology_step(dt, boundaries);
        calc_pressure_layer(boundaries);
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
    /// @param pressure_initial Начальное значение давления для метода Ньютона
    void solve(const boundaries_type& initial_conditions, double volumetric_flow_initial = 0.2)
    {
        // Количество точек
        size_t n = pipe.profile.get_point_count();

        // Инициализация реологии
        auto& current = buffer.current();
        for (double& density : current.density_std.value) {
            density = initial_conditions.density;
        }

        //// Начальный гидравлический расчет
        calc_pressure_layer(initial_conditions, volumetric_flow_initial);
    }
private:
    /// @brief Проводится расчёт шага движения партии
    /// @param dt Временной шаг моделирования
    /// @param boundaries Краевые условия
    void make_rheology_step(double dt, const boundaries_type& boundaries) {
        // Сохраняем расход из текущего слоя до advance()
        double volumetric_flow = buffer.current().std_volumetric_flow;
        
        advance(); // Сдвигаем текущий и предыдущий слои

        // Преобразуем граничные условия в формат endogenous_values_t
        pde_solvers::endogenous_values_t endogenous_boundaries;
        endogenous_boundaries.density_std.value = boundaries.density;
        endogenous_boundaries.density_std.confidence = true;

        iso_nonbarotropic_pipe_solver_t solver(pipe, buffer);
        solver.transport_step(dt, volumetric_flow, endogenous_boundaries);
    }

    /// @brief Рассчёт профиля давления методом Ньютона над Эйлером (задача PP)
    /// @param boundaries Краевые условия
    /// @param volumetric_flow_initial Начальное значение расхода для метода Ньютона
    void calc_pressure_layer(const boundaries_type& boundaries, double volumetric_flow_initial) {
        iso_nonbarotropic_pipe_solver_t solver(pipe, buffer);
        solver.hydro_solve_PP(boundaries.pressure_in, boundaries.pressure_out);
    }
public:
    /// @brief Рассчёт шага моделирования, включающий в себя расчёт шага движения партии и гидравлический расчёт
    /// buffer.current после вызова содержит свежерасчитанный слой
    /// @param dt временной шаг моделирования
    /// @param boundaries Краевые условие
    void step(double dt, const boundaries_type& boundaries) {
        make_rheology_step(dt, boundaries);
        double volumetric_flow_initial = buffer.previous().std_volumetric_flow;
        calc_pressure_layer(boundaries, volumetric_flow_initial);
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
        throw std::runtime_error("not implemented");
        //auto& current = buffer.current();

        //// Проверяем наличие данных о плотности
        //if (current.density_std.value.empty()) {
        //    throw std::runtime_error("density profile is empty");
        //}

        //// Создаем граничные условия для задачи PP
        //iso_nonbarotropic_pipe_PP_task_boundaries_t boundaries;
        //boundaries.pressure_in = pressure_input;
        //boundaries.pressure_out = pressure_output;
        //boundaries.density = current.density_std.value[0];

        //// TODO: задавать начальное приближение расхода (брать из настроек солвера?)
        //double volumetric_flow_initial = 0.2;
        //int euler_direction = +1;
        //iso_nonbaro_impulse_equation_t pipeModel(pipe, current, volumetric_flow_initial, euler_direction);

        //// Создаем объект класса для расчета невязки при решении PP задачи методом Ньютона
        //solve_condensate_PP<iso_nonbaro_impulse_equation_t, iso_nonbarotropic_pipe_PP_task_boundaries_t, layer_type> solver_pp(
        //    pipeModel, boundaries, current);

        //fixed_solver_parameters_t<1, 0, golden_section_search> parameters;
        //parameters.residuals_norm = 0.1; // погрешность 0.1 Па
        //parameters.argument_increment_norm = 0;
        //parameters.residuals_norm_allow_early_exit = true;

        //// Создание структуры для записи результатов расчета
        //fixed_solver_result_t<1> result;
        //fixed_newton_raphson<1>::solve_dense(solver_pp, { volumetric_flow_initial }, parameters, &result);

        //// Обновляем расход в текущем слое
        //current.std_volumetric_flow = result.argument;

        //return result.argument;
    }

    /// @brief Решение гидравлической задачи QP
    virtual double hydro_solve_QP(double volumetric_flow, double pressure_output) override {
        auto& current = buffer.current();

        // Проверяем наличие данных о плотности
        if (current.density_std.value.empty()) {
            throw std::runtime_error("iso_nonbarotropic_pipe_solver_t::hydro_solve_QP: density profile is empty");
        }

        // Рассчитываем профиль давления методом Эйлера в обратном направлении (от выхода ко входу)
        std::vector<double>& p_profile = current.pressure;
        int euler_direction = -1;
        iso_nonbaro_improver_impulse_equation_t pipeModel(pipe, current, volumetric_flow, euler_direction);
        solve_euler<1>(pipeModel, euler_direction, pressure_output, &p_profile);

        // Обновляем расход в текущем слое
        current.std_volumetric_flow = volumetric_flow;

        // Возвращаем давление на входе (первый элемент профиля)
        return p_profile.front();
    }

    /// @brief Решение гидравлической задачи PQ
    /// @return Давление на выходе, Па
    virtual double hydro_solve_PQ(double volumetric_flow, double pressure_in) override {
        auto& current = buffer.current();

        // Проверяем наличие данных о плотности
        if (current.density_std.value.empty()) {
            throw std::runtime_error("iso_nonbaro_improver_pipe_solver_t::hydro_solve_PQ: density profile is empty");
        }
        if (current.improver_concentration.value.empty()) {
            throw std::runtime_error("iso_nonbaro_improver_pipe_solver_t::hydro_solve_PQ: concentration profile is empty");
        }

        // Рассчитываем профиль давления методом Эйлера
        std::vector<double>& p_profile = current.pressure;
        int euler_direction = +1; // Задаем направление для Эйлера
        iso_nonbaro_improver_impulse_equation_t pipeModel(pipe, current, volumetric_flow, euler_direction);
        solve_euler<1>(pipeModel, euler_direction, pressure_in, &p_profile);

        // Обновляем расход в текущем слое
        current.std_volumetric_flow = volumetric_flow;

        // Возвращаем давление на выходе (последний элемент профиля)
        return p_profile.back();
    }

    /// @brief Выполнение транспортного шага (расчет движения партий)
    virtual void transport_step(double dt, double volumetric_flow, const pde_solvers::endogenous_values_t& boundaries) override {
        buffer.current().std_volumetric_flow = volumetric_flow;

        if (std::isinf(dt)) {
            {
                auto value_buffer_wrapper = buffer.get_buffer_wrapper(
                    &iso_nonbaro_improver_pipe_endogenious_layer_t::get_density_wrapper);
                std::vector<double>& val = value_buffer_wrapper.current().vars;
                std::fill(val.begin(), val.end(), boundaries.density_std.value);

                auto confidence_buffer_wrapper = buffer.get_buffer_wrapper(
                    &iso_nonbaro_improver_pipe_endogenious_layer_t::get_density_confidence_wrapper);
                std::vector<double>& conf = confidence_buffer_wrapper.current().vars;
                std::fill(conf.begin(), conf.end(), boundaries.density_std.confidence);
            }
            {
                auto concentration_buffer_wrapper = buffer.get_buffer_wrapper(
                    &iso_nonbaro_improver_pipe_endogenious_layer_t::get_improver_concentration_wrapper);
                std::vector<double>& concentration_val = concentration_buffer_wrapper.current().vars;
                std::fill(concentration_val.begin(), concentration_val.end(), boundaries.improver.value);

                auto confidence_concentration_wrapper = buffer.get_buffer_wrapper(
                    &iso_nonbaro_improver_pipe_endogenious_layer_t::get_density_confidence_wrapper);
                std::vector<double>& concentration_conf = confidence_concentration_wrapper.current().vars;
                std::fill(concentration_conf.begin(), concentration_conf.end(), boundaries.improver.confidence);
            }

        }
        else {
            // считаем партии с помощью QUICKEST-ULTIMATE
            pde_solvers::pipe_advection_pde_t advection_model(pde_solvers::circle_area(pipe.wall.diameter),
                volumetric_flow, pipe.profile.coordinates);
            {
                auto buffer_wrapper = buffer.get_buffer_wrapper(
                    &iso_nonbaro_improver_pipe_endogenious_layer_t::get_density_wrapper);
                pde_solvers::quickest_ultimate_fv_solver solver(advection_model, buffer_wrapper);
                solver.step(dt, boundaries.density_std.value, boundaries.density_std.value);
            }
            {
                auto buffer_wrapper = buffer.get_buffer_wrapper(
                    &iso_nonbaro_improver_pipe_endogenious_layer_t::get_density_confidence_wrapper);
                pde_solvers::quickest_ultimate_fv_solver solver(advection_model, buffer_wrapper);
                solver.step(dt, boundaries.density_std.confidence, boundaries.density_std.confidence);
            }

            {
                auto buffer_wrapper = buffer.get_buffer_wrapper(
                    &iso_nonbaro_improver_pipe_endogenious_layer_t::get_improver_concentration_wrapper);
                pde_solvers::quickest_ultimate_fv_solver solver(advection_model, buffer_wrapper);
                solver.step(dt, boundaries.improver.value, boundaries.improver.value);
            }
            {
                auto buffer_wrapper = buffer.get_buffer_wrapper(
                    &iso_nonbaro_improver_pipe_endogenious_layer_t::get_improver_concentration_confidence_wrapper);
                pde_solvers::quickest_ultimate_fv_solver solver(advection_model, buffer_wrapper);
                solver.step(dt, boundaries.improver.confidence, boundaries.improver.confidence);
            }

        }
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

