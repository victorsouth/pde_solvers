#pragma once

namespace pde_solvers {
;

/// @brief Профиль параметров для конденсатопровода (без температуры и ПТП)
struct condensate_pipe_layer {
    /// @brief Номинальный объемный расход
    double std_volumetric_flow{ std::numeric_limits<double>::quiet_NaN() };
    /// @brief Профиль давления
    std::vector<double> pressure;
    /// @brief Профиль плотности с достоверностью
    confident_layer_t density;
    /// @brief Профиль вспомогательных расчетов для метода конечных объемов (и для вязкости, и для плотности)
    quickest_ultimate_fv_solver_traits<1>::specific_layer specific;
    /// @brief Инициализация профилей
    /// @param point_count Количество точек
    condensate_pipe_layer(size_t point_count)
        : density(point_count, 850.0)
        , specific(point_count)
        , pressure(point_count)
    {
    }

    /// @brief Подготовка плотности для расчета методом конечных объемов 
    /// @param layer Слой
    /// @return Обертка над составным слоем
    static quickest_ultimate_fv_wrapper<1> get_density_wrapper(condensate_pipe_layer& layer)
    {
        return quickest_ultimate_fv_wrapper<1>(layer.density.value, layer.specific);
    }

};

/// @brief Структура, содержащая в себе краевые условия задачи PQ
struct condensate_pipe_PQ_task_boundaries_t {
    /// @brief Изначальный объемный расход
    double volumetric_flow;
    /// @brief Изначальное давление на входе
    double pressure_in;
    /// @brief Изначальная плотность на входе
    double density;
    /// @brief Создание структуры со значениями по умолчанию
    static condensate_pipe_PQ_task_boundaries_t default_values() {
        condensate_pipe_PQ_task_boundaries_t result;
        result.volumetric_flow = 0.2;
        result.pressure_in = 6e6;
        result.density = 850;
        return result;
    }
};


class condensate_pipe_PQ_task_t {
public:
    /// @brief Тип слоя
    using layer_type = condensate_pipe_layer;
    /// @brief Тип буфера
    using buffer_type = ring_buffer_t<layer_type>;
    /// @brief Тип граничных условий
    using boundaries_type = condensate_pipe_PQ_task_boundaries_t;
private:
    // Модель трубы
    condensate_pipe_properties_t pipe;
    // Создаётся буфер, тип слоя которого определяется в зависимости от типа солвера
    buffer_type buffer;
public:
    /// @brief Конструктор
    /// @param pipe Модель трубопровода
    condensate_pipe_PQ_task_t(const condensate_pipe_properties_t& pipe)
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
        for (double& density : current.density.value) {
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
        size_t n = pipe.profile.get_point_count();
        std::vector<double>Q_profile(n, boundaries.volumetric_flow); // задаем по трубе новый расход из временного ряда

        advance(); // Сдвигаем текущий и предыдущий слои

        buffer.current().std_volumetric_flow = boundaries.volumetric_flow; 

        // считаем партии с помощью QUICKEST-ULTIMATE
        PipeQAdvection advection_model(pipe, Q_profile);

        // Шаг по плотности
        auto density_wrapper = buffer.get_buffer_wrapper(
            &condensate_pipe_layer::get_density_wrapper);
        quickest_ultimate_fv_solver solver_rho(advection_model, density_wrapper);
        solver_rho.step(dt, boundaries.density, boundaries.density);
    }

    /// @brief Рассчёт профиля давления методом Эйлера (задача PQ)
    /// @param boundaries Краевые условия
    void calc_pressure_layer(const boundaries_type& boundaries) {

        auto& current = buffer.current();

        std::vector<double>& p_profile = current.pressure;
        int euler_direction = +1; // Задаем направление для Эйлера

        condensate_pipe_PQ_parties_t pipeModel(pipe, current.density.value, boundaries.volumetric_flow, euler_direction);
        solve_euler<1>(pipeModel, euler_direction, boundaries.pressure_in, &p_profile);
    }
public:
    /// @brief Рассчёт шага моделирования, включающий в себя расчёт шага движения партии и гидравлический расчёт
    /// Функция делат сдвиг буфера (advance) так, что buffer.current после вызова содержит свежерасчитанный слой
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
    condensate_pipe_layer& get_current_layer() {
        return buffer.current();
    }
};

/// @brief Структура, содержащая в себе краевые условия задачи PP
struct condensate_pipe_PP_task_boundaries_t {
    /// @brief Изначальное давление на входе
    double pressure_in;
    /// @brief Изначальное давление на выходе
    double pressure_out;
    /// @brief Изначальная плотность на входе
    double density;
    /// @brief Создание структуры со значениями по умолчанию
    static condensate_pipe_PP_task_boundaries_t default_values() {
        condensate_pipe_PP_task_boundaries_t result;
        result.pressure_out = 0.6e6;
        result.pressure_in = 6e6;
        result.density = 850;
        return result;
    }
};


/// @brief класс для нахождения расхода Q для задачи PP с помощью метода Ньютона
/// @tparam boundaries_type класс граничных условий
/// @tparam layer_type класс уровней в buffer
template <typename BoundariesType, typename LayerType>
class solve_condensate_PP : public fixed_system_t<1> {
    using fixed_system_t<1>::var_type;
private:
    /// @brief слой расчета
    LayerType& current_layer;
    /// @brief ГУ
    const BoundariesType& bound;
    /// @brief свойства трубы
    const condensate_pipe_properties_t& pipe;

public:
    solve_condensate_PP(const condensate_pipe_properties_t& pipe, const BoundariesType& bound, LayerType& current_layer)
        : pipe(pipe)
        , bound(bound)
        , current_layer(current_layer)
    {
    }

    /// @brief функция невязки для решения методом Ньютона
    /// @param x - неизвестное (для задачи PP является расходом)
    /// @return 
    virtual double residuals(const double& x) {
        auto& current = current_layer;

        std::vector<double>& p_profile = current.pressure;
        int euler_direction = +1; // Задаем направление для Эйлера
        condensate_pipe_PQ_parties_t pipeModel(pipe, current.density.value, x, euler_direction);
        solve_euler<1>(pipeModel, euler_direction, bound.pressure_in, &p_profile);

        return p_profile.back() - bound.pressure_out;
    }

    /// @brief переопределяем целевую функцию, чтобы был модуль невязок
    /// @param r 
    /// @return 
    virtual double objective_function(const var_type& r) const override {
        return std::abs(r);
    }

};


class condensate_pipe_PP_task_t {
public:
    /// @brief Тип слоя
    using layer_type = condensate_pipe_layer;
    /// @brief Тип буфера
    using buffer_type = ring_buffer_t<layer_type>;
    /// @brief Тип граничных условий
    using boundaries_type = condensate_pipe_PP_task_boundaries_t;
private:
    // Модель трубы
    condensate_pipe_properties_t pipe;
    // Создаётся буфер, тип слоя которого определяется в зависимости от типа солвера
    buffer_type buffer;


public:
    /// @brief Конструктор
    /// @param pipe Модель трубопровода
    condensate_pipe_PP_task_t(const condensate_pipe_properties_t& pipe)
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

        // Инициализация начального профиля плотности (не важно, ячейки или точки)
        for (double& density : current.density.value) {
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
        size_t n = pipe.profile.get_point_count();
        std::vector<double>Q_profile(n, buffer.current().std_volumetric_flow); // задаем по трубе новый расход из временного ряда

        advance(); // Сдвигаем текущий и предыдущий слои

        // считаем партии с помощью QUICKEST-ULTIMATE
        PipeQAdvection advection_model(pipe, Q_profile);

        // Шаг по плотности
        auto density_wrapper = buffer.get_buffer_wrapper(
            &condensate_pipe_layer::get_density_wrapper);
        quickest_ultimate_fv_solver solver_rho(advection_model, density_wrapper);
        solver_rho.step(dt, boundaries.density, boundaries.density);
    }

    /// @brief Рассчёт профиля давления методом Ньютона над Эйлером (задача PP)
    /// @param boundaries Краевые условия
    /// @param pressure_initial Начальное значение расхода для метода Ньютона
    void calc_pressure_layer(const boundaries_type& boundaries, double volumetric_flow_initial) {

        auto& current = buffer.current();

        // создаем объект класса для расчета невязки при решении PP задачи методом Ньютона
        solve_condensate_PP<boundaries_type, layer_type> test = solve_condensate_PP(pipe, boundaries, current);
        fixed_solver_parameters_t<1, 0, golden_section_search> parameters;
        parameters.residuals_norm = 0.1; // погрешность 0.1 Па
        parameters.argument_increment_norm = 0;
        parameters.residuals_norm_allow_early_exit = true;
        // Создание структуры для записи результатов расчета
        fixed_solver_result_t<1> result;
        fixed_newton_raphson<1>::solve_dense(test, { volumetric_flow_initial }, parameters, &result);
        current.std_volumetric_flow = result.argument;
    }
public:
    /// @brief Рассчёт шага моделирования, включающий в себя расчёт шага движения партии и гидравлический расчёт
    /// Функция делат сдвиг буфера (advance) так, что buffer.current после вызова содержит свежерасчитанный слой
    /// @param dt временной шаг моделирования
    /// @param boundaries Краевые условие
    void step(double dt, const boundaries_type& boundaries) {
        make_rheology_step(dt, boundaries);
        // берем давление на выходе из предыдущего слоя как начальное для нового расчета
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
    condensate_pipe_layer& get_current_layer() {
        return buffer.current();
    }
};

/// @brief Солвер квазистационарного гидравлического расчета для конденсатопровода
class iso_nonbarotropic_pipe_solver_t : public pipe_solver_interface_t {
public:
    /// @brief Тип слоя
    using layer_type = condensate_pipe_layer;
    /// @brief Тип буфера
    using buffer_type = ring_buffer_t<layer_type>;
    /// @brief Тип параметров трубы
    using pipe_parameters_type = condensate_pipe_properties_t;

private:
    /// @brief Ссылка на свойства конденсатопровода
    const condensate_pipe_properties_t& pipe;
    /// @brief Ссылка на буфер слоев
    buffer_type& buffer;

public:
    /// @brief Конструктор
    /// @param pipe Ссылка на свойства конденсатопровода
    /// @param buffer Ссылка на буфер слоев
    iso_nonbarotropic_pipe_solver_t(
        const condensate_pipe_properties_t& pipe,
        buffer_type& buffer)
        : pipe(pipe)
        , buffer(buffer)
    {
    }

    /// @brief Решение гидравлической задачи PP (заданы давления на входе и выходе, найти расход)
    virtual double hydro_solve_PP(double pressure_input, double pressure_output) override {
        auto& current = buffer.current();

        // Проверяем наличие данных о плотности
        if (current.density.value.empty()) {
            throw std::runtime_error("density profile is empty");
        }

        // Создаем граничные условия для задачи PP
        condensate_pipe_PP_task_boundaries_t boundaries;
        boundaries.pressure_in = pressure_input;
        boundaries.pressure_out = pressure_output;
        boundaries.density = current.density.value[0];

        // Используем начальное значение расхода из текущего слоя, если оно есть
        if (std::isnan(current.std_volumetric_flow)) {
            throw std::runtime_error("std_volumetric_flow is not initialized");
        }
        double volumetric_flow_initial = current.std_volumetric_flow;

        // Создаем объект класса для расчета невязки при решении PP задачи методом Ньютона
        solve_condensate_PP<condensate_pipe_PP_task_boundaries_t, layer_type> solver_pp(
            pipe, boundaries, current);
            
        fixed_solver_parameters_t<1, 0, golden_section_search> parameters;
        parameters.residuals_norm = 0.1; // погрешность 0.1 Па
        parameters.argument_increment_norm = 0;
        parameters.residuals_norm_allow_early_exit = true;
            
        // Создание структуры для записи результатов расчета
        fixed_solver_result_t<1> result;
        fixed_newton_raphson<1>::solve_dense(solver_pp, { volumetric_flow_initial }, parameters, &result);
            
        // Обновляем расход в текущем слое
        current.std_volumetric_flow = result.argument;
            
        return result.argument;
    }

    /// @brief Решение гидравлической задачи QP (заданы расход и давление на выходе, найти давление на входе)
    virtual double hydro_solve_QP(double volumetric_flow, double pressure_output) override {
        auto& current = buffer.current();

        // Проверяем наличие данных о плотности
        if (current.density.value.empty()) {
            throw std::runtime_error("iso_nonbarotropic_pipe_solver_t::hydro_solve_QP: density profile is empty");
        }

        // Рассчитываем профиль давления методом Эйлера в обратном направлении (от выхода ко входу)
        std::vector<double>& p_profile = current.pressure;
        int euler_direction = -1; // Задаем обратное направление для Эйлера (от выхода ко входу)
        condensate_pipe_PQ_parties_t pipeModel(pipe, current.density.value, volumetric_flow, euler_direction);
        solve_euler<1>(pipeModel, euler_direction, pressure_output, &p_profile);

        // Обновляем расход в текущем слое
        current.std_volumetric_flow = volumetric_flow;

        // Возвращаем давление на входе (первый элемент профиля)
        return p_profile.front();
    }

    /// @brief Решение гидравлической задачи PQ (заданы расход и давление на входе, найти давление на выходе)
    /// @return Давление на выходе, Па
    virtual double hydro_solve_PQ(double volumetric_flow, double pressure_in) override {
        auto& current = buffer.current();

        // Проверяем наличие данных о плотности
        if (current.density.value.empty()) {
            throw std::runtime_error("iso_nonbarotropic_pipe_solver_t::hydro_solve_PQ: density profile is empty");
        }

        // Рассчитываем профиль давления методом Эйлера
        std::vector<double>& p_profile = current.pressure;
        int euler_direction = +1; // Задаем направление для Эйлера
        condensate_pipe_PQ_parties_t pipeModel(pipe, current.density.value, volumetric_flow, euler_direction);
        solve_euler<1>(pipeModel, euler_direction, pressure_in, &p_profile);

        // Обновляем расход в текущем слое
        current.std_volumetric_flow = volumetric_flow;

        // Возвращаем давление на выходе (последний элемент профиля)
        return p_profile.back();
    }

    /// @brief Вычисление якобиана для задачи PP
    /// Использует формулу переворота производной: dQ/dP = 1/(dP/dQ)
    /// @param pressure_input Давление на входе, Па
    /// @param pressure_output Давление на выходе, Па
    /// @return Массив из двух элементов: [dQ/dP_in, dQ/dP_out]
    virtual std::array<double, 2> hydro_solve_PP_jacobian(double pressure_input, double pressure_output) override {
        // Вычисляем базовое решение - расход при заданных давлениях
        double Q_base = hydro_solve_PP(pressure_input, pressure_output);
            
        // Малое приращение для численного дифференцирования (0.1% от расхода)
        const double eps = std::max(1e-6, std::abs(Q_base) * 1e-3);
            
        // Вычисляем производную перепада давления по расходу dP/dQ
        // Вычисляем давление на выходе при базовом и увеличенном расходе
        double P_out_base = hydro_solve_PQ(Q_base, pressure_input);
        double P_out_plus = hydro_solve_PQ(Q_base + eps, pressure_input);
            
        // Производная давления на выходе по расходу: dP_out/dQ
        double dP_out_dQ = (P_out_plus - P_out_base) / eps;
            
        // Производная перепада давления по расходу: dP/dQ = -dP_out/dQ
        double dP_dQ = -dP_out_dQ;
            
        // Используем формулу переворота производной:
        double dQ_dP_in = 1.0 / dP_dQ;
        double dQ_dP_out = -1.0 / dP_dQ;
            
        return { dQ_dP_in, dQ_dP_out };
    }

    /// @brief Выполнение транспортного шага (расчет движения партий)
    virtual void transport_step(double dt, double volumetric_flow, const pde_solvers::endogenous_values_t& boundaries) override {
        size_t n = pipe.profile.get_point_count();
        std::vector<double> Q_profile(n, volumetric_flow); // задаем по трубе расход

        buffer.advance(+1); // Сдвигаем текущий и предыдущий слои

        auto& current = buffer.current();
        current.std_volumetric_flow = volumetric_flow;

        // Считаем партии с помощью QUICKEST-ULTIMATE
        PipeQAdvection advection_model(pipe, Q_profile);

        // Шаг по плотности
        auto density_wrapper = buffer.get_buffer_wrapper(
            &condensate_pipe_layer::get_density_wrapper);
        quickest_ultimate_fv_solver solver_rho(advection_model, density_wrapper);
            
        // Используем плотность из граничных условий, если она задана и достоверна
        double density_in;
        if (boundaries.density.confidence) {
            density_in = boundaries.density.value;
        } else if (!current.density.value.empty()) {
            density_in = current.density.value[0];
        } else {
            throw std::runtime_error("density is not available in boundaries and density profile is empty");
        }
        double density_out = density_in; // Для простоты используем то же значение на выходе
            
        solver_rho.step(dt, density_in, density_out);
    }

    /// @brief Геттер для текущего слоя
    condensate_pipe_layer& get_current_layer() {
        return buffer.current();
    }
};

}

