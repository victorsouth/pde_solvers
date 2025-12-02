#pragma once

namespace pde_solvers {

/// @brief Профиль параметров для конденсатопровода (без температуры и ПТП)
struct condensate_pipe_layer {
    /// @brief Номинальный объемный расход
    double std_volumetric_flow{ std::numeric_limits<double>::quiet_NaN() };
    /// @brief Профиль давления
    std::vector<double> pressure;
    /// @brief Профиль плотности
    std::vector<double> density;
    /// @brief Профиль вспомогательных расчетов для метода конечных объемов (и для вязкости, и для плотности)
    quickest_ultimate_fv_solver_traits<1>::specific_layer specific;
    /// @brief Инициализация профилей
    /// @param point_count Количество точек
    condensate_pipe_layer(size_t point_count)
        : density(point_count - 1)
        , specific(point_count)
        , pressure(point_count)
    {
    }

    /// @brief Подготовка плотности для расчета методом конечных объемов 
    /// @param layer Слой
    /// @return Обертка над составным слоем
    static quickest_ultimate_fv_wrapper<1> get_density_wrapper(condensate_pipe_layer& layer)
    {
        return quickest_ultimate_fv_wrapper<1>(layer.density, layer.specific);
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
        for (double& density : current.density) {
            density = initial_conditions.density;
        }

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

        condensate_pipe_PQ_parties_t pipeModel(pipe, current.density, boundaries.volumetric_flow, euler_direction);
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
template <typename boundaries_type, typename layer_type>
class Solve_condesate_PP : public fixed_system_t<1> {
    using fixed_system_t<1>::var_type;
private:
    /// @brief слой расчета
    layer_type* current_layer;
    /// @brief ГУ
    boundaries_type bound;
    /// @brief свойства трубы
    const condensate_pipe_properties_t pipe;

public:
    Solve_condesate_PP(const condensate_pipe_properties_t& pipe)
        : pipe(pipe)
        , current_layer(nullptr)
    {
    }

    /// @brief функция невязки для решения методом Ньютона
    /// @param x - неизвестное (для задачи PP является расходом)
    /// @return 
    var_type residuals(const var_type& x) {
        auto& current = *current_layer;

        std::vector<double>& p_profile = current.pressure;
        int euler_direction = +1; // Задаем направление для Эйлера
        condensate_pipe_PQ_parties_t pipeModel(pipe, current.density, x, euler_direction);
        solve_euler<1>(pipeModel, euler_direction, bound.pressure_in, &p_profile);

        return
        {
            p_profile.back() - bound.pressure_out
        };
    }

    /// @brief задание новых ГУ и нового слоя расчета 
    /// @param new_bound новые ГУ
    /// @param new_current_layer новый слой расчета 
    void set_parametrs(const boundaries_type& new_bound, layer_type& new_current_layer) {
        bound = new_bound;
        current_layer = &new_current_layer;
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
    /// @brief объект для решения PP задачи методом Ньютона
    Solve_condesate_PP<boundaries_type, layer_type> test;


public:
    /// @brief Конструктор
    /// @param pipe Модель трубопровода
    condensate_pipe_PP_task_t(const condensate_pipe_properties_t& pipe)
        : pipe(pipe)
        , buffer(2, pipe.profile.get_point_count())
        , test(pipe)
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
        for (double& density : current.density) {
            density = initial_conditions.density;
        }

        //// Начальный гидравлический расчет
        calc_pressure_layer(initial_conditions);
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
    void calc_pressure_layer(const boundaries_type& boundaries) {

        auto& current = buffer.current();

        test.set_parametrs(boundaries, current);
        fixed_solver_parameters_t<1, 0> parameters;
        parameters.argument_increment_norm = 1e-6;
        // Создание структуры для записи результатов расчета
        fixed_solver_result_t<1> result;
        fixed_newton_raphson<1>::solve_dense(test, { 0.2 }, parameters, &result);
        current.std_volumetric_flow = result.argument;
        // !!!!!!!! Метод Ньютона считает с заданной точностью, поэтому давление на выходе немного отличается !!!!!!!
        current.pressure.back() = boundaries.pressure_out;
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

}

