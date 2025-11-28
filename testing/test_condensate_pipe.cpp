#define GTEST_BREAK_ON_FAILURE 1
#define GTEST_CATCH_EXCEPTIONS 0
#define GTEST_HAS_SEH 0
#define _VARIADIC_MAX 10 /* for gtest */
#include "gtest/gtest.h"

#include <fixed/fixed.h>
#include <pde_solvers/pde_solvers.h>

namespace pde_solvers {
;

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


struct condensate_pipe_properties_t : public pipe_properties_t {
    double kinematic_viscosity{1e-7};
};


/// @brief Уравнение сохранения импульса для трубы с учетом движения партий
class condensate_pipe_PQ_parties_t : public ode_t<1>
{
public:
    using ode_t<1>::equation_coeffs_type;
    using ode_t<1>::right_party_type;
    using ode_t<1>::var_type;
protected:
    const std::vector<double>& rho_profile;
    const condensate_pipe_properties_t& pipe;
    const double flow;
    const int solver_direction;
public:
    /// @brief Констуктор уравнения трубы
    /// @param pipe Ссылка на сущность трубы
    /// @param oil Ссылка на сущность нефти
    /// @param flow Объемный расход
    /// @param solver_direction Направление расчета по Эйлеру, должно обязательно совпадать с параметром солвера Эйлера
    condensate_pipe_PQ_parties_t(const condensate_pipe_properties_t& pipe,
        const std::vector<double>& rho_profile, 
        double flow,
        int solver_direction)
        : pipe(pipe)
        , rho_profile(rho_profile)
        , flow(flow)
        , solver_direction(solver_direction)
    {
    }

    /// @brief Возвращает известную уравнению сетку
    virtual const std::vector<double>& get_grid() const override {
        return pipe.profile.coordinates;
    }

    /// @brief Возвращает значение правой части ДУ
    /// @param grid_index Обсчитываемый индекс расчетной сетки
    /// @param point_vector Начальные условия
    /// @return Значение правой части ДУ в точке point_vector
    virtual right_party_type ode_right_party(
        size_t grid_index, const var_type& point_vector) const override
    {

        /// Обработка индекса в случае расчетов на границах трубы
        /// Чтобы не выйти за массив высот, будем считать dz/dx в соседней точке
        size_t rheo_index = grid_index;

        if (pipe.profile.get_point_count() == rho_profile.size())
        {
            // Случай расчета партий в точках (например для метода характеристик)
            if (solver_direction == +1)
                rheo_index += 1;
            else
                rheo_index -= 1;
        }
        else
        {
            // Случай расчета партий в ячейках (например для quickest ultimate) 
            rheo_index = solver_direction == +1
                ? grid_index
                : grid_index - 1;
        }
        double rho = rho_profile[rheo_index];
        double S_0 = pipe.wall.getArea();
        double v = flow / (S_0);
        double Re = v * pipe.wall.diameter / pipe.kinematic_viscosity;
        double lambda = pipe.resistance_function(Re);
        double tau_w = lambda / 8 * rho * v * abs(v);
        double height_derivative = pipe.profile.get_height_derivative(grid_index, solver_direction);
        double result = -4 * tau_w / pipe.wall.diameter - rho * M_G * height_derivative;
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


}

// TODO: Тест на расчет гидравлики
// TODO: Тест на расчет 

/// @brief Пример тестирования неизотермического расчета (Удалить !!!)
/// @param  
/// @param  
TEST(nonisothermal, example)
{
    pde_solvers::pipe_noniso_properties_t pipe;
    double x0 = 0;
    double xl = 1e5;
    double d = 1;
    size_t n = 100;

    pipe.profile = pde_solvers::pipe_profile_t::create(n, x0, xl, 0, 0, 10e6);
    pipe.wall.diameter = d;

    pde_solvers::oil_parameters_t oil;
    oil.viscosity.nominal_viscosity = 6e-7;

    pde_solvers::noniso_qsm_model_type model_type =
        pde_solvers::noniso_qsm_model_type::Isothemal;

    pde_solvers::qsm_noniso_TP_task_boundaries_t model(
        pipe,       // передаём как const ref
        oil,        // передаём как const ref
        model_type  // передаём значение
    );

    pde_solvers::nonisothermal_quasistatic_PQ_task_boundaries_t_p initial_conditions =
        pde_solvers::nonisothermal_quasistatic_PQ_task_boundaries_t_p::default_values();

    double dt = 10;
    model.solve(initial_conditions);
    model.step(dt, initial_conditions);

}

TEST(Lurie, qp)
{
    pde_solvers::pipe_properties_t model;
    double x0 = 0;
    double xl = 80000; // 80 км
    size_t n = 1;
    double z0 = 100;
    double zn = 50;
    model.profile = pde_solvers::pipe_profile_t::create(n, x0, xl, z0, zn, 10e6);

    model.wall.diameter = 0.7;
    model.wall.wallThickness = 10e-3; // 10 мм
    model.wall.equivalent_roughness = 0.015 * 1e-3;


    pde_solvers::isothermal_quasistatic_PQ_task_t Solver (model, pde_solvers::QuasistaticModelType::Stationary);
    pde_solvers::isothermal_quasistatic_PQ_task_boundaries_t initial_conditions;
    initial_conditions.pressure_in = 5.65 * 1e6;
    initial_conditions.density = 870;
    initial_conditions.volumetric_flow = 3500.0 / 3600.0;
    initial_conditions.viscosity = 15 * 1e-6; //15 сСт
    Solver.solve(initial_conditions);
    double relation = ((600000 - Solver.get_buffer().get_layers()[0].pressure[n]) / 600000.0);
    EXPECT_TRUE(relation < 0.011);


}

TEST(Lurie, pp)
{
    pde_solvers::pipe_properties_t model;
    double x0 = 0;
    double xl = 80000; // 80 км
    size_t n = 1;
    double z0 = 100;
    double zn = 50;
    model.profile = pde_solvers::pipe_profile_t::create(n, x0, xl, z0, zn, 10e6);

    model.wall.diameter = 0.7;
    model.wall.wallThickness = 10e-3; // 10 мм
    model.wall.equivalent_roughness = 0.015 * 1e-3;


    pde_solvers::isothermal_quasistatic_PP_task_t Solver(model, pde_solvers::QuasistaticModelType::Stationary);
    pde_solvers::isothermal_quasistatic_PP_task_boundaries_t initial_conditions;
    initial_conditions.pressure_in = 5.65 * 1e6;
    initial_conditions.density = 870;
    initial_conditions.pressure_out = 0.6 * 1e6;
    initial_conditions.viscosity = 15 * 1e-6; //15 сСт
    Solver.solve(initial_conditions);
    double relation = ((0.9722 - Solver.get_buffer().get_layers()[0].volume) / 0.9722);
    EXPECT_TRUE(relation < 0.01);

    Solver.step(10, initial_conditions);


}

TEST(Lurie, Newton)
{
    class sample_system : public fixed_system_t<2>
    {
        using fixed_system_t<2>::var_type;

    public:
        // Задание функции невязок
        var_type residuals(const var_type& x) {
            return
            {
                pow(x[0] - 2, 3) - 8.0,
                pow(x[1] - 1, 3) - 27.0
            };
        }
    };

    double eps = 1e-6;
    const size_t razmernost = 2;
    // Создание экземпляра класса, который и будет решаемой системой
    sample_system test;
    // Задание настроек решателя по умолчанию
    fixed_solver_parameters_t<razmernost, 0> parameters; // <2, 0> <размерность, >
    parameters.argument_increment_norm = eps;
    // Создание структуры для записи результатов расчета
    fixed_solver_result_t<razmernost> result;
    // Решение системы нелинейныйх уравнений <2> с помощью решателя Ньютона - Рафсона
    // { 0, 0 } - Начальное приближение
    fixed_newton_raphson<razmernost>::solve_dense(test, { 0, 0 }, parameters, &result);
	// Проверка результатов
    EXPECT_NEAR(result.argument[0], 4.0, 1e-6);
	EXPECT_NEAR(result.argument[1], 4.0, 1e-6);
}

//TEST(CondensatePipe, Develop) {
//
//}