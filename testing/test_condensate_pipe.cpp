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


/// @brief Уравнение сохранения импульса для трубы с учетом движения партий
class condensate_pipe_PP_parties_t : public ode_t<1>
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
    condensate_pipe_PP_parties_t(const condensate_pipe_properties_t& pipe,
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
        condensate_pipe_PP_parties_t pipeModel(pipe, current.density, x, euler_direction);
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


class CondensateQPTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Создаем простую равномерную трубу для тестов
        pde_solvers::pipe_profile_t profile;
        const size_t points_count = 10;
        const double length = 1000.0; // 1 км

        // Создаем равномерную сетку
        for (size_t i = 0; i < points_count; ++i) {
            double x = i * length / (points_count - 1);
            double h = 0.0; // Нулевой уклон для упрощения
            profile.coordinates.push_back(x);
            profile.heights.push_back(h);
        }

        pipe.profile = profile;
        pipe.wall.diameter = 0.5; // 500 мм
        pipe.wall.wallThickness = 10e-3; // 0.1 мм
        pipe.kinematic_viscosity = 1e-6; // 1 сСт
    }

    pde_solvers::condensate_pipe_properties_t pipe;
};


/// @brief 1. Тест на начальный стационарный расчет
TEST_F(CondensateQPTest, InitialSteadyStateSolution) {

    //Arrange
    pde_solvers::condensate_pipe_PQ_task_t task(pipe);

    auto initial_conditions = pde_solvers::condensate_pipe_PQ_task_boundaries_t::default_values();
    initial_conditions.volumetric_flow = 0.3; // м³/с
    initial_conditions.pressure_in = 5e6; // 5 МПа
    initial_conditions.density = 800.0; // кг/м³

    // Act
    task.solve(initial_conditions);
    auto& layer = task.get_current_layer();

    //Assert
    
    // Проверяем инициализацию плотности
    for (const auto& density : layer.density) {
        EXPECT_NEAR(density, 800.0, 1e-6);
    }

    // Проверяем расчет давления
    ASSERT_FALSE(layer.pressure.empty());

    // Давление на входе должно быть равно заданному
    EXPECT_NEAR(layer.pressure.front(), 5e6, 1e-6);

    // Давление должно уменьшаться вдоль трубы из-за гидравлических потерь
    for (size_t i = 1; i < layer.pressure.size(); ++i) {
        EXPECT_LT(layer.pressure[i], layer.pressure[i - 1]);
    }

    // Проверяем, что давление положительное
    for (const auto& pressure : layer.pressure) {
        EXPECT_GT(pressure, 0.0);
    }
}

/// @brief 2. Тест на зависимость потерь давления от расхода
TEST_F(CondensateQPTest, PressureLossVsFlowRate) {
    //Arrange
    pde_solvers::condensate_pipe_PQ_task_t task(pipe);

    std::vector<double> flows = { 0.1, 0.3, 0.5, 0.7 }; // м³/с
    std::vector<double> pressure_drops;

    auto initial_conditions = pde_solvers::condensate_pipe_PQ_task_boundaries_t::default_values();
    initial_conditions.pressure_in = 5e6;
    initial_conditions.density = 850.0;
    
    //Act
    for (double flow : flows) {
        initial_conditions.volumetric_flow = flow;
        task.solve(initial_conditions);

        auto& layer = task.get_current_layer();
        double pressure_drop = layer.pressure.front() - layer.pressure.back();
        pressure_drops.push_back(pressure_drop);

        // Для отладки
        std::cout << "Flow: " << flow << " m³/s, Pressure drop: "
            << pressure_drop << " Pa" << std::endl;
    }

    //Assert
    // Потери давления должны увеличиваться с ростом расхода
    for (size_t i = 1; i < pressure_drops.size(); ++i) {
        EXPECT_GT(pressure_drops[i], pressure_drops[i - 1]);
    }

}

/// @brief 3(Аналогичен тесту 2, можно удалить). Тест на влияние плотности на потери давления
TEST_F(CondensateQPTest, DensityEffectOnPressureLoss) {
    //Arrange
    pde_solvers::condensate_pipe_PQ_task_t task(pipe);

    std::vector<double> densities = { 700.0, 850.0, 1000.0 }; // кг/м³
    std::vector<double> pressure_drops;

    auto initial_conditions = pde_solvers::condensate_pipe_PQ_task_boundaries_t::default_values();
    initial_conditions.volumetric_flow = 0.3;
    initial_conditions.pressure_in = 5e6;

    //Act
    for (double density : densities) {
        initial_conditions.density = density;
        task.solve(initial_conditions);

        auto& layer = task.get_current_layer();
        double pressure_drop = layer.pressure.front() - layer.pressure.back();
        pressure_drops.push_back(pressure_drop);

        std::cout << "Density: " << density << " kg/m³, Pressure drop: "
            << pressure_drop << " Pa" << std::endl;
    }
    //Assert
    // Потери давления должны увеличиваться с ростом плотности
    for (size_t i = 1; i < pressure_drops.size(); ++i) {
        EXPECT_GT(pressure_drops[i], pressure_drops[i - 1]);
    }

    // Проверяем линейную зависимость потерь от плотности
    for (size_t i = 0; i < densities.size(); ++i) {
        double expected_ratio = densities[i] / densities[0];
        double actual_ratio = pressure_drops[i] / pressure_drops[0];

        EXPECT_NEAR(actual_ratio, expected_ratio, expected_ratio * 0.1);
    }
}

/// @brief 4(можно удалить). Тест на влияние уклона трубы
TEST_F(CondensateQPTest, PipeSlopeEffect) {
    //Arrange
    // Создаем трубу с уклоном
    pde_solvers::pipe_profile_t sloped_profile;
    const size_t points_count = 10;
    const double length = 1000.0;

    for (size_t i = 0; i < points_count; ++i) {
        double x = i * length / (points_count - 1);
        double h = 10.0 * x / length; // Уклон 1%
        sloped_profile.coordinates.push_back(x);
        sloped_profile.heights.push_back(h);
    }

    pipe.profile = sloped_profile;
    pde_solvers::condensate_pipe_PQ_task_t task(pipe);

    // Для сравнения создаем горизонтальную трубу
    pde_solvers::pipe_profile_t flat_profile;
    for (size_t i = 0; i < points_count; ++i) {
        double x = i * length / (points_count - 1);
        flat_profile.coordinates.push_back(x);
        flat_profile.heights.push_back(0.0);
    }

    pde_solvers::condensate_pipe_properties_t flat_pipe = pipe;
    flat_pipe.profile = flat_profile;
    pde_solvers::condensate_pipe_PQ_task_t flat_task(flat_pipe);

    auto initial_conditions = pde_solvers::condensate_pipe_PQ_task_boundaries_t::default_values();
    initial_conditions.volumetric_flow = 0.3;
    initial_conditions.pressure_in = 5e6;
    initial_conditions.density = 850.0;

    //Act
    task.solve(initial_conditions);

    auto& layer = task.get_current_layer();

    // При движении вверх по уклону перепад давления должен быть больше
    // чем при движении по горизонтали (дополнительные гравитационные потери)

    flat_task.solve(initial_conditions);

    double sloped_drop = layer.pressure.front() - layer.pressure.back();
    double flat_drop = flat_task.get_current_layer().pressure.front()
        - flat_task.get_current_layer().pressure.back();

    //Assert
    // Для трубы с уклоном вверх потери должны быть больше
    EXPECT_GT(sloped_drop, flat_drop);
}

/// @brief 5(можно удалить).Тест на корректность метода Эйлера
TEST_F(CondensateQPTest, EulerMethodConsistency) {
    //Arrange
    pde_solvers::condensate_pipe_PQ_task_t task(pipe);

    auto initial_conditions = pde_solvers::condensate_pipe_PQ_task_boundaries_t::default_values();

    //Act
    task.solve(initial_conditions);
    auto& initial_layer = task.get_current_layer();

    // Запоминаем начальное давление
    std::vector<double> initial_pressure = initial_layer.pressure;

    // Делаем шаг по времени с теми же условиями
    // Для стационарных условий профиль давления не должен меняться
    double dt = 10.0; // 10 секунд
    task.step(dt, initial_conditions);

    auto& new_layer = task.get_current_layer();

    //Assert
    // Проверяем, что давление изменилось минимально (только из-за численных погрешностей)
    double max_diff = 0.0;
    for (size_t i = 0; i < initial_pressure.size(); ++i) {
        double diff = std::abs(initial_pressure[i] - new_layer.pressure[i]);
        max_diff = std::max(max_diff, diff);
        EXPECT_LT(diff, 1e-3); // Допустимая погрешность 1 Па
    }

}

/// @brief 6. Тест на численную устойчивость
TEST_F(CondensateQPTest, NumericalStability) {
    //Arrange
    pde_solvers::condensate_pipe_PQ_task_t task(pipe);

    auto initial_conditions = pde_solvers::condensate_pipe_PQ_task_boundaries_t::default_values();
    task.solve(initial_conditions);

    // Выполняем много шагов и проверяем стабильность
    const int steps = 100;
    double dt = 1.0; // 1 секунда

    std::vector<double> inlet_pressures;
    std::vector<double> outlet_pressures;

    for (int i = 0; i < steps; ++i) {
        //Act
        task.step(dt, initial_conditions);
        auto& layer = task.get_current_layer();

        inlet_pressures.push_back(layer.pressure.front());
        outlet_pressures.push_back(layer.pressure.back());

        //Assert
        // Проверяем, что давление остается в физических пределах
        EXPECT_GT(layer.pressure.front(), 0.0);
        EXPECT_GT(layer.pressure.back(), 0.0);
    }

    //Assert
    // При постоянных граничных условиях давления должны стабилизироваться
    double inlet_variation = 0.0;
    double outlet_variation = 0.0;

    for (int i = 1; i < steps; ++i) {
        inlet_variation += std::abs(inlet_pressures[i] - inlet_pressures[i - 1]);
        outlet_variation += std::abs(outlet_pressures[i] - outlet_pressures[i - 1]);
    }

    inlet_variation /= (steps - 1);
    outlet_variation /= (steps - 1);

    // Вариации должны быть малы
    EXPECT_LT(inlet_variation, 1.0); // Менее 1 Па в среднем
    EXPECT_LT(outlet_variation, 1.0);
}

/// @brief 7.Тест на крайние случаи
TEST_F(CondensateQPTest, ZeroVolume) {
    // Очень маленький расход
    {
        //Arrange
        pde_solvers::condensate_pipe_PQ_task_t task(pipe);
        auto conditions = pde_solvers::condensate_pipe_PQ_task_boundaries_t::default_values();
        conditions.volumetric_flow = 0; // нулевой расход

        EXPECT_NO_THROW(task.solve(conditions));

        auto& layer = task.get_current_layer();
        // При нулевом расходе перепад давления должен быть минимальным
        double pressure_drop = layer.pressure.front() - layer.pressure.back();
        EXPECT_EQ(pressure_drop, 0); // 
    }
}


class CondensatePPTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Arrange: Создаем простую равномерную трубу для тестов
        pde_solvers::pipe_profile_t profile;
        const size_t points_count = 10;
        const double length = 1000.0; // 1 км

        // Создаем равномерную сетку
        for (size_t i = 0; i < points_count; ++i) {
            double x = i * length / (points_count - 1);
            double h = 0.0; // Нулевой уклон для упрощения
            profile.coordinates.push_back(x);
            profile.heights.push_back(h);
        }

        pipe.profile = profile;
        pipe.wall.diameter = 0.5; // 500 мм
        pipe.wall.wallThickness = 0.0001; // 0.1 мм
        pipe.kinematic_viscosity = 1e-6; // 1 сСт
    }

    pde_solvers::condensate_pipe_properties_t pipe;
};



/// @brief 8. Тест на начальный стационарный расчет PP задачи
TEST_F(CondensatePPTest, InitialSteadyStateSolution) {
    // Arrange
    pde_solvers::condensate_pipe_PP_task_t task(pipe);
    auto initial_conditions = pde_solvers::condensate_pipe_PP_task_boundaries_t::default_values();
    initial_conditions.pressure_in = 5e6;
    initial_conditions.pressure_out = 4e6;
    initial_conditions.density = 850.0;

    // Act
    task.solve(initial_conditions);
    auto& layer = task.get_current_layer();

    // Assert
    // Проверяем инициализацию плотности
    for (const auto& density : layer.density) {
        EXPECT_NEAR(density, 850.0, 1e-6);
    }

    // Проверяем, что расход был рассчитан
    EXPECT_FALSE(std::isnan(layer.std_volumetric_flow));
    EXPECT_GT(layer.std_volumetric_flow, 0.0);

    // Проверяем расчет давления
    ASSERT_FALSE(layer.pressure.empty());
    EXPECT_NEAR(layer.pressure.front(), 5e6, 1e-6);
    EXPECT_NEAR(layer.pressure.back(), 4e6, 1e-6);

    // Давление должно уменьшаться вдоль трубы
    for (size_t i = 1; i < layer.pressure.size(); ++i) {
        EXPECT_LT(layer.pressure[i], layer.pressure[i - 1]);
    }
}

/// @brief 9(можно удалить). Тест на зависимость расхода от перепада давления
TEST_F(CondensatePPTest, FlowRateVsPressureDrop) {
    // Arrange
    pde_solvers::condensate_pipe_PP_task_t task(pipe);
    auto base_conditions = pde_solvers::condensate_pipe_PP_task_boundaries_t::default_values();
    base_conditions.pressure_in = 5e6;
    base_conditions.density = 850.0;

    std::vector<double> pressure_drops = { 0.5e6, 1e6, 2e6, 3e6 };
    std::vector<double> calculated_flows;

    // Act
    for (double pressure_drop : pressure_drops) {
        base_conditions.pressure_out = base_conditions.pressure_in - pressure_drop;
        task.solve(base_conditions);
        calculated_flows.push_back(task.get_current_layer().std_volumetric_flow);
    }

    // Assert
    // Расход должен увеличиваться с увеличением перепада давления
    for (size_t i = 1; i < calculated_flows.size(); ++i) {
        EXPECT_GT(calculated_flows[i], calculated_flows[i - 1]);
    }

    // Для турбулентного режима расход ~ sqrt(ΔP)
    for (size_t i = 0; i < pressure_drops.size(); ++i) {
        double expected_ratio = std::sqrt(pressure_drops[i] / pressure_drops[0]);
        double actual_ratio = calculated_flows[i] / calculated_flows[0];
        EXPECT_NEAR(actual_ratio, expected_ratio, expected_ratio * 0.3);
    }
}

/// @brief 10(можно удалить). Тест на взаимную проверку PP и PQ задач
TEST_F(CondensatePPTest, PPandPQMutualVerification) {
    // Arrange
    pde_solvers::condensate_pipe_PP_task_t pp_task(pipe);
    pde_solvers::condensate_pipe_PQ_task_t pq_task(pipe);

    auto pp_conditions = pde_solvers::condensate_pipe_PP_task_boundaries_t::default_values();
    pp_conditions.pressure_in = 5e6;
    pp_conditions.pressure_out = 4e6;
    pp_conditions.density = 850.0;

    // Act: 1. Рассчитываем PP задачу
    pp_task.solve(pp_conditions);
    double calculated_flow = pp_task.get_current_layer().std_volumetric_flow;

    // Act: 2. Используем найденный расход в PQ задаче
    auto pq_conditions = pde_solvers::condensate_pipe_PQ_task_boundaries_t::default_values();
    pq_conditions.pressure_in = pp_conditions.pressure_in;
    pq_conditions.volumetric_flow = calculated_flow;
    pq_conditions.density = pp_conditions.density;

    pq_task.solve(pq_conditions);
    double pq_pressure_out = pq_task.get_current_layer().pressure.back();

    double relative = abs(pq_pressure_out - pp_conditions.pressure_out) / pp_conditions.pressure_out;
    // Assert: Давление на выходе должно совпадать
    EXPECT_NEAR(relative, 0, 1e-6);
}

/// @brief 11. Тест на обратную проверку PQ -> PP
TEST_F(CondensatePPTest, PQandPPMutualVerificationReverse) {
    // Arrange
    pde_solvers::condensate_pipe_PQ_task_t pq_task(pipe);
    pde_solvers::condensate_pipe_PP_task_t pp_task(pipe);

    auto pq_conditions = pde_solvers::condensate_pipe_PQ_task_boundaries_t::default_values();
    pq_conditions.pressure_in = 5e6;
    pq_conditions.volumetric_flow = 0.3;
    pq_conditions.density = 850.0;

    // Act: 1. Рассчитываем PQ задачу
    pq_task.solve(pq_conditions);
    double calculated_pressure_out = pq_task.get_current_layer().pressure.back();

    // Act: 2. Используем найденное давление в PP задаче
    auto pp_conditions = pde_solvers::condensate_pipe_PP_task_boundaries_t::default_values();
    pp_conditions.pressure_in = pq_conditions.pressure_in;
    pp_conditions.pressure_out = calculated_pressure_out;
    pp_conditions.density = pq_conditions.density;

    pp_task.solve(pp_conditions);
    double pp_flow = pp_task.get_current_layer().std_volumetric_flow;
    double relative = abs(pp_flow - pq_conditions.volumetric_flow) / pq_conditions.volumetric_flow;
    // Assert: Расход должен совпадать
    EXPECT_NEAR(relative, 0, 1e-5);
}

/// @brief 12.(удалить) Тест на согласованность PP и PQ при изменении параметров
TEST_F(CondensatePPTest, PPandPQConsistencyParameterSweep) {
    // Arrange
    std::vector<double> test_flows = { 0.1, 0.2, 0.3, 0.4 };
    std::vector<double> test_densities = { 700.0, 850.0, 1000.0 };
    const double pressure_in = 5e6;
    const double tolerance = 0.02; // 2%

    // Act & Assert
    for (double flow : test_flows) {
        for (double density : test_densities) {
            // Arrange для каждой итерации
            pde_solvers::condensate_pipe_PQ_task_t pq_task(pipe);
            pde_solvers::condensate_pipe_PP_task_t pp_task(pipe);

            auto pq_conditions = pde_solvers::condensate_pipe_PQ_task_boundaries_t::default_values();
            pq_conditions.pressure_in = pressure_in;
            pq_conditions.volumetric_flow = flow;
            pq_conditions.density = density;

            // Act: PQ расчет
            pq_task.solve(pq_conditions);
            double pq_pressure_out = pq_task.get_current_layer().pressure.back();

            auto pp_conditions = pde_solvers::condensate_pipe_PP_task_boundaries_t::default_values();
            pp_conditions.pressure_in = pressure_in;
            pp_conditions.pressure_out = pq_pressure_out;
            pp_conditions.density = density;

            // Act: PP расчет
            pp_task.solve(pp_conditions);
            double pp_flow = pp_task.get_current_layer().std_volumetric_flow;

            // Assert
            EXPECT_NEAR(pp_flow, flow, flow * tolerance)
                << "Flow=" << flow << " m³/s, Density=" << density << " kg/m³";
        }
    }
}

/// @brief 13.(Удалить) Тест на влияние шага по времени
TEST_F(CondensatePPTest, TimeStepSensitivity) {
    // Arrange
    pde_solvers::condensate_pipe_PP_task_t base_task(pipe);
    auto initial_conditions = pde_solvers::condensate_pipe_PP_task_boundaries_t::default_values();
    initial_conditions.pressure_in = 5e6;
    initial_conditions.pressure_out = 4e6;
    initial_conditions.density = 850.0;

    base_task.solve(initial_conditions);
    double initial_flow = base_task.get_current_layer().std_volumetric_flow;

    std::vector<double> time_steps = { 0.1, 1.0, 10.0 };
    std::vector<double> flow_changes;

    // Act
    for (double dt : time_steps) {
        pde_solvers::condensate_pipe_PP_task_t task(pipe);
        task.solve(initial_conditions);
        task.step(dt, initial_conditions);

        double new_flow = task.get_current_layer().std_volumetric_flow;
        double flow_change = std::abs(new_flow - initial_flow) / initial_flow;
        flow_changes.push_back(flow_change);
    }

    // Assert
    for (double flow_change : flow_changes) {
        EXPECT_LT(flow_change, 0.01);
    }

    for (size_t i = 1; i < flow_changes.size(); ++i) {
        EXPECT_NEAR(flow_changes[i], flow_changes[0], 1e-6);
    }
}

/// @brief 14.Тест на численную устойчивость PP задачи
TEST_F(CondensatePPTest, NumericalStabilityPP) {
    // Arrange
    pde_solvers::condensate_pipe_PP_task_t task(pipe);
    auto initial_conditions = pde_solvers::condensate_pipe_PP_task_boundaries_t::default_values();
    initial_conditions.pressure_in = 5e6;
    initial_conditions.pressure_out = 4.5e6;
    initial_conditions.density = 850.0;

    const int steps = 10;
    const double dt = 10.0;
    std::vector<double> flow_history;

    // Act: Начальный расчет
    task.solve(initial_conditions);
    flow_history.push_back(task.get_current_layer().std_volumetric_flow);

    // Act: Много шагов по времени
    for (int i = 0; i < steps; ++i) {
        task.step(dt, initial_conditions);
        flow_history.push_back(task.get_current_layer().std_volumetric_flow);
    }

    // Assert
    // Проверяем, что расход стабилизировался
    double max_variation = 0.0;
    for (int i = 1; i < flow_history.size(); ++i) {
        max_variation = std::max(max_variation,
            std::abs(flow_history[i] - flow_history[i - 1]));
    }

    EXPECT_LT(max_variation, 1e-6);
}

/// @brief 15. Тест на влияние плотности на расход при фиксированном перепаде
TEST_F(CondensatePPTest, DensityEffectOnFlowRate) {
    // Arrange
    pde_solvers::condensate_pipe_PP_task_t task(pipe);
    auto base_conditions = pde_solvers::condensate_pipe_PP_task_boundaries_t::default_values();
    base_conditions.pressure_in = 5e6;
    base_conditions.pressure_out = 4e6;

    std::vector<double> densities = { 700.0, 850.0, 1000.0 };
    std::vector<double> calculated_flows;

    // Act
    for (double density : densities) {
        base_conditions.density = density;
        task.solve(base_conditions);
        calculated_flows.push_back(task.get_current_layer().std_volumetric_flow);
    }

    // Assert
    // Расход должен уменьшаться с увеличением плотности
    for (size_t i = 1; i < calculated_flows.size(); ++i) {
        EXPECT_LT(calculated_flows[i], calculated_flows[i - 1]);
    }
}

/// @brief 16. Тест на крайние случаи PP задачи
TEST_F(CondensatePPTest, EdgeCasesPP_ZeroPressureDrop) {
    // Arrange
    pde_solvers::condensate_pipe_PP_task_t task(pipe);
    auto conditions = pde_solvers::condensate_pipe_PP_task_boundaries_t::default_values();
    conditions.pressure_in = 5e6;
    conditions.pressure_out = 5e6; 
    conditions.density = 850.0;

    // Act
    task.solve(conditions);
    auto& layer = task.get_current_layer();

    // Assert
    EXPECT_EQ(layer.std_volumetric_flow, 0);
}

/// @brief 17. (удалить) Тест на крайние случаи PP задачи - большой перепад
TEST_F(CondensatePPTest, EdgeCasesPP_LargePressureDrop) {
    // Arrange
    pde_solvers::condensate_pipe_PP_task_t task(pipe);
    auto conditions = pde_solvers::condensate_pipe_PP_task_boundaries_t::default_values();
    conditions.pressure_in = 10e6;
    conditions.pressure_out = 1e6; // 9 МПа перепада
    conditions.density = 850.0;

    // Act
    task.solve(conditions);
    auto& layer = task.get_current_layer();

    // Assert
    EXPECT_GT(layer.std_volumetric_flow, 0.0);

    for (size_t i = 1; i < layer.pressure.size(); ++i) {
        EXPECT_LT(layer.pressure[i], layer.pressure[i - 1]);
    }
}

/// @brief 18.(удалить) Тест на некорректные входные данные PP задачи
TEST_F(CondensatePPTest, EdgeCasesPP_InvalidDensity) {
    // Arrange
    pde_solvers::condensate_pipe_PP_task_t task(pipe);
    auto conditions = pde_solvers::condensate_pipe_PP_task_boundaries_t::default_values();
    conditions.pressure_in = 5e6;
    conditions.pressure_out = 4e6;
    conditions.density = 0.0; // Нулевая плотность

    // Act & Assert
    EXPECT_ANY_THROW(task.solve(conditions));
}

/// @brief 19.(удалить) Тест на метод Ньютона
TEST_F(CondensatePPTest, NewtonMethodConvergence) {
    // Arrange
    auto conditions = pde_solvers::condensate_pipe_PP_task_boundaries_t::default_values();
    conditions.pressure_in = 5e6;
    conditions.pressure_out = 4e6;
    conditions.density = 850.0;

    // Act & Assert для разных сценариев
    {
        // Тест 1: Нормальные условия
        pde_solvers::condensate_pipe_PP_task_t task1(pipe);
        EXPECT_NO_THROW(task1.solve(conditions));

        auto& layer1 = task1.get_current_layer();
        EXPECT_GT(layer1.std_volumetric_flow, 0.0);
        EXPECT_LT(layer1.std_volumetric_flow, 10.0);
    }

    {
        // Тест 2: Близкие давления
        auto conditions2 = conditions;
        conditions2.pressure_out = 4.9e6;

        pde_solvers::condensate_pipe_PP_task_t task2(pipe);
        EXPECT_NO_THROW(task2.solve(conditions2));

        auto& layer2 = task2.get_current_layer();
        EXPECT_GT(layer2.std_volumetric_flow, 0.0);
    }
}

/// @brief 20.(похож на предыдущий, удалить) Тест на физическую корректность результатов
TEST_F(CondensatePPTest, PhysicalCorrectness) {
    // Arrange
    pde_solvers::condensate_pipe_PP_task_t task(pipe);
    auto conditions = pde_solvers::condensate_pipe_PP_task_boundaries_t::default_values();
    conditions.pressure_in = 5e6;
    conditions.pressure_out = 3e6;
    conditions.density = 850.0;

    // Act
    task.solve(conditions);
    auto& layer = task.get_current_layer();

    // Assert
    // 1. Расход должен быть положительным
    EXPECT_GT(layer.std_volumetric_flow, 0.0);

    // 2. Давление должно монотонно уменьшаться
    for (size_t i = 1; i < layer.pressure.size(); ++i) {
        EXPECT_LT(layer.pressure[i], layer.pressure[i - 1]);
    }

    // 3. Все давления должны быть положительными
    for (const auto& pressure : layer.pressure) {
        EXPECT_GT(pressure, 0.0);
    }

    // 4. Плотности должны быть положительными
    for (const auto& density : layer.density) {
        EXPECT_GT(density, 0.0);
    }
}

/// @brief 21. (удалить) Тест на производительность PP задачи
TEST_F(CondensatePPTest, PerformanceBenchmark) {
    // Arrange
    auto conditions = pde_solvers::condensate_pipe_PP_task_boundaries_t::default_values();
    conditions.pressure_in = 5e6;
    conditions.pressure_out = 4e6;
    conditions.density = 850.0;

    const int iterations = 100;
    auto start_time = std::chrono::high_resolution_clock::now();

    // Act
    for (int i = 0; i < iterations; ++i) {
        pde_solvers::condensate_pipe_PP_task_t task(pipe);
        task.solve(conditions);
        // Просто проверяем, что расчет выполняется без ошибок
        EXPECT_GT(task.get_current_layer().std_volumetric_flow, 0.0);
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(
        end_time - start_time).count();

    // Assert
    double average_time = static_cast<double>(duration) / iterations;
    std::cout << "Average calculation time: " << average_time << " ms" << std::endl;

    // Для длинных труб или сложных расчетов устанавливаем разумный лимит
    EXPECT_LT(average_time, 1000.0); // Менее 1 секунды на расчет
}

/// @brief 22.(удалить) Тест с разными шагами по времени
TEST_F(CondensatePPTest, DifferentTimeSteps) {

    // Arrange
    pde_solvers::condensate_pipe_PQ_task_boundaries_t pq_conditions;
    pq_conditions.pressure_in = 5e6;
    pq_conditions.volumetric_flow = 0.3;
    pq_conditions.density = 850.0;

    
    // Разные шаги по времени
    std::vector<double> time_steps = { 1.0, 5.0, 10.0, 30.0, 60.0 };

    for (double dt : time_steps) {

        // Act 
        // 1. PQ расчет
        pde_solvers::condensate_pipe_PQ_task_t pq_task(pipe);
        pq_task.solve(pq_conditions);
        double initial_pressure_out = pq_task.get_current_layer().pressure.back();

        // 2. Один шаг PQ
        pq_task.step(dt, pq_conditions);
        double pressure_after_step = pq_task.get_current_layer().pressure.back();

        // 3. PP расчет с полученным давлением
        pde_solvers::condensate_pipe_PP_task_boundaries_t pp_conditions;
        pp_conditions.pressure_in = pq_conditions.pressure_in;
        pp_conditions.pressure_out = pressure_after_step;
        pp_conditions.density = pq_conditions.density;

        pde_solvers::condensate_pipe_PP_task_t pp_task(pipe);
        pp_task.solve(pp_conditions);
        double pp_flow = pp_task.get_current_layer().std_volumetric_flow;

        // 4. Проверка
        double error = abs(pq_conditions.volumetric_flow - pp_flow) / pq_conditions.volumetric_flow; 

        // Assert
        EXPECT_LT(error, 1e-4);
    }
}
}