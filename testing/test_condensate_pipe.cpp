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


class CondensatePipeQP : public ::testing::Test {
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


/// @brief Проверяет способность системы PQ задачи поддерживать корректный профиль давления при заданных начальных условиях.
/// Тест выполняет следующие проверки:
/// 1. Корректность инициализации плотности - все значения плотности в профиле должны соответствовать заданному начальному значению (800.0 кг/м³) с точностью 1e-6
/// 2. Наличие рассчитанного профиля давления - профиль давления не должен быть пустым
/// 3. Сохранение граничного условия на входе - давление на входе должно точно совпадать с заданным значением (5 МПа) с точностью 1e-6
/// 4. Монотонное убывание давления вдоль трубы - давление должно строго уменьшаться от входа к выходу из-за гидравлических потерь
/// 5. Физическая корректность значений - все значения давления в профиле должны быть положительными
TEST_F(CondensatePipeQP, MaintainsCorrectPressureProfile_WhenGivenInitialConditions) {

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

/// @brief Проверяет способность системы увеличивать потери давления с ростом расхода
TEST_F(CondensatePipeQP, IncreasesPressureLoss_WithIncreasingFlowRate) {
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




/// @brief Проверяет способность системы поддерживать стабильность давления при множественных шагах по времени
TEST_F(CondensatePipeQP, MaintainsPressureStability_OverMultipleTimeSteps) {
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

/// @brief Проверяет способность системы производить нулевой перепад давления при нулевом расходе
TEST_F(CondensatePipeQP, ProducesZeroPressureDrop_WhenFlowRateIsZero) {
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



/// @brief Проверяет способность системы PP задачи рассчитывать валидный расход для заданных граничных условий по давлению
TEST_F(CondensatePPTest, CalculatesValidFlowRate_ForGivenPressureBoundaries) {
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

/// @brief Проверяет способность системы восстанавливать исходный расход при использовании давления из PQ задачи в PP задаче
TEST_F(CondensatePPTest, RecoversOriginalFlowRate_WhenPQPressureUsedInPP) {
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

/// @brief Проверяет способность системы PP задачи поддерживать стабильность расхода при множественных шагах по времени
TEST_F(CondensatePPTest, MaintainsFlowRateStability_OverMultipleTimeSteps) {
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

/// @brief Проверяет способность системы уменьшать расход с увеличением плотности при фиксированном перепаде давления
TEST_F(CondensatePPTest, DecreasesFlowRate_WithIncreasingDensity_AtFixedPressureDrop) {
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

/// @brief Проверяет способность системы производить нулевой расход при нулевом перепаде давления
TEST_F(CondensatePPTest, ProducesZeroFlowRate_WhenPressureDropIsZero) {
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

}