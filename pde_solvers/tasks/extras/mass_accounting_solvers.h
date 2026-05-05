#pragma once

namespace pde_solvers {
;

/// @brief Свойства конденсатопровода с поддержкой расчета массы
struct iso_nonbaro_pipe_mass_accounting_properties_t : public iso_nonbaro_pipe_properties_t {
    /// @brief Флаг расчета массы вещества в трубе.
    bool calculate_mass{ false };
    /// @brief Создание структуры со значениями по умолчанию
    static iso_nonbaro_pipe_mass_accounting_properties_t default_values() {
        return iso_nonbaro_pipe_mass_accounting_properties_t(
            iso_nonbaro_pipe_properties_t::default_values());
    }
    /// @brief Конструктор из базовых свойств трубы
    explicit iso_nonbaro_pipe_mass_accounting_properties_t(
        const iso_nonbaro_pipe_properties_t& base_properties)
        : iso_nonbaro_pipe_properties_t(base_properties)
    {
    }
    /// @brief Конструктор из JSON данных
    /// Подразумевается, что при загрузке параметров объекта поле calculate_mass берется из таска.
    /// Поэтому переиспользуем конструктор предка.
    iso_nonbaro_pipe_mass_accounting_properties_t(const pde_solvers::pipe_json_data& json_data)
        : iso_nonbaro_pipe_properties_t(json_data)
    {
    }
};

/// @brief Расширенный расчетный слой с учетом расчета массы в трубе
struct iso_nonbaro_pipe_mass_accounting_layer_t : public iso_nonbaro_pipe_layer_t {
    /// @brief Масса вещества в трубе, кг
    double mass{ std::numeric_limits<double>::quiet_NaN() };
    /// @brief Инициализация слоя
    iso_nonbaro_pipe_mass_accounting_layer_t(size_t point_count)
        : iso_nonbaro_pipe_layer_t(point_count)
    {
    }
};

/// @brief Расширенный солвер с учетом расчета массы в трубе
class iso_nonbaro_pipe_mass_accounting_solver_t
    : public iso_nonbaro_pipe_solver_templated_t<iso_nonbaro_pipe_mass_accounting_layer_t> {
public:
    /// @brief Тип слоя (наследник от слоя основного солвера)
    using layer_type = iso_nonbaro_pipe_mass_accounting_layer_t;
    /// @brief Тип буфера
    using buffer_type = typename iso_nonbaro_pipe_solver_templated_t<layer_type>::buffer_type;
    /// @brief Тип параметров трубы
    using pipe_parameters_type = iso_nonbaro_pipe_mass_accounting_properties_t;

private:
    using base_solver_type = iso_nonbaro_pipe_solver_templated_t<layer_type>;
    /// @brief Параметры трубы с флагом расчета массы
    const pipe_parameters_type& mass_accounting_pipe;

public:
    /// @brief Фиктивный констуктор для совместмости с селектором рассчитываемых свойств
    iso_nonbaro_pipe_mass_accounting_solver_t(
        const pipe_parameters_type& pipe,
        buffer_type& buffer,
        const pde_solvers::endogenous_selector_t& endogenous_selector)
        : base_solver_type(pipe, buffer, endogenous_selector)
        , mass_accounting_pipe(pipe)
    {
    }

    /// @brief Конструктор
    iso_nonbaro_pipe_mass_accounting_solver_t(
        const pipe_parameters_type& pipe,
        buffer_type& buffer)
        : base_solver_type(pipe, buffer)
        , mass_accounting_pipe(pipe)
    {
    }

    /// @brief Вычисляет массу вещества в трубе по текущему слою
    double calculate_mass() const {
        const layer_type& current_layer = buffer.current();
        return calculate_isothermal_incompressible_fluid_mass_in_rigid_pipe(
            mass_accounting_pipe.profile.coordinates,
            current_layer.density_std.value,
            mass_accounting_pipe.wall.getArea());
    }

    /// @brief Выполнение транспортного шага (расчет движения партий)
    virtual void transport_step(double dt, double volumetric_flow, 
        const pde_solvers::endogenous_values_t& boundaries) override 
    {
        base_solver_type::transport_step(dt, volumetric_flow, boundaries);
        if (mass_accounting_pipe.calculate_mass) {
            buffer.current().mass = calculate_mass();
        }
    }

    /// @brief Транспортное решение при бесконечном dt (заполнение трубы граничными значениями)
    virtual void transport_solve(double volumetric_flow, 
        const pde_solvers::endogenous_values_t& boundaries) override 
    {
        base_solver_type::transport_solve(volumetric_flow, boundaries);
        if (mass_accounting_pipe.calculate_mass) {
            buffer.current().mass = calculate_mass();
        }
    }
};

/// @brief Расширенный расчетный слой с учетом массы для небаротропной трубы с присадкой
struct iso_nonbaro_improver_pipe_mass_accounting_layer_t : public iso_nonbaro_improver_pipe_layer_t {
    /// @brief Масса вещества в трубе, кг
    double mass{ std::numeric_limits<double>::quiet_NaN() };
    /// @brief Инициализация слоя
    iso_nonbaro_improver_pipe_mass_accounting_layer_t(size_t point_count)
        : iso_nonbaro_improver_pipe_layer_t(point_count)
    {
    }
};

/// @brief Свойства конденсатопровода с присадкой и поддержкой расчета массы
struct iso_nonbaro_improver_pipe_mass_accounting_properties_t : public iso_nonbaro_improver_pipe_properties_t {
    /// @brief Флаг расчета массы вещества в трубе.
    bool calculate_mass{ false };
    /// @brief Конструктор из JSON данных
    /// Подразумевается, что при загрузке параметров объекта поле calculate_mass берется из таска.
    /// Поэтому переиспользуем конструктор предка.
    iso_nonbaro_improver_pipe_mass_accounting_properties_t(const pde_solvers::pipe_json_data& json_data)
        : iso_nonbaro_improver_pipe_properties_t(json_data)
    {
    }
};

/// @brief Расширенный солвер для небаротропной трубы с присадкой и учетом массы
class iso_nonbaro_improver_pipe_mass_accounting_solver_t
    : public iso_nonbaro_improver_pipe_solver_templated_t<iso_nonbaro_improver_pipe_mass_accounting_layer_t> {
public:
    /// @brief Тип слоя (наследник от слоя основного солвера)
    using layer_type = iso_nonbaro_improver_pipe_mass_accounting_layer_t;
    /// @brief Тип буфера
    using buffer_type = ring_buffer_t<iso_nonbaro_improver_pipe_mass_accounting_layer_t>;
    /// @brief Тип параметров трубы
    using pipe_parameters_type = iso_nonbaro_improver_pipe_mass_accounting_properties_t;

private:
    using base_solver_type = iso_nonbaro_improver_pipe_solver_templated_t<layer_type>;
    /// @brief Параметры трубы с флагом расчета массы
    const pipe_parameters_type& mass_accounting_pipe;

public:
    /// @brief Констуктор с селектором рассчитываемых свойств
    iso_nonbaro_improver_pipe_mass_accounting_solver_t(
        const pipe_parameters_type& pipe,
        buffer_type& buffer,
        const pde_solvers::endogenous_selector_t& endogenous_selector)
        : base_solver_type(pipe, buffer, endogenous_selector)
        , mass_accounting_pipe(pipe)
    {
    }

    /// @brief Конструктор
    iso_nonbaro_improver_pipe_mass_accounting_solver_t(
        const pipe_parameters_type& pipe,
        buffer_type& buffer)
        : base_solver_type(pipe, buffer)
        , mass_accounting_pipe(pipe)
    {
    }

    /// @brief Вычисляет массу вещества в трубе по текущему слою
    double calculate_mass() const {
        const layer_type& current_layer = buffer.current();
        return calculate_isothermal_incompressible_fluid_mass_in_rigid_pipe(
            mass_accounting_pipe.profile.coordinates,
            current_layer.density_std.value,
            mass_accounting_pipe.wall.getArea());
    }

    /// @brief Выполнение транспортного шага (расчет движения партий)
    virtual void transport_step(double dt, double volumetric_flow, 
        const pde_solvers::endogenous_values_t& boundaries) override {
        base_solver_type::transport_step(dt, volumetric_flow, boundaries);
        if (mass_accounting_pipe.calculate_mass) {
            buffer.current().mass = calculate_mass();
        }
    }

    /// @brief Транспортное решение при бесконечном dt (заполнение трубы граничными значениями)
    virtual void transport_solve(double volumetric_flow, 
        const pde_solvers::endogenous_values_t& boundaries) override {
        base_solver_type::transport_solve(volumetric_flow, boundaries);
        if (mass_accounting_pipe.calculate_mass) {
            buffer.current().mass = calculate_mass();
        }
    }
};

}
