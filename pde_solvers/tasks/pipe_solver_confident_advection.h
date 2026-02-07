#pragma once

namespace pde_solvers {
;

/// @brief Адвекционное решение при бесконечном dt (заполнение трубы граничным значением)
/// @tparam BufferType Тип буфера с методом get_buffer_wrapper
/// @tparam WrapperGetter Тип геттера обёртки слоя
template <typename BufferType, typename WrapperGetter>
void solve_advection(const endogenous_confident_value_t& boundary_value,
    BufferType& buffer,
    WrapperGetter value_getter,
    WrapperGetter confidence_getter)
{
    auto value_buffer_wrapper = buffer.get_buffer_wrapper(value_getter);
    std::vector<double>& val = value_buffer_wrapper.current().vars;
    std::fill(val.begin(), val.end(), boundary_value.value);

    auto confidence_buffer_wrapper = buffer.get_buffer_wrapper(confidence_getter);
    std::vector<double>& conf = confidence_buffer_wrapper.current().vars;
    std::fill(conf.begin(), conf.end(), boundary_value.confidence);
}

/// @brief Шаг адвекции одного параметра (конечный шаг по времени)
/// @tparam BufferType Тип буфера с методом get_buffer_wrapper
/// @tparam WrapperGetter Тип геттера обёртки слоя
template <typename BufferType, typename WrapperGetter>
void step_advection(double dt, double volumetric_flow,
    const endogenous_confident_value_t& boundary_value,
    double pipe_diameter,
    const std::vector<double>& profile_coordinates,
    BufferType& buffer,
    WrapperGetter value_getter,
    WrapperGetter confidence_getter)
{
    if (!std::isfinite(dt)) {
        throw std::runtime_error("step_advection: dt must be finite; use solve_advection for steady-state");
    }
    pipe_advection_pde_t advection_model(circle_area(pipe_diameter),
        volumetric_flow, profile_coordinates);
    {
        auto buffer_wrapper = buffer.get_buffer_wrapper(value_getter);
        quickest_ultimate_fv_solver solver(advection_model, buffer_wrapper);
        solver.step(dt, boundary_value.value, boundary_value.value);
    }
    {
        auto buffer_wrapper = buffer.get_buffer_wrapper(confidence_getter);
        quickest_ultimate_fv_solver solver(advection_model, buffer_wrapper);
        solver.step(dt, boundary_value.confidence, boundary_value.confidence);
    }
}

}
