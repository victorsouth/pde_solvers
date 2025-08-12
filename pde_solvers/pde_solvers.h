#pragma once

#define M_PI       3.14159265358979323846
#define M_G        9.81
#define KELVIN_OFFSET 273.15
#ifndef ATMOSPHERIC_PRESSURE 
#define ATMOSPHERIC_PRESSURE 101325.0
#endif
#define TECHNICAL_ATMOSPHERE 98097.0



#include <pde_solvers/timeseries.h>

#include "core/ring_buffer.h"
#include "core/differential_equation.h"
#include "core/profile_structures.h"

#include "solvers/moc_solver.h"
#include "solvers/ode_solver.h"
#include "solvers/godunov_solver.h"
#include "solvers/quick_solver.h"


#include "pipe/oil.h"
#include "pipe/pipe_hydraulic_computations.h"
#include "pipe/pipe_hydraulic_struct.h"
#include "pipe/pipe_hydraulic_pde.h"
#include "pipe/pipe_profile_utils.h"
#include "pipe/pipe_advection_pde.h"
#include "pipe/pipe_advection_solver.h"

#include "pipe/heat_transfer.h"
#include "pipe/pipe_heat_struct.h"
#include "pipe/pipe_dynamic_soil.h"
#include "pipe/pipe_dynamic_soil_multizone.h"
#include "pipe/pipe_heat_computations.h"
#include "pipe/pipe_heat_pde.h"
#include "pipe/pipe_heat_util.h"



#include "solvers/diffusion_solver.h" // нужно инклудить после объявления трубы и проч.

#include "tasks/pipe_heat_task_util.h"
#include "tasks/endogenous_values.h"

#include "tasks/isothermal_quasistatic_task.h"
#include "tasks/isothermal_quasistatic_ident.h"
#include "tasks/nonisothermal_quasistatic_task.h"
#include "tasks/nonisothermal_quasistatic_task_p.h"
#include "tasks/nonisothermal_quasistatic_ident.h"

namespace pde_solvers
{
;

/// @brief Интерпретирует степень достоверность как булев флаг достоверности
inline bool discriminate_confidence_level(double confidence_level)
{
    return confidence_level > 0.95; // с константой экспериментируем
}

/// @brief Расчетный слой и его код достоверности
struct confident_layer_t {
    /// @brief Сам расчетный параметр
    std::vector<double> value;
    /// @brief Код достоверности - считается численным методом, 
    /// поэтому не булевый а вещественный
    std::vector<double> confidence;
    /// @brief Принимает количество точек, инициализирует количество ячеек
    confident_layer_t(size_t point_count, double initial_value)
        : value(point_count - 1, initial_value)
        , confidence(point_count - 1, 0.0)
    { }
    bool is_confident_layer() const {
        for (double confidence_value : confidence) {
            if (discriminate_confidence_level(confidence_value) == false)
                return false;
        }
        return true;
    }
};

/// @brief Расчетные целевые профили по трубе 
struct pipe_endogenous_variable_layer_t
{
    /// @brief Номинальный объемный расход
    double std_volumetric_flow{ std::numeric_limits<double>::quiet_NaN() };
    /// @brief Массовый расход
    std::vector<double> mass_flow;
    /// @brief Скорость потока
    std::vector<double> velocity;
    /// @brief Профиль плотности
    confident_layer_t density;
    /// @brief Профиль вязкости (рабочая при изотермическом расчете)
    confident_layer_t viscosity;
    /// @brief Серосодержанние
    confident_layer_t sulfur;
    /// @brief Концентрация ПТП
    confident_layer_t improver;
    /// @brief Температура
    confident_layer_t temperature;
    /// @brief Вязкость при 0С (сортовая при неизотермическом расчете)
    confident_layer_t viscosity0;
    /// @brief Вязкость при 20С (сортовая при неизотермическом расчете)
    confident_layer_t viscosity20;
    /// @brief Вязкость при 50С (сортовая при неизотермическом расчете)
    confident_layer_t viscosity50;
    /// @brief Принимает количество точек, инициализирует количество ячеек
    pipe_endogenous_variable_layer_t(size_t point_count)
        : mass_flow(point_count - 1, 0.0)
        , velocity(point_count - 1, 0.0)
        , density(point_count - 1, 860)
        , viscosity(point_count - 1, 1e-6)
        , sulfur(point_count - 1, 1e-3)
        , improver(point_count - 1, 0.0)
        , temperature(point_count - 1, 300)
        , viscosity0(point_count - 1, 0.0)
        , viscosity20(point_count - 1, 0.20e-6)
        , viscosity50(point_count - 1, 0.50e-6)
    {

    }
};

/// @brief Сбрасывает код достоверности в "ложь"
inline void reset_confidence(pipe_endogenous_variable_layer_t* layer) {
    auto invalidate = [](confident_layer_t& parameter_layer) {
        std::fill(parameter_layer.confidence.begin(), parameter_layer.confidence.end(), false);
        };

    invalidate(layer->density);
    invalidate(layer->viscosity);
    invalidate(layer->sulfur);
    invalidate(layer->improver);
    invalidate(layer->temperature);
    invalidate(layer->viscosity0);
    invalidate(layer->viscosity20);
    invalidate(layer->viscosity50);
}

/// @brief Расчетный профиль на квикест
struct pipe_endogenous_calc_layer_t : public pipe_endogenous_variable_layer_t
{
    /// @brief Профиль вспомогательных расчетов для метода конечных объемов (и для вязкости, и для плотности)
    quickest_ultimate_fv_solver_traits<1>::specific_layer specific;
    /// @brief Инициализация профилей
    /// @param point_count Количество точек
    pipe_endogenous_calc_layer_t(size_t point_count)
        : pipe_endogenous_variable_layer_t(point_count)
        , specific(point_count)
    {
    }

    /// @brief Подготовка плотности для расчета методом конечных объемов   
    static quickest_ultimate_fv_wrapper<1> get_density_wrapper(
        pipe_endogenous_calc_layer_t& layer)
    {
        return quickest_ultimate_fv_wrapper<1>(layer.density.value, layer.specific);
    }
    static quickest_ultimate_fv_wrapper<1> get_density_confidence_wrapper(
        pipe_endogenous_calc_layer_t& layer)
    {
        return quickest_ultimate_fv_wrapper<1>(layer.density.confidence, layer.specific);
    }
    /// @brief Подготовка вязкости для расчета методом конечных объемов 
    static quickest_ultimate_fv_wrapper<1> get_viscosity_wrapper(
        pipe_endogenous_calc_layer_t& layer)
    {
        return quickest_ultimate_fv_wrapper<1>(layer.viscosity.value, layer.specific);
    }
    static quickest_ultimate_fv_wrapper<1> get_viscosity_confidence_wrapper(
        pipe_endogenous_calc_layer_t& layer)
    {
        return quickest_ultimate_fv_wrapper<1>(layer.viscosity.confidence, layer.specific);
    }
};


}


