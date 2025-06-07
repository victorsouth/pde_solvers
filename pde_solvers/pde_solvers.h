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

#include "solvers/diffusion_solver.h" // нужно инклудить после объявления трубы и проч.

#include "tasks/isothermal_quasistatic_task.h"
#include "tasks/isothermal_quasistatic_ident.h"

namespace pde_solvers
{
;

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
    std::vector<double> density;
    /// @brief Профиль вязкости (рабочая при изотермическом расчете)
    std::vector<double> viscosity;
    /// @brief Серосодержанние
    std::vector<double> sulfur;
    /// @brief Концентрация ПТП
    std::vector<double> improver;
    /// @brief Температура
    std::vector<double> temperature;
    /// @brief Вязкость при 0С (сортовая при неизотермическом расчете)
    std::vector<double> viscosity0;
    /// @brief Вязкость при 20С (сортовая при неизотермическом расчете)
    std::vector<double> viscosity20;
    /// @brief Вязкость при 50С (сортовая при неизотермическом расчете)
    std::vector<double> viscosity50;
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
    static quickest_ultimate_fv_wrapper<1> get_density_wrapper(pipe_endogenous_calc_layer_t& layer)
    {
        return quickest_ultimate_fv_wrapper<1>(layer.density, layer.specific);
    }
    /// @brief Подготовка вязкости для расчета методом конечных объемов 
    static quickest_ultimate_fv_wrapper<1> get_viscosity_wrapper(pipe_endogenous_calc_layer_t& layer)
    {
        return quickest_ultimate_fv_wrapper<1>(layer.viscosity, layer.specific);
    }
};


}


