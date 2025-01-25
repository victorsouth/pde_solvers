#pragma once

namespace pde_solvers {
;

/*
   Модели теплообмена для трубы без учета тепловой динамики грунта
   Реализован теплообмен через теплопроводность (твердые тела) 
     и лучистый перенос (вакуумная изоляция)
   Параметры трубы (в частности) уже заданы в PipeParameters
*/

/// @brief Теплоперенос через теплопроводность твердых тел. 
/// Изоляция однослойная (из одного материала), либо уже предварительно эквивалентированная
struct thermal_conductivity_parameters_t {
    /// @brief Толщина цилиндрической стенки от внутреннего диаметра 
    double thickness;
    /// @brief Теплопроводность цилиндрической стенки
    double conductivity;
    /// @brief Коэффициент теплопроводности UA
    /// @param diameter_inner Внутренний диаметр трубы
    double get_heat_transfer_coefficient(double diameter_inner) const {
        const double& lambda = conductivity;
        double d_in = diameter_inner;
        double d_out = d_in + 2 * thickness;
        double Kt = 2 * lambda / (d_in  * log(d_out / d_in));
        return Kt;
    }
    /// @brief Тепловой поток теплопроводности
    double get_heat_flow(double temperature_inner, 
        double temperature_outer, double diameter_inner) const 
    {
        const double Kt = get_heat_transfer_coefficient(diameter_inner);
        const double& Tin = temperature_inner;
        const double& Tout = temperature_outer;
        return -Kt * (Tin - Tout);
    }
};

/// @brief Параметры лучистого переноса тепла между концентрическими трубами
struct radiative_transport_parameters_t {
    /// @brief Диаметр внешней трубы
    double thickness;
    /// @brief Степень черноты внутренней трубы
    double emissivity_factor_inner{ 0.075 };
    /// @brief Степень черноты внешней трубы
    double emissivity_factor_outer{ 0.075 };
    /// @brief Приведенная степень черноты
    double get_reduced_emissivity_factor(double diameter_inner) const {
        double diameter_outer = diameter_inner + 2 * thickness;
        double denum =
            1 / emissivity_factor_inner +
            diameter_inner / diameter_outer * (1 / emissivity_factor_outer - 1);
        double result = 1 / denum;
        return result;
    }
    /// @brief Удельный тепловой поток
    /// @param temperature_inner 
    /// @return 
    double get_heat_flow(double temperature_inner,
        double temperature_outer, double diameter_inner) const 
    {
        const double& Tinner = temperature_inner;
        const double& Touter = temperature_outer;
        constexpr double boltzmann_constant = 5.76e-8;
        double eps_red = get_reduced_emissivity_factor(diameter_inner);
        double q = -eps_red * boltzmann_constant * (pow(Tinner, 4) - pow(Touter, 4));
        return q;
    }
};

/// @brief Базовая модель теплообмена для трубы. 
/// Характер теплообмен по длине трубы неизменен, зонирования нет!
struct heat_model_general_t {
    /// @brief Температура окружающей средой
    double ambient_temperature{ 273 };
    virtual double get_heat_flow(double T_fluid, double diameter_inner) const = 0;
};

/// @brief Комбинированная модель теплообмена: 
/// теплопроводность и лучистая
struct heat_model_conductivity_radiative_t : public heat_model_general_t {
    /// @brief Параметры теплопроводности
    thermal_conductivity_parameters_t conductivity;
    /// @brief Параметры лучистого теплообмена
    radiative_transport_parameters_t radiative;
    /// @brief Доля теплопроводности в общем теплообмене (из теплопроводности и лучистого)
    double conductivity_frac{ 0.04 };
    virtual double get_heat_flow(double T_fluid, double diameter_inner) const override {
        double q_cond = conductivity.get_heat_flow(T_fluid, ambient_temperature, diameter_inner);
        double q_rad = radiative.get_heat_flow(T_fluid, ambient_temperature, diameter_inner);
        double q = q_cond * conductivity_frac + q_rad * (1 - conductivity_frac);
        return q;
    }
};

}