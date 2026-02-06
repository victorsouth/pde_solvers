#pragma once

namespace pde_solvers {
;

/// @brief Площадь сечения трубы
inline double circle_area(double diameter) {
    return M_PI * diameter * diameter / 4;
}

/// @brief Гидравлическое сопротивление по Шифринсону
/// \param reynolds_number
/// \param relative_roughness
/// \return
inline double hydraulic_resistance_shifrinson(double reynolds_number, double relative_roughness)
{
    return 0.11 * pow(relative_roughness, 0.25);
}

/// @brief Гидравлическое сопротивление по Альтушулю
/// @param reynolds_number 
/// @param relative_roughness 
/// @return 
inline double hydraulic_resistance_altshul(double reynolds_number, double relative_roughness)
{
    return 0.11 * pow(relative_roughness + 68 / reynolds_number, 0.25);
}


/// @brief Расчет гидравлического сопротивления в широком диапазоне чисел Рейнольдса
/// Квадратичное трение по формуле Исаева [Морозова, Коршак, ф-ла (1)]
/// @param reynolds_number 
/// @param relative_roughness 
/// @return 
inline double hydraulic_resistance_isaev(double reynolds_number, double relative_roughness) {
    const double Re = fabs(reynolds_number);
    const double& Ke = relative_roughness;

    double lam;

    // Формулы по РД-75.180.00-КТН-258-10. 
    // В конце нет форулы для больших чисел Рейнольдса

    if (Re < 1) {
        lam = 64; // Стокс при Re = 1
    }
    else if (Re < 2320)
    {
        lam = 64 / Re; // Стокс 
    }
    else if (Re < 4000)
    {
        // Стокс + Блазиус, сглаженный переход
        double gm = 1 - exp(-0.002 * (Re - 2320));
        lam = 64 / Re * (1 - gm) + 0.3164 / pow(Re, 0.25) * gm;
    }
    //else if(Re < 1e4)//min(1e5,27/pow(Ke,1.143))  27/pow(Ke,1.143)   1e4   //10/Ke
    //{
    //    lam=0.3164/pow(Re,0.25);
    //}
    else if (Re < 560 / Ke)
    {
        // Исаев по [Морозова, Коршак], ф-ла (1)
        lam = 1.0 / fixed_solvers::sqr(-1.8 * log10(6.8 / Re + pow(Ke / 3.7, 1.1)));
    }
    else
    {
        // Шифринсон
        lam = 0.11 * pow(Ke, 0.25);
    }

    return lam;
}

/// @brief Стационарный расчет трубопровода по граничным давлениям для заданной системы 
/// уравнений методом Эйлера
/// @tparam PipeModel 
/// @param model 
/// @param Pin 
/// @param Pout 
/// @param layer 
/// @return 
template <typename PipeModel>
inline double solve_pipe_PP(PipeModel& model, double Pin, double Pout,
    profile_wrapper<double, 2>* layer)
{
    auto g = [&](double G)
    {
        solve_euler_corrector<2>(model, -1, { Pout, G }, layer);
        double Pin_calc = layer->profile(0).front();
        return Pin - Pin_calc;
    };
    fixed_scalar_wrapper_t f(g, 1e-3);

    fixed_solver_parameters_t<1, 0> parameters;
    parameters.constraints.relative_boundary = 50; // ограничение на шаг по расходу
    fixed_solver_result_t<1> result;
    fixed_newton_raphson<1>::solve_dense(f, { 0 }, parameters, &result);

    return result.argument;
}


/// @brief класс для нахождения расхода Q для задачи PP с помощью метода Ньютона
/// @tparam PipeEquationType тип уравнения импульса (iso_nonbaro_impulse_equation_t или iso_nonbaro_improver_impulse_equation_t)
/// @tparam BoundariesType класс граничных условий
/// @tparam LayerType класс уровней в buffer
template <typename PipeEquationType, typename BoundariesType, typename LayerType>
class impulse_solver_PP : public fixed_system_t<1> {
    /// @brief Тип переменной для системы уравнений
    using fixed_system_t<1>::var_type;
private:
    /// @brief ссылка на уравнение импульса (расход обновляется в residuals)
    PipeEquationType& pipe_model;
    /// @brief слой расчета
    LayerType& current_layer;
    /// @brief ГУ
    const BoundariesType& bound;

public:
    /// @brief Конструктор класса для решения задачи PP методом Ньютона
    /// @param pipe_model Ссылка на экземпляр уравнения импульса конденсатопровода
    /// @param bound Граничные условия задачи PP
    /// @param current_layer Текущий расчетный слой
    impulse_solver_PP(PipeEquationType& pipe_model,
        const BoundariesType& bound, LayerType& current_layer)
        : pipe_model(pipe_model)
        , bound(bound)
        , current_layer(current_layer)
    {
    }
    /// @brief Невязка по давлению как функция от расхода
    virtual double residuals(const double& std_volumetric_flow) {
        pipe_model.set_flow(std_volumetric_flow);
        std::vector<double>& p_profile = current_layer.pressure;

        int euler_direction = pipe_model.get_solver_direction();
        if (euler_direction > 0) {
            solve_euler<1>(pipe_model, euler_direction, bound.pressure_in, &p_profile);
            return p_profile.back() - bound.pressure_out;
        }
        else {
            solve_euler<1>(pipe_model, euler_direction, bound.pressure_out, &p_profile);
            return p_profile.front() - bound.pressure_in;
        }
    }

    /// @brief переопределяем целевую функцию, чтобы был модуль невязок
    virtual double objective_function(const var_type& r) const override {
        return std::abs(r);
    }
};



}