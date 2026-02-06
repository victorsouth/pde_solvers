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

/// @brief PQ/QP задача на уравнении импульса, численный расчет Эйлером
/// @tparam PipeEquationType Тип уравнения импульса
/// @param pipe Параметры трубы для уравнения импульса
/// @param current_layer Текущий слой с гидравлическими параметрами и эндогенными параметрами
/// @param std_volumetric_flow Объемный расход
/// @param bound_pressure Граничное давление (граница зависит от euler_direction)
/// @param euler_direction Определяет сторону, откуда идет расчет и соответственно границу, 
/// для которой задано граничное давление. 
/// Если +1 расчет слева направо, задано P_in
/// Если -1 расчет справла налево, задано P_out
/// @return Давление на конечной границе
template <typename PipeEquationType>
double rigorous_impulse_solve_QP(
    const typename PipeEquationType::pipe_properties_type& pipe,
    typename PipeEquationType::layer_type& current_layer,
    double std_volumetric_flow, double bound_pressure, int euler_direction) {

    current_layer.std_volumetric_flow = std_volumetric_flow;

    PipeEquationType pipeModel(pipe, current_layer, std_volumetric_flow, euler_direction);
    std::vector<double>& p_profile = current_layer.pressure;
    if (euler_direction > 0) {
        solve_euler<1>(pipeModel, euler_direction, 
            bound_pressure /* имеет смысл входного давления*/, 
            &p_profile);
        return p_profile.back();
    }
    else {
        solve_euler<1>(pipeModel, euler_direction, 
            bound_pressure /* имеет смысл выходного давления */, 
            &p_profile);
        return p_profile.front();
    }
}

/// @brief Расчетчик PP задачи методом Ньютона 
/// поверх PQ задачи на уравнении импульса (в свою очередь рассчитываемой Эйлером)
/// @tparam PipeEquationType тип уравнения импульса (должен иметь pipe_properties_type и layer_type)
template <typename PipeEquationType>
class rigorous_impulse_solver_PP : public fixed_system_t<1> {
private:
    /// @brief Ссылка на свойства трубы
    const typename PipeEquationType::pipe_properties_type& pipe;
    /// @brief слой расчета
    typename PipeEquationType::layer_type& current_layer;
    /// @brief Давление на входе, Па
    double bound_pressure_in;
    /// @brief Давление на выходе, Па
    double bound_pressure_out;

public:
    /// @brief Конструктор класса для решения задачи PP методом Ньютона
    /// @param pipe Ссылка на свойства конденсатопровода
    /// @param current_layer Текущий расчетный слой
    /// @param bound_pressure_in Давление на входе, Па
    /// @param bound_pressure_out Давление на выходе, Па
    rigorous_impulse_solver_PP(const typename PipeEquationType::pipe_properties_type& pipe,
        typename PipeEquationType::layer_type& current_layer, 
        double bound_pressure_in, double bound_pressure_out)
        : pipe(pipe)
        , bound_pressure_in(bound_pressure_in)
        , bound_pressure_out(bound_pressure_out)
        , current_layer(current_layer)
    {
    }
    /// @brief Невязка по давлению как функция от расхода
    /// Метод Эйлера интегрирует против направления расхода (от выхода при Q>=0, от входа при Q<0)
    virtual double residuals(const double& std_volumetric_flow) {
        // Направление Эйлера противоположно направлению расхода 
        // (обязательно для расчета самотеков, для напорного течения хуже не будет)
        int euler_direction = (std_volumetric_flow >= 0) ? -1 : +1;
        if (euler_direction > 0) {
            double calc_pressure_out = rigorous_impulse_solve_QP<PipeEquationType>(pipe, current_layer,
                std_volumetric_flow, bound_pressure_in, euler_direction);
            return calc_pressure_out - bound_pressure_out;
        }
        else {
            double calc_pressure_in = rigorous_impulse_solve_QP<PipeEquationType>(pipe, current_layer,
                std_volumetric_flow, bound_pressure_out, euler_direction);
            return calc_pressure_in - bound_pressure_in;
        }
    }

    /// @brief переопределяем целевую функцию, чтобы был модуль невязок
    virtual double objective_function(const double& r) const override {
        return std::abs(r);
    }
    /// @brief Расчет PP задачи методом Ньютона
    /// Если численный результат result = nullptr, то кидаем исключение, когда/если метод не сойдется
    double solve(double volumetric_flow_initial, fixed_solver_result_t<1>* result = nullptr)
    {
        fixed_solver_parameters_t<1, 0, golden_section_search> parameters;
        parameters.residuals_norm = 0.1; // погрешность 0.1 Па
        parameters.argument_increment_norm = 0;
        parameters.residuals_norm_allow_early_exit = true;

        // Создание структуры для записи результатов расчета
        if (result == nullptr) {
            fixed_solver_result_t<1> result_carrier;
            fixed_newton_raphson<1>::solve_dense(*this, { volumetric_flow_initial },
                parameters, &result_carrier);

            if (result_carrier.result_code == numerical_result_code_t::Converged) {
                return result_carrier.argument;
            }
            else {
                throw std::runtime_error("Solve PP not converged");
            }
        }
        else {
            fixed_newton_raphson<1>::solve_dense(*this, { volumetric_flow_initial },
                parameters, result);
            return result->result_code == numerical_result_code_t::Converged
                ? result->argument
                : std::numeric_limits<double>::quiet_NaN();
        }
    }

};



}