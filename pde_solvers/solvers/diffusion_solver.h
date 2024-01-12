#pragma once

/// @brief Солвер физической диффузии при движении партий
/// Дидковская Новый метод расчета многопродуктовых магистральных трубопроводов 2018 ф-ла #7
class diffusion_transport_solver
{
private:
    const pipe_properties_t& pipe;
    const oil_parameters_t& oil;

public:
    diffusion_transport_solver(
        const pipe_properties_t& pipe,
        const oil_parameters_t& oil
    )
        : pipe(pipe)
        , oil(oil)
    {

    }
private:
    /// @brief Подинтегральное выражение
    /// Дидковская Новый метод расчета многопродуктовых магистральных
    /// трубопроводов 2018 ф-ла #7
    /// @param ts Верхний предел интеграла (безразмерное время)
    /// @param tau Текущее безразмерное время
    /// @param xs Безразмерная координата, для которой в конечном итоге будет построена функция параметра от безразмерного времени
    /// @param C_tau Значение параметра для текущего безразмерного времени
    /// @param Pe Число Пекле Pe = v*L/K
    /// @return Значение подинтегрального выражения
    static double function_under_integral2(double ts, double tau, double xs, double C_tau, double Pe)
    {
        //Исходная формула
        double R = (xs / (pow(ts - tau, 1.5))) * (exp((-Pe / 4) * (pow(xs - (ts - tau), 2) / (ts - tau)))) * C_tau;
        return R;
    }

    /// @brief Расчет величины целевого параметра для фиксированного момента времени для заданной координаты
    /// Численное интегрирование формулы (7) из 
    /// Дидковская Новый метод расчета многопродуктовых магистральных трубопроводов 2018
    /// @param ts Заданное время в безразмерной форме
    /// @param xs Заданная координата в безразмерной форме
    /// @param Pe Число Пекле
    /// @param h Шаг расчета интеграла по безразмерному времени
    /// @param N Количество точек краевых условий на входе, используемых в расчете
    /// @param input Краевые условия на входе
    /// @return Значение целевого параметра в точке (ts, xs)
    static double get_C_x_t2(double ts, double xs, double Pe, double h, size_t N, const vector<double>& input)
    {
        double C = 0;
        if (N > input.size())
        {
            throw std::runtime_error("Wrong input length");
        }

        for (size_t i = 0; i < N - 1; i++)
        {
            double tau = (i + 1) * h;
            double tau_p = (i)*h;
            double I2 = function_under_integral2(ts, tau, xs, input[i], Pe);
            double I1 = function_under_integral2(ts, tau_p, xs, input[i], Pe);
            C = C + (I2 + I1) * h / 2; // метод трапеций
            //C = C + I1 * h; // метод прямоугольников, только для 
        }
        C = C * sqrt(Pe) / (2 * sqrt(M_PI));
        return C;
    }


    static double calc_diffusive_transport(double t, double x, double delta_t,
        double v, double L, double K,
        const vector<double>& input)
    {
        double T = L / v;
        double Pe = v * L / K;

        double ts = t / T; // безразмерное время
        double xs = x / L; // безразмерная координата

        double delta_ts = delta_t / T; // шаг интегрирования в безразмерном времени

        size_t N = static_cast<size_t>(t / delta_t + 0.5);

        return get_C_x_t2(ts, xs, Pe, delta_ts, N, input);
    }


public:

    static vector<double> create_boundary(double initial_value, double final_value,
        size_t length, size_t start_change, size_t end_change)
    {
        vector<double> result(length);
        for (size_t index = 0; index < start_change; ++index) {
            result[index] = initial_value;
        }

        end_change = std::max(start_change, end_change);

        for (size_t index = start_change; index < end_change; ++index) {
            result[index] = initial_value + (final_value - initial_value) * (index - start_change) / (end_change - start_change);
        }

        for (size_t index = end_change; index < length; ++index) {
            result[index] = final_value;
        }
        return result;

    }



    static double calc_diffusion_coefficient(
        const pipe_properties_t& pipe, const oil_parameters_t& oil, double v
    )
    {
        double nu = oil.viscosity.nominal_viscosity;
        double Re = v * pipe.wall.diameter / nu;
        double relative_roughness = pipe.wall.relativeRoughness();
        double lambda = hydraulic_resistance_altshul(Re, relative_roughness);
        double K = 3.211 * sqrt(lambda) * v * pipe.wall.diameter;//0.1
        return K;
    }

    /// @brief Расчет временного ряда параметра на выходе трубопровода при движении партий
    /// Поскольку расчет идет очень медленно, предусмотрена возможность 
    /// задания произвольных моментов времени для выходных параметров (метод это позволяет)
    /// Коэффициент продольного перемешивания K зависит от скорости, поэтому считается тут внутри
    /// @param t_output Моменты времени, для которых считается параметр на выходе трубопровода
    /// @param delta_t Период дискретизации входного временного ряда, он же используется в интеграле
    /// @param input Временной ряд параметра на входе
    /// @param v Скорость потока
    /// @param use_offset_trick Использовать обход проблемы с нулевыми начальными условиями 
    /// (это корректно, т.к. модель линейна)
    /// @return Выходной временной ряд параметра для моментов времени t_output
    vector<double> solve(
        const vector<double>& t_output,
        double delta_t, vector<double> input,
        double v, bool use_offset_trick)
    {
        double pipe_length = pipe.profile.getLength();
        double K = calc_diffusion_coefficient(pipe, oil, v);

        double offset = 0;
        if (use_offset_trick) {
            offset = input[0];
            // переводим input в приращения относительно начального смещения
            for (double& in : input) {
                in -= offset;
            }
        }


        vector<double> output(t_output.size());

#pragma omp parallel for
        for (int i = 0; i < static_cast<int>(t_output.size()); i++)
        {
            double ti = t_output[i];
            double x = pipe_length;
            double value_at_ti = calc_diffusive_transport(ti, x, delta_t, v, pipe_length, K, input);
            output[i] = value_at_ti;
        }

        if (use_offset_trick) {
            // учитываем, что 0 по выходу соответствует величине offset
            for (double& out : output) {
                out += offset;
            }
        }

        return output;
    }

};
