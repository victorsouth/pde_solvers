#pragma once
//
//VectorXd rosenbrock_function_terms(const VectorXd& value)
//{
//    VectorXd result(2);
//    result(0) = 10 * (value(1) - value(0) * value(0));
//    result(1) = 1 - value(0);
//    return result;
//};
//
//double rosenbrock_function(const VectorXd& value)
//{
//    //return 100*sqr(value(1) - sqr(value(0))) + sqr(1-value(0));
//    return rosenbrock_function_terms(value).squaredNorm();
//}
//
//
//
///// @brief Нелинейная система уравнений
//struct equation_system_t {
//    /// @brief Параметр численного расчета производной
//    double epsilon;
//    /// @brief Способ расчета производной - односторонняя или двустороння разность
//    bool one_sided_derivatives;
//    /// @brief Используется при расчете якобиана
//    bool allow_derivative_dimesion_brutal_fix;
//    /// @brief Трудно понять, зачем это нужно. Признак какой-то?
//    bool performing_jacobian_numerical_calculation;
//
//public:
//    equation_system_t(double _epsilon = 1e-8);
//    virtual ~equation_system_t() = default;
//
//public:
//    virtual VectorXd residuals(const VectorXd& originals) = 0;
//    virtual double operator()(const VectorXd& originals);
//
//    virtual void modify_model(const VectorXd& originals);
//
//    virtual vector<Eigen::Triplet<double>> jacobian_column(const VectorXd& argument, size_t column);
//    virtual SparseMatrix<double> jacobian(const VectorXd& argument);
//
//    virtual vector<double> get_gains() const;
//    virtual size_t dimension() const = 0;
//    virtual size_t argument_dimension() const = 0;
//
//    virtual VectorXd estimation() const;
//    virtual void update_estimation(const VectorXd& argument);
//public:
//    virtual bool is_linear() const;
//    virtual void get_linear_problem(MatrixXd* A, VectorXd* b);
//
//    /// @brief Кастомный критерий завершения расчета
//    virtual bool custom_stop_criteria(const VectorXd& current_residuals);
//    /// @brief Кастомный способ исследования регулировка шага
//    /// (добавлен для расчетов алгоритмов ТГЦ)
//    /// По умолчанию ничего не делает
//    /// @param argument Текущий аргумент
//    /// @param argument_increment Приращение шага метода Ньютона-Рафсона
//    virtual void custom_line_research(
//        const VectorXd& argument, const VectorXd& argument_increment);
//
//};
//
//
//class scalar_function_t : public equation_system_t {
//private:
//public:
//    scalar_function_t() = default;
//    virtual double operator()(const VectorXd& originals) = 0;
//    virtual VectorXd residuals(const Eigen::VectorXd&);
//    virtual size_t dimension() const;
//    //virtual VectorXd gradient(const VectorXd& originals, double* function_value = nullptr) = 0;
//    virtual vector<double> get_gains() const;
//    virtual size_t argument_dimension() const = 0;
//    virtual VectorXd estimation() const;
//    virtual void update_estimation(const VectorXd& argument);
//};
//
//
///// @brief Настройка солвера для систем уравнений и для оптимизации (что не очень хорошо, разделить бы)
//struct equation_solver_parameters_t {
//    double line_search_max_step{ 1 }; // полезно уменьшать <1, чтобы не было зацикливания
//    double line_search_min_step{ 5e-2 };
//    double line_search_step_divider{ 1.5 };
//    //solver_constraints_t constraints;
//    VectorXd initial_argument;
//    bool allow_exceptions{ false };
//    bool allow_random_step{ false };
//    bool step_criteria_assuming_search_step{ false };
//    /// @brief использовать квадратичное программирование в задачах оптимизации
//    bool optimize_quadprog{ false };
//    bool argument_history{ false };
//    bool argument_increment_history{ false };
//    /// @brief Использовать относительное приращение при расчете минимального критерия выхода
//    bool argument_increment_relative{ false };
//    double argument_increment_norm{ 1e-4 };
//    size_t iteration_count{ 100 };
//};
//
//
///// @brief Результаты численного расчета
//struct equation_solver_result_t {
//
//    /// @brief Невязка
//    VectorXd residuals;
//    /// @brief Аргумент функции
//    VectorXd argument;
//
//    /// @brief Параметры, с которым расчет был запущен (зачем??? 5.12.2022)
//    equation_solver_parameters_t settings;
//
//    /// @brief Балл за расчет
//    int score{ 0 };
//    /// @brief Флаг наличия сходимости
//    bool converged{ false };
//    bool argument_increment_criteria{ false }; // минимальное приращение спуска
//    bool custom_stop_criteria{ false }; // критерий останова, зависящий от задачи
//    bool enough_iterations_criteria{ false }; // пройдено заданное количество итераций
//    bool line_search_criteria{ false }; // выход при неуспешном проведении линейного поиска
//    double argument_increment_metric{ 0 };
//
//    double solution_time{ std::numeric_limits<double>::infinity() };
//    double solution_metric{ std::numeric_limits<double>::infinity() };
//
//    bool has_any_stop_criteria() const;
//    wstring verbalize_steps() const;
//    /// @brief Вывод собранных данных по исследованию целевой функции
//    /// Выводит шаг и затем значения ц.ф. в диапазоне alpha [0, 1]
//    /// Формат CSV ";"
//    void print_target_function(std::wostream& s) const;
//
//    size_t get_minimum_step_count() const;
//};
//
//
//VectorXd optimize_gauss_newton(scalar_function_t& function, 
//    const VectorXd& initial_argument, 
//    equation_solver_parameters_t solver_parameters, 
//    equation_solver_result_t* result)
//{
//    size_t n = initial_argument.size();
//
//    if (n == 0) {
//        result->converged = true;
//        return initial_argument;
//    }
//
//    result->argument = initial_argument;
//    VectorXd& argument = result->argument;
//    VectorXd& r = result->residuals;
//    vector<double> gains = function.get_gains();
//    VectorXd argument_increment, search_direction;
//
//    r = function.residuals(argument);
//
//    //result->analysis.learning_curve.push_back(r.norm());
//
//    result->custom_stop_criteria = function.custom_stop_criteria(r);
//    if (result->custom_stop_criteria) {
//        result->converged = true;
//        return argument;
//    }
//
//    for (size_t iteration = 0; iteration < solver_parameters.iteration_count; ++iteration)
//    {
//        SparseMatrix<double> Jsparse = function.jacobian(argument);
//
//
// /*       if (solver_parameters.optimize_quadprog) {
//            // для сравнения
//            JacobiSVD<MatrixXd> svd_solver(J, ComputeThinU | ComputeThinV);
//            VectorXd search_direction_unconstrained = svd_solver.solve(-r);
//            callback.after_search_direction(iteration, &argument, &r, &Jsparse, &search_direction_unconstrained);
//
//            MatrixXd H = J.transpose() * J;
//            VectorXd f = r.transpose() * J;
//
//            MatrixXd A;  VectorXd b;
//            solver_parameters.constraints.get_inequalities_constraints(argument, &A, &b);
//
//            // см. формулы в "Идентификация - quadprog.docx"
//            MatrixXd A_increment = A;
//            VectorXd b_increment = -A * argument + b;
//            solve_quadprog_inequalities(H, f, A_increment, b_increment,
//                VectorXd::Zero(argument.size()), &search_direction);
//
//            trim_increment(solver_parameters.constraints.boundaries, search_direction);
//            trim_increment_relative(solver_parameters.constraints.relative_boundaries, search_direction);
//
//
//            callback.after_search_direction_trim(iteration, &argument, &r, &Jsparse, &search_direction);
//        }
//        else */
//
//        {
//            // Наш Якобиан J получен из ряда Тейлора r(x0 + dx) = r(x0) + J(x0)dx, 
//            // линеаризованная задача МНК выглядит так: ||r(x0) + J(x0)dx|| -> min по dx
//            // Алгоритм JacobiSVD подразумевает задачу  ||Y - Ax|| -> min по x
//            // В итоге, либо передавать -J, либо -r. Последнее вычислительно быстрее
//            JacobiSVD<MatrixXd> svd_solver(J, ComputeThinU | ComputeThinV);
//            search_direction = svd_solver.solve(-r);
//
//            //trim_increment_min(solver_parameters.constraints.minimum, argument, search_direction);
//            //trim_increment_max(solver_parameters.constraints.maximum, argument, search_direction);
//            //trim_increment(solver_parameters.constraints.boundaries, search_direction);
//            //trim_increment_relative(solver_parameters.constraints.relative_boundaries, search_direction);
//
//        }
//
//        directed_function_t directed_function(function, argument, search_direction);
//        const double min_step = 1e-3;//IDENT_DERIVATIVE_INCREMENT;
//        bool found_result;
//        double search_step = line_search_zero_order_real_min(directed_function, 1.0, min_step, 2, &found_result, &callback);
//        if (!found_result) {
//            callback.after_criteria(true, false, false);
//            result->line_search_criteria = true;
//            result->converged = false;
//            break;
//        }
//
//        argument_increment = search_step * search_direction;
//        argument += argument_increment;
//
//        r = function.residuals(argument);
//
//
//        double argument_increment_metric = solver_parameters.step_criteria_assuming_search_step
//            ? argument_increment_factor(gains, argument_increment)
//            : argument_increment_factor(gains, search_direction);
//
//        result->argument_increment_criteria =
//            argument_increment_metric < solver_parameters.argument_increment_norm;
//        result->custom_stop_criteria = function.custom_stop_criteria(r);
//        result->enough_iterations_criteria = iteration + 1 == solver_parameters.iteration_count;
//
//        if (result->has_any_stop_criteria()) {
//            result->converged = true;
//            break;
//        }
//    }
//
//    return argument;
//}

#include <Eigen/Dense>
using namespace Eigen;

/// @brief Система алгебраических уравнений - базовый класс
//template <std::ptrdiff_t ResidualsDimension, std::ptrdiff_t ArgumentDimension>
class fixed_least_squares_function_t {
protected:
    /// @brief Относительное (!) приращение для расчета производных
    double epsilon{ 1e-6 };
public:
    //typedef typename fixed_system_types<ArgumentDimension>::var_type argument_type;
    //typedef typename fixed_system_types<ResidualsDimension>::var_type residuals_type;
    //typedef array<argument_type, Dimension> matrix_type;
    //typedef typename fixed_system_types<Dimension>::equation_coeffs_type matrix_value;

    typedef VectorXd argument_type;
    typedef VectorXd residuals_type;
    typedef MatrixXd matrix_type;

public:
    /// @brief Расчет целевой функции по аргументу
    double operator()(const argument_type& x) {
        residuals_type r = residuals(x);
        double sum_of_squares = r.squaredNorm();
        return sum_of_squares;
    }

    /// @brief Невязки системы уравнений
    virtual residuals_type residuals(const argument_type& x) = 0;
    /// @brief Якобиан системы уравнений
    virtual matrix_type jacobian_dense(const argument_type& x) {
        return jacobian_dense_numeric(x);
    }

    inline matrix_type jacobian_dense_numeric(const argument_type& x)
    {
        argument_type arg = x;

        matrix_type J;

        for (int arg_index = 0; arg_index < x.size(); ++arg_index) {
            double e = epsilon * std::max(1.0, abs(arg[arg_index]));
            arg[arg_index] = x[arg_index] + e;
            residuals_type f_plus = residuals(arg);
            arg[arg_index] = x[arg_index] - e;
            residuals_type f_minus = residuals(arg);
            arg[arg_index] = x[arg_index];

            residuals_type Jcol = (f_plus - f_minus) / (2 * e);

            if (arg_index == 0) {
                J = matrix_type(Jcol.size(), x.size());
            }
            for (size_t row = 0; row < static_cast<size_t>(x.size()); ++row) {
                J(row,arg_index) = Jcol[row];
            }
        }
        return J;
    }


    /// @brief Специфический критерий успешного завершения расчета
    /// @param r Текущее значения невязок
    /// @param x Текущее значение аргумента
    /// @return Флаг успешного завершения
    virtual bool custom_success_criteria(const residuals_type& r, const argument_type& x)
    {
        return true;
    }


};


//template <size_t Dimension>
class fixed_optimize_gauss_newton {
public:
    //typedef typename fixed_system_types<Dimension>::var_type var_type;
    //typedef typename fixed_system_types<Dimension>::right_party_type function_type;
    //typedef typename fixed_system_types<Dimension>::equation_coeffs_type equation_coeffs_type;

    typedef VectorXd argument_type;
    typedef VectorXd residuals_type;
    typedef MatrixXd matrix_type;
private:

    /// @brief Проверка значения на Nan/infinite для скалярного случая
    /// @param value Проверяемое значение
    /// @return true/false
    static inline bool has_not_finite(const double value)
    {
        if (!std::isfinite(value)) {
            return true;
        }
        return false;
    }

    /// @brief Проверка значения на Nan/infinite для векторного случая
    /// @param value Проверяемое значение
    /// @return true/false
    template <typename Container>
    static inline bool has_not_finite(const Container& values)
    {
        for (const double value : values) {
            if (has_not_finite(value))
                return true;
        }
        return false;
    }
private:
    static double argument_increment_factor(
        const argument_type& argument, const argument_type& argument_increment)
    {
        return argument_increment.norm() / argument.size();

    }
private:
    /// @brief Проведение процедуры линейного поиска по заданному алгоритму
    /// @tparam LineSearch Алгоритм линейного поиска
    /// @param line_search_parameters параметры линейного поиска
    /// @param residuals Невязки
    /// @param argument Текущий аргумент, относительного которого делается приращение
    /// @param r Текущая невязка в системе уравнений для быстрого расчета ц.ф.
    /// @param p Приращение аргумента
    /// @return Величина шага. По соглашению если алгоритм линейного поиска не сошелся, будет NaN
    template <typename LineSearch>
    static double perform_line_search(
        const typename LineSearch::parameters_type& line_search_parameters,
        fixed_least_squares_function_t& function,
        const argument_type& argument, const residuals_type& r, const argument_type& p)
    {
        auto directed_function = [&](double step) {
            return function(argument + step * p);
            };

        // Диапазон поиска, значения функции на границах диапазона
        double a = 0;
        double b = line_search_parameters.maximum_step;
        double function_a = directed_function(a);
        double function_b = directed_function(b);

        auto [search_step, elapsed_iterations] = LineSearch::search(
            line_search_parameters,
            directed_function, a, b, function_a, function_b);
        return search_step;
    }

public:
    /// @brief Запуск численного метода
    /// @tparam LineSearch Алгоритм регулировки шага поиска
    /// @param residuals Функция невязок
    /// @param initial_argument Начальное приближение
    /// @param solver_parameters Настройки поиска
    /// @param result Результаты расчета
    template <
        typename LineSearch = divider_search>
    static void optimize(
        fixed_least_squares_function_t& function,
        const argument_type& initial_argument,
        const fixed_solver_parameters_t<-1, 0, LineSearch>& solver_parameters,
        fixed_solver_result_t<-1>* result
    )
    {
        size_t n = initial_argument.size();

        if (n == 0) {
            result->result_code = numerical_result_code_t::Converged;
            return;
        }

        result->argument = initial_argument;
        VectorXd& argument = result->argument;
        VectorXd& r = result->residuals;
        
        VectorXd argument_increment, search_direction;

        r = function.residuals(argument);

        bool custom_stop_criteria = function.custom_success_criteria(r, argument);
        /*if (custom_stop_criteria) {
            result->result_code = numerical_result_code_t::Converged;
            return;
        }*/

        for (size_t iteration = 0; iteration < solver_parameters.iteration_count; ++iteration)
        {
            MatrixXd J = function.jacobian_dense(argument);

            /*       if (solver_parameters.optimize_quadprog) {
                       // для сравнения
                       JacobiSVD<MatrixXd> svd_solver(J, ComputeThinU | ComputeThinV);
                       VectorXd search_direction_unconstrained = svd_solver.solve(-r);
                       callback.after_search_direction(iteration, &argument, &r, &Jsparse, &search_direction_unconstrained);

                       MatrixXd H = J.transpose() * J;
                       VectorXd f = r.transpose() * J;

                       MatrixXd A;  VectorXd b;
                       solver_parameters.constraints.get_inequalities_constraints(argument, &A, &b);

                       // см. формулы в "Идентификация - quadprog.docx"
                       MatrixXd A_increment = A;
                       VectorXd b_increment = -A * argument + b;
                       solve_quadprog_inequalities(H, f, A_increment, b_increment,
                           VectorXd::Zero(argument.size()), &search_direction);

                       trim_increment(solver_parameters.constraints.boundaries, search_direction);
                       trim_increment_relative(solver_parameters.constraints.relative_boundaries, search_direction);


                       callback.after_search_direction_trim(iteration, &argument, &r, &Jsparse, &search_direction);
                   }
                   else */

            {
                // Наш Якобиан J получен из ряда Тейлора r(x0 + dx) = r(x0) + J(x0)dx, 
                // линеаризованная задача МНК выглядит так: ||r(x0) + J(x0)dx|| -> min по dx
                // Алгоритм JacobiSVD подразумевает задачу  ||Y - Ax|| -> min по x
                // В итоге, либо передавать -J, либо -r. Последнее вычислительно быстрее
                JacobiSVD<MatrixXd> svd_solver(J, ComputeThinU | ComputeThinV);
                search_direction = svd_solver.solve(-r);

                //trim_increment_min(solver_parameters.constraints.minimum, argument, search_direction);
                //trim_increment_max(solver_parameters.constraints.maximum, argument, search_direction);
                //trim_increment(solver_parameters.constraints.boundaries, search_direction);
                //trim_increment_relative(solver_parameters.constraints.relative_boundaries, search_direction);

            }


            double search_step = perform_line_search<LineSearch>(
                solver_parameters.line_search, function, argument, r, search_direction);

            if (!std::isfinite(search_step)) {
                result->result_code = numerical_result_code_t::LineSearchFailed;
                break;
            }

            argument_increment = search_step * search_direction;
            argument += argument_increment;

            r = function.residuals(argument);


            // Проверка критерия выхода по малому относительному приращению
            double argument_increment_metric = solver_parameters.step_criteria_assuming_search_step
                ? argument_increment_factor(argument, argument_increment)
                : argument_increment_factor(argument, search_direction);
            bool argument_increment_criteria =
                argument_increment_metric < solver_parameters.argument_increment_norm;
            bool custom_criteria = false;
            if (custom_criteria || argument_increment_criteria) {
                result->result_code = numerical_result_code_t::Converged;
                break;
            }

            //double argument_increment_metric = solver_parameters.step_criteria_assuming_search_step
            //    ? argument_increment_factor(gains, argument_increment)
            //    : argument_increment_factor(gains, search_direction);

            //result->argument_increment_criteria =
            //    argument_increment_metric < solver_parameters.argument_increment_norm;
            //result->custom_stop_criteria = function.custom_stop_criteria(r);
            //result->enough_iterations_criteria = iteration + 1 == solver_parameters.iteration_count;

            //if (result->has_any_stop_criteria()) {
            //    result->converged = true;
            //    break;
            //}
        }
    }

};



TEST(OptimiseGaussNewton, Test1)
{
    /// @brief J(x1, x2) = (x1 - 2)^2 + (x2 - 1)^2 
    class simple_sum_of_squares_function : public fixed_least_squares_function_t
    {        
    public:
        VectorXd residuals(const VectorXd& x) {
            VectorXd result(2);
            result[0] = x(0) - 2.0;
            result[1] = x(1) - 1.0;
            return result;
        }
    };

    simple_sum_of_squares_function function;
    fixed_solver_parameters_t<-1, 0, golden_section_search> parameters;
    fixed_solver_result_t<-1> result;

    VectorXd initial = VectorXd::Zero(2);

    fixed_optimize_gauss_newton::optimize(function, initial, parameters, &result);

    //fixed_newton_raphson<2>::solve_dense(test, { 0, 0 }, parameters, &result);

}
