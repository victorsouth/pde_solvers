#pragma once
//
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

//!!!!!!!!!!!!!!!
///// @brief Система алгебраических уравнений - базовый класс
////template <std::ptrdiff_t ResidualsDimension, std::ptrdiff_t ArgumentDimension>
//class fixed_least_squares_function_t {
//protected:
//    /// @brief Относительное (!) приращение для расчета производных
//    double epsilon{ 1e-6 };
//public:
//    //typedef typename fixed_system_types<ArgumentDimension>::var_type argument_type;
//    //typedef typename fixed_system_types<ResidualsDimension>::var_type residuals_type;
//    //typedef array<argument_type, Dimension> matrix_type;
//    //typedef typename fixed_system_types<Dimension>::equation_coeffs_type matrix_value;
//
//    typedef VectorXd argument_type;
//    typedef VectorXd residuals_type;
//    typedef MatrixXd matrix_type;
//
//public:
//    /// @brief Расчет целевой функции по аргументу
//    double operator()(const argument_type& x) {
//        residuals_type r = residuals(x);
//        double sum_of_squares = r.squaredNorm();
//        return sum_of_squares;
//    }
//
//    /// @brief Невязки системы уравнений
//    virtual residuals_type residuals(const argument_type& x) = 0;
//    /// @brief Якобиан системы уравнений
//    virtual matrix_type jacobian_dense(const argument_type& x) {
//        return jacobian_dense_numeric(x);
//    }
//
//    inline matrix_type jacobian_dense_numeric(const argument_type& x)
//    {
//        argument_type arg = x;
//
//        matrix_type J;
//
//        for (int arg_index = 0; arg_index < x.size(); ++arg_index) {
//            double e = epsilon * std::max(1.0, abs(arg[arg_index]));
//            arg[arg_index] = x[arg_index] + e;
//            residuals_type f_plus = residuals(arg);
//            arg[arg_index] = x[arg_index] - e;
//            residuals_type f_minus = residuals(arg);
//            arg[arg_index] = x[arg_index];
//
//            residuals_type Jcol = (f_plus - f_minus) / (2 * e);
//
//            if (arg_index == 0) {
//                J = matrix_type(Jcol.size(), x.size());
//            }
//            for (size_t row = 0; row < static_cast<size_t>(Jcol.size()); ++row) {
//                J(row,arg_index) = Jcol[row];
//            }
//        }
//        return J;
//    }
//
//
//    /// @brief Специфический критерий успешного завершения расчета
//    /// @param r Текущее значения невязок
//    /// @param x Текущее значение аргумента
//    /// @return Флаг успешного завершения
//    virtual bool custom_success_criteria(const residuals_type& r, const argument_type& x)
//    {
//        return true;
//    }
//
//
//};

//!!!!!!!!!!!!!!!
//template <size_t Dimension>
//class fixed_optimize_gauss_newton {
//public:
//    //typedef typename fixed_system_types<Dimension>::var_type var_type;
//    //typedef typename fixed_system_types<Dimension>::right_party_type function_type;
//    //typedef typename fixed_system_types<Dimension>::equation_coeffs_type equation_coeffs_type;
//
//    typedef VectorXd argument_type;
//    typedef VectorXd residuals_type;
//    typedef MatrixXd matrix_type;
//private:
//
//    /// @brief Проверка значения на Nan/infinite для скалярного случая
//    /// @param value Проверяемое значение
//    /// @return true/false
//    static inline bool has_not_finite(const double value)
//    {
//        if (!std::isfinite(value)) {
//            return true;
//        }
//        return false;
//    }
//
//    /// @brief Проверка значения на Nan/infinite для векторного случая
//    /// @param value Проверяемое значение
//    /// @return true/false
//    template <typename Container>
//    static inline bool has_not_finite(const Container& values)
//    {
//        for (const double value : values) {
//            if (has_not_finite(value))
//                return true;
//        }
//        return false;
//    }
//private:
//    static double argument_increment_factor(
//        const argument_type& argument, const argument_type& argument_increment)
//    {
//        return argument_increment.norm() / argument.size();
//
//    }
//private:
//    /// @brief Проведение процедуры линейного поиска по заданному алгоритму
//    /// @tparam LineSearch Алгоритм линейного поиска
//    /// @param line_search_parameters параметры линейного поиска
//    /// @param residuals Невязки
//    /// @param argument Текущий аргумент, относительного которого делается приращение
//    /// @param r Текущая невязка в системе уравнений для быстрого расчета ц.ф.
//    /// @param p Приращение аргумента
//    /// @return Величина шага. По соглашению если алгоритм линейного поиска не сошелся, будет NaN
//    template <typename LineSearch>
//    static double perform_line_search(
//        const typename LineSearch::parameters_type& line_search_parameters,
//        fixed_least_squares_function_t& function,
//        const argument_type& argument, const residuals_type& r, const argument_type& p)
//    {
//        auto directed_function = [&](double step) {
//            return function(argument + step * p);
//            };
//
//        // Диапазон поиска, значения функции на границах диапазона
//        double a = 0;
//        double b = line_search_parameters.maximum_step;
//        double function_a = directed_function(a);
//        double function_b = directed_function(b);
//
//        auto [search_step, elapsed_iterations] = LineSearch::search(
//            line_search_parameters,
//            directed_function, a, b, function_a, function_b);
//        return search_step;
//    }
//
//public:
//
//    /// @brief Запуск численного метода
//    /// @tparam LineSearch Алгоритм регулировки шага поиска
//    /// @param residuals Функция невязок
//    /// @param initial_argument Начальное приближение
//    /// @param solver_parameters Настройки поиска
//    /// @param result Результаты расчета
//    template <
//        typename LineSearch = divider_search>
//    static void optimize(
//        fixed_least_squares_function_t& function,
//        const argument_type& initial_argument,
//        const fixed_solver_parameters_t<-1, 0, LineSearch>& solver_parameters,
//        fixed_solver_result_t<-1>* result
//    )
//    {
//        size_t n = initial_argument.size();
//
//        if (n == 0) {
//            result->result_code = numerical_result_code_t::Converged;
//            return;
//        }
//
//        result->argument = initial_argument;
//        VectorXd& argument = result->argument;
//        VectorXd& r = result->residuals;
//        
//        VectorXd argument_increment, search_direction;
//
//        r = function.residuals(argument);
//
//        bool custom_stop_criteria = function.custom_success_criteria(r, argument);
//        /*if (custom_stop_criteria) {
//            result->result_code = numerical_result_code_t::Converged;
//            return;
//        }*/
//
//        for (size_t iteration = 0; iteration < solver_parameters.iteration_count; ++iteration)
//        {
//            MatrixXd J = function.jacobian_dense(argument);
//
//            /*       if (solver_parameters.optimize_quadprog) {
//                       // для сравнения
//                       JacobiSVD<MatrixXd> svd_solver(J, ComputeThinU | ComputeThinV);
//                       VectorXd search_direction_unconstrained = svd_solver.solve(-r);
//                       callback.after_search_direction(iteration, &argument, &r, &Jsparse, &search_direction_unconstrained);
//
//                       MatrixXd H = J.transpose() * J;
//                       VectorXd f = r.transpose() * J;
//
//                       MatrixXd A;  VectorXd b;
//                       solver_parameters.constraints.get_inequalities_constraints(argument, &A, &b);
//
//                       // см. формулы в "Идентификация - quadprog.docx"
//                       MatrixXd A_increment = A;
//                       VectorXd b_increment = -A * argument + b;
//                       solve_quadprog_inequalities(H, f, A_increment, b_increment,
//                           VectorXd::Zero(argument.size()), &search_direction);
//
//                       trim_increment(solver_parameters.constraints.boundaries, search_direction);
//                       trim_increment_relative(solver_parameters.constraints.relative_boundaries, search_direction);
//
//
//                       callback.after_search_direction_trim(iteration, &argument, &r, &Jsparse, &search_direction);
//                   }
//                   else */
//
//            {
//                // Наш Якобиан J получен из ряда Тейлора r(x0 + dx) = r(x0) + J(x0)dx, 
//                // линеаризованная задача МНК выглядит так: ||r(x0) + J(x0)dx|| -> min по dx
//                // Алгоритм JacobiSVD подразумевает задачу  ||Y - Ax|| -> min по x
//                // В итоге, либо передавать -J, либо -r. Последнее вычислительно быстрее
//                JacobiSVD<MatrixXd> svd_solver(J, ComputeThinU | ComputeThinV);
//                search_direction = svd_solver.solve(-r);
//
//                //trim_increment_min(solver_parameters.constraints.minimum, argument, search_direction);
//                //trim_increment_max(solver_parameters.constraints.maximum, argument, search_direction);
//                //trim_increment(solver_parameters.constraints.boundaries, search_direction);
//                //trim_increment_relative(solver_parameters.constraints.relative_boundaries, search_direction);
//
//            }
//
//
//            double search_step = perform_line_search<LineSearch>(
//                solver_parameters.line_search, function, argument, r, search_direction);
//
//            if (!std::isfinite(search_step)) {
//                result->result_code = numerical_result_code_t::LineSearchFailed;
//                break;
//            }
//
//            argument_increment = search_step * search_direction;
//            argument += argument_increment;
//
//            r = function.residuals(argument);
//
//
//            // Проверка критерия выхода по малому относительному приращению
//            double argument_increment_metric = solver_parameters.step_criteria_assuming_search_step
//                ? argument_increment_factor(argument, argument_increment)
//                : argument_increment_factor(argument, search_direction);
//            bool argument_increment_criteria =
//                argument_increment_metric < solver_parameters.argument_increment_norm;
//            bool custom_criteria = false;
//            if (custom_criteria || argument_increment_criteria) {
//                result->result_code = numerical_result_code_t::Converged;
//                break;
//            }
//
//            //double argument_increment_metric = solver_parameters.step_criteria_assuming_search_step
//            //    ? argument_increment_factor(gains, argument_increment)
//            //    : argument_increment_factor(gains, search_direction);
//
//            //result->argument_increment_criteria =
//            //    argument_increment_metric < solver_parameters.argument_increment_norm;
//            //result->custom_stop_criteria = function.custom_stop_criteria(r);
//            //result->enough_iterations_criteria = iteration + 1 == solver_parameters.iteration_count;
//
//            //if (result->has_any_stop_criteria()) {
//            //    result->converged = true;
//            //    break;
//            //}
//        }
//    };
//
//};


TEST(OptimizeGaussNewton, Test1)
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

    VectorXd initial = VectorXd::Zero(2);
    simple_sum_of_squares_function function;

    fixed_solver_parameters_t<-1, 0, golden_section_search> parameters;
    fixed_solver_result_t<-1> result;

    fixed_optimize_gauss_newton::optimize(function, initial, parameters, &result);
    ASSERT_EQ(result.result_code, numerical_result_code_t::Converged);
    ASSERT_NEAR(result.argument(0), 2.0, parameters.argument_increment_norm);
    ASSERT_NEAR(result.argument(1), 1.0, parameters.argument_increment_norm);
}

TEST(OptimizeGaussNewton, RosenbrokFunction)
{
    /// @brief J(x1, x2) = (x1 - 2)^2 + (x2 - 1)^2 
    class rosenbrock_function_t : public fixed_least_squares_function_t
    {
    public:
        VectorXd residuals(const VectorXd& x) {
            VectorXd result(2);
            result[0] = 10 * (x(1) - x(0) * x(0));
            result[1] = 1 - x(0);
            return result;
        }
    };

    VectorXd initial = VectorXd::Zero(2);
    rosenbrock_function_t function;

    fixed_solver_parameters_t<-1, 0, golden_section_search> parameters;
    fixed_solver_result_t<-1> result;

    fixed_optimize_gauss_newton::optimize(function, initial, parameters, &result);
    ASSERT_EQ(result.result_code, numerical_result_code_t::Converged);
    ASSERT_NEAR(result.argument(0), 1.0, parameters.argument_increment_norm);
    ASSERT_NEAR(result.argument(1), 1.0, parameters.argument_increment_norm);
};

template <typename LayerType>
class data_processor_t {
public:
    virtual void process_data(size_t step_index, const LayerType& layer) = 0;
};


class isothermal_qsm_simulation_result_collector_t
    : public data_processor_t<density_viscosity_quasi_layer<true>>
{
public:
    typedef density_viscosity_quasi_layer<true> layer_type;
protected:
    vector<double> pipe_pressure_out;
public:
    isothermal_qsm_simulation_result_collector_t(const vector<double>& times)
        : pipe_pressure_out(times.size(), std::numeric_limits<double>::quiet_NaN())
    {

    }
    virtual void process_data(size_t step_index, 
        const density_viscosity_quasi_layer<true>& layer) override
    {
        // at() - проверяет выход за границы массива
        //pipe_pressure_out.at(step_index) = layer.pressure.back();
        pipe_pressure_out[step_index] = layer.pressure.back();
    }
    const vector<double>& get_pressure_out_calculated() const {
        return pipe_pressure_out;
    }
};

/// @brief Пакетный изотермический квазистатический расчет с предподсчитанным временем
/// делает статический расчет task.solve, а затем столько раз task.step, сколько временных меток в times
/// @tparam Solver МХ или QUICKEST
/// @tparam LayerType Точки под МХ или ячейки под QUICKEST
template <typename Solver, typename LayerType>
inline void isothermal_quasistatic_batch(
    isothermal_quasistatic_PQ_task_t<Solver>& task,
    const vector<double>& times,
    const vector<vector<double>>& boundary_timeseries,
    data_processor_t<LayerType>* data_processor
)
{
    isothermal_quasistatic_PQ_task_boundaries_t initial_boundaries(boundary_timeseries[0]);
    task.solve(initial_boundaries);
    data_processor->process_data(0, task.get_buffer().current());

    for (size_t step_index = 1; step_index < times.size(); step_index++)
    {
        double time_step = times[step_index] - times[step_index - 1];
        isothermal_quasistatic_PQ_task_boundaries_t boundaries(boundary_timeseries[step_index]);

        task.step(time_step, boundaries);

        data_processor->process_data(step_index, task.get_buffer().current());
    }
};

TEST(TimeSeries, PrepareTimeSeries)
{
    // Параметры трубы
    pipe_properties_t pipe;
    // Путь к результатам ресёрча
    string path;
    // Путь к реальным данным с трубопровода
    std::string folder = "../research/2024-08-quasistationary-with-real-data/data/";
    // Временные ряды краевых условий
    vector<pair<vector<time_t>, vector<double>>> control_tag_data;
    // Временные ряды эталонных данных
    vector<pair<vector<time_t>, vector<double>>> etalon_tag_data;
    // Задаём период
    string start_period = "02.08.2021 00:00:00";
    string end_period = "02.08.2021 02:00:00";

    // Указываем имя файла и желаемый шаг новой сетки
    string file_name = folder + "coord_heights.csv";
    //Желаемый шаг
    double desired_dx = 200;

    // Создаём новый профиль с постоянным шагом
    pipe.profile = pipe_profile_uniform::get_uniform_profile_from_csv(desired_dx, file_name);
    pipe.wall.diameter = 1;

    vector<pair<string, string>>parameters =
    {
        { folder + "Q_in", "m3/h-m3/s"s },
        { folder + "p_in", "MPa"s },
        { folder + "rho_in", "kg/m3"s },
        { folder + "visc_in", "mm^2/s-m^2/s"s },
        { folder + "p_out", "MPa"s}

    };

    // Считываем временные ряды параметров
    csv_multiple_tag_reader tags(parameters);
    etalon_tag_data = { control_tag_data.back() };
    control_tag_data.pop_back();

    // Помещаем временные ряды в вектор
    vector_timeseries_t control_parameters_time_series(control_tag_data);
    vector_timeseries_t etalon_parameters_time_series(etalon_tag_data);

    double step = 600;

    time_t start_period_time = max(control_parameters_time_series.get_start_date(), etalon_parameters_time_series.get_start_date());
    time_t end_period_time = min(control_parameters_time_series.get_end_date(), etalon_parameters_time_series.get_end_date());
    time_t duration = (end_period_time - start_period_time);

    size_t dots_count = static_cast<size_t>(ceil(duration / step) + 0.00001);

    vector<double> times(dots_count);
    vector<vector<double>> control_data(dots_count);
    vector<vector<double>> etalon_pressure(dots_count);

    for (size_t i = 0; i < dots_count; i++)
    {
        times[i] = step * i;
        time_t t = start_period_time + static_cast<time_t>(times[i] + 0.5);

        control_data[i] = control_parameters_time_series(t);
        etalon_pressure[i] = control_parameters_time_series(t);

    }
    
};

pipe_properties_t read_profile_data()
{
    pipe_properties_t pipe;
    // Путь к реальным данным с трубопровода
    std::string folder = "../research/2024-08-quasistationary-with-real-data/data/";

    // Указываем имя файла и желаемый шаг новой сетки
    string file_name = folder + "coord_heights.csv";
    //Желаемый шаг
    double desired_dx = 200;

    // Создаём новый профиль с постоянным шагом
    pipe.profile = pipe_profile_uniform::get_uniform_profile_from_csv(desired_dx, file_name);
    pipe.wall.diameter = 1;
    return pipe;
};


inline pipe_properties_t prepare_pipe(const std::string& path)
{
    // Указываем имя файла и желаемый шаг новой сетки
    //string file_name = folder + "coord_heights.csv";
    //Желаемый шаг
    std::string folder = path + "coord_heights.csv";
    double desired_dx = 200;
    pipe_properties_t pipe;

    // Создаём новый профиль с постоянным шагом
    pipe.profile = pipe_profile_uniform::get_uniform_profile_from_csv(desired_dx, folder);
    pipe.wall.diameter = 1;

    return pipe;
};

inline std::tuple<vector<double>, vector<vector<double>>, vector<double>> prepare_real_data(const std::string& path_to_real_data)
{
    // Временные ряды краевых условий
    vector<pair<vector<time_t>, vector<double>>> control_tag_data;
    // Временные ряды эталонных данных
    vector<pair<vector<time_t>, vector<double>>> etalon_tag_data;
    // Задаём период
    string start_period = "01.08.2021 00:00:00";
    string end_period = "01.09.2021 00:00:00";

    vector<pair<string, string>>parameters =
    {
        { path_to_real_data + "Q_in", "m3/h-m3/s"s },
        { path_to_real_data + "p_in", "MPa"s },
        { path_to_real_data + "rho_in", "kg/m3"s },
        { path_to_real_data + "visc_in", "mm^2/s-m^2/s"s },
        { path_to_real_data + "p_out", "MPa"s}

    };

    // Считываем временные ряды параметров
    csv_multiple_tag_reader tags(parameters);
    control_tag_data = tags.read_csvs(start_period, end_period);
    etalon_tag_data = { control_tag_data.back() };
    control_tag_data.pop_back();

    // Помещаем временные ряды в вектор
    vector_timeseries_t control_parameters_time_series(control_tag_data);
    vector_timeseries_t etalon_parameters_time_series(etalon_tag_data);

    double step = 60;

    time_t start_period_time = max(control_parameters_time_series.get_start_date(), etalon_parameters_time_series.get_start_date());
    time_t end_period_time = min(control_parameters_time_series.get_end_date(), etalon_parameters_time_series.get_end_date());
    time_t duration = (end_period_time - start_period_time);

    size_t dots_count = static_cast<size_t>(ceil(duration / step) + 0.00001);

    vector<double>  times = vector<double>(dots_count);
    vector<vector<double>> control_data = vector<vector<double>>(dots_count);
    vector<double> etalon_pressure = vector<double>(dots_count);

    for (size_t i = 0; i < dots_count; i++)
    {
        times[i] = step * i;
        time_t t = start_period_time + static_cast<time_t>(times[i] + 0.5);

        control_data[i] = control_parameters_time_series(t);
        etalon_pressure[i] = etalon_parameters_time_series(t).front();

    };

    return std::make_tuple(std::move(times), std::move(control_data), std::move(etalon_pressure));
}

//void print_diff_pressure_before_after(const double d_before, const double d_after)
//{
//    python_printer printer;
//
//    printer.print_profiles<double>(static_cast<time_t>(0),
//        times,
//        vector<vector<double>>{ calc_vector_residuals(d_before), calc_vector_residuals(d_after) },
//        "time,time,diff_press_before,diff_press_after",
//        folder + "diff_press.csv");
//}


struct ident_parameters_isothermal_qsm {
    double diameter_adaptation{ 1 };
    double friction_adaptation{ 1 };
    void set_adaptation(const pipe_properties_t& nominal_pipe, pipe_properties_t* pipe_to_ident) const 
    {
        pipe_to_ident->wall.diameter = diameter_adaptation * nominal_pipe.wall.diameter;
        pipe_to_ident->wall.resistance_function_adaptation = friction_adaptation;
    }
    
};

struct ident_isothermal_qsm_pipe_settings {
    bool ident_diameter{false};
    bool ident_friction{false};
    ident_parameters_isothermal_qsm unpack_ident_parameters(const VectorXd& packed_ident_parameters) {
        ident_parameters_isothermal_qsm result;
        size_t index = 0;
        if (ident_diameter) {
            result.diameter_adaptation = packed_ident_parameters(index++);
        }
        if (ident_friction) {
            result.friction_adaptation = packed_ident_parameters(index++);
        }
        return result;
    }
    void check_parameters() const {
        // если ничего не задано или наоборот, все задано, кинуть исключение
        // throw std::runtime_error("");
    }
};


class ident_isothermal_qsm_pipe_diameter_t : public fixed_least_squares_function_t
{
    const ident_isothermal_qsm_pipe_settings settings;
    const pipe_properties_t pipe_nominal;
    pipe_properties_t pipe_to_ident;

    const vector<double>& times;
    const vector<vector<double>>& control_data;
    const vector<double>& etalon_pressure;
public:
    ident_isothermal_qsm_pipe_diameter_t(
        const ident_isothermal_qsm_pipe_settings& settings,
        const pipe_properties_t& pipe, const vector<double>& times, const vector<vector<double>>& control_data, const vector<double>& etalon_pressure)
        : settings(settings)
        , pipe_nominal{ pipe }
        , pipe_to_ident{ pipe }
        , times{ times }
        , control_data{ control_data }
        , etalon_pressure{ etalon_pressure }
    {
        settings.check_parameters();
    }
protected:
    vector<double> calc_vector_residuals(const ident_parameters_isothermal_qsm& ident_parameters)
    {
        isothermal_qsm_simulation_result_collector_t collector(times);
        ident_parameters.set_adaptation(pipe_nominal, &pipe_to_ident);

        isothermal_quasistatic_PQ_task_t<quickest_ultimate_fv_solver> task(pipe_to_ident);
        isothermal_quasistatic_batch<quickest_ultimate_fv_solver, isothermal_qsm_simulation_result_collector_t::layer_type>(
            task,
            times,
            control_data,
            &collector
        );

        const vector<double>& calc_pressure = collector.get_pressure_out_calculated();
        vector<double> simulation_result(times.size());

        std::transform(calc_pressure.begin(), calc_pressure.end(), etalon_pressure.begin(), simulation_result.begin(),
            [](double etalon, double calc) { return etalon - calc; });


        return simulation_result;
    }
public:
    virtual VectorXd residuals(const VectorXd& d) override {

        ident_parameters_isothermal_qsm ident_parameters;
        ident_parameters.diameter_adaptation = d(0);


        vector<double> simulation_result = calc_vector_residuals(ident_parameters);

        Eigen::Map<VectorXd> result(simulation_result.data(), simulation_result.size());

        return result;
    }
    double ident(fixed_solver_result_t<-1>* result = nullptr, fixed_solver_result_analysis_t<-1>* analysis = nullptr) {
        fixed_solver_parameters_t<-1, 0, golden_section_search> parameters;

        fixed_solver_result_t<-1> result_local;
        if (result == nullptr) {
            result = &result_local;
        }

        VectorXd initial_d(1); 
        initial_d(0) = 1;
        fixed_optimize_gauss_newton::optimize(*this, initial_d, parameters, result, analysis);

        if (result->result_code == numerical_result_code_t::Converged) {
            return result->argument(0);
        }
        else {
            return std::numeric_limits<double>::quiet_NaN();
        }
    }
};

void print_j_d(const double d, const double J, const std::string folder)
{
    python_printer printer;
    printer.print_profiles<double>(static_cast<time_t>(0),
        { d },
        vector<vector<double>>{ {J} },
        "time,d,J",
        folder + "j.csv");
}


TEST(OptimiseGaussNewton, PipeIdentification)
{
    // Путь к реальным данным с трубопровода
    std::string data_path = "../research/2024-08-quasistationary-with-real-data/data/";

    pipe_properties_t pipe = prepare_pipe(data_path);
    auto [times, control_data, etalon_pressure] = prepare_real_data(data_path);

    ident_isothermal_qsm_pipe_settings ident_settings;
    ident_settings.ident_diameter = true;

    ident_isothermal_qsm_pipe_diameter_t test_ident(ident_settings, pipe, times, control_data, etalon_pressure);

    fixed_solver_result_t<-1> result;
    fixed_solver_result_analysis_t<-1> analysis;

    double result_d = test_ident.ident(&result, &analysis);
}

TEST(OptimiseGaussNewton, PipeIdentificationWithPrinter)
{
    std::string path = prepare_research_folder_for_qsm_model();


    //string folder = prepare_research_folder_for_qsm_model();

    //VectorXd initial_d = VectorXd::Zero(1);
    //initial_d(0) = 1;

    //// Путь к реальным данным с трубопровода
    //std::string data_path = "../research/2024-08-quasistationary-with-real-data/data/";
    //pipe_properties_t pipe = prepare_pipe(data_path);
    //auto [times, control_data, etalon_pressure] = prepare_real_data(data_path);

    //ident_isothermal_qsm_pipe_diameter_t test_ident;
    ////VectorXd residuals = test_ident.residuals(initial_d);



    //test_ident.print_diff_pressure_before_after(initial_d(0), result.argument(0));

}