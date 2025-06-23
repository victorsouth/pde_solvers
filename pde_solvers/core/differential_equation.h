#pragma once

namespace pde_solvers {

/// @brief Дифференциальное уравнение
template <size_t Dimension>
class differential_equation_t
{
public:
    /// @brief Возвращает известную уравнению сетку
    virtual const std::vector<double>& get_grid() const = 0;
};




/// @brief Обыкновенное дифференциальное уравнение - базовый класс
template <size_t Dimension>
class ode_t : public differential_equation_t<Dimension> {
public:
    /// @brief Тип переменной уравнения
    typedef typename fixed_system_types<Dimension>::var_type var_type;
    /// @brief Тип правой части уравнения
    typedef typename fixed_system_types<Dimension>::right_party_type right_party_type;
    /// @brief Тип матрицы коэффициентов
    typedef typename fixed_system_types<Dimension>::equation_coeffs_type equation_coeffs_type;
public:
    /// @brief Возвращает правую часть системы уравнений dx/dt = f(x)
    /// @param grid_index 
    /// @param point_vector 
    /// @return Значение правой части ОДУ
    virtual right_party_type ode_right_party(
        size_t grid_index, const var_type& point_vector) const = 0;
};



/// @brief Уравнение в частных производных
template <size_t Dimension>
class pde_t : public ode_t<Dimension>
{
public:
    /// @brief Тип матрицы коэффициентов
    using typename ode_t<Dimension>::equation_coeffs_type;
    /// @brief Тип правой части уравнения
    using typename ode_t<Dimension>::right_party_type;
    /// @brief Тип переменной уравнения
    using typename ode_t<Dimension>::var_type;
public:
    /// @brief Получение матрицы коэффициентов уравнений
    /// @param grid_index Индекс точки в расчетной сетке
    /// @param point_vector Вектор переменных в точке
    /// @return Матрица коэффициентов уравнений
    virtual equation_coeffs_type getEquationsCoeffs(
        size_t grid_index, const var_type& point_vector) const = 0;
    /// @brief Получение обратной матрицы коэффициентов уравнений
    /// @param grid_index Индекс точки в расчетной сетке
    /// @param point_vector Вектор переменных в точке
    /// @return Обратная матрица коэффициентов уравнений
    virtual equation_coeffs_type getEquationsCoeffsInv(
        size_t grid_index, const var_type& point_vector) const = 0;
    /// @brief Получение вектора правой части системы (источниковый член)
    /// @param grid_index Индекс точки в расчетной сетке
    /// @param point_vector Вектор переменных в точке
    /// @return Вектор правой части системы
    virtual right_party_type getSourceTerm(
        size_t grid_index, const var_type& point_vector) const = 0;
    /// @brief Получение собственных чисел и соответствующих им ЛЕВЫХ собственных векторов
    /// @param index Индекс точки в расчетной сетке
    /// @param u Вектор переменных в точке
    /// @return Пара: Список собственных чисел, список собственных векторов
    virtual std::pair<var_type, equation_coeffs_type> GetLeftEigens(
        size_t index, const var_type& u) const = 0;
    /// @brief Вычисление правых собственных векторов и собственных значений
    /// @return Пара: вектор собственных значений и матрица правых собственных векторов
    virtual std::pair<var_type, equation_coeffs_type> GetRightEigens(
        size_t index, const var_type& u) const = 0;
    /// @brief Получение правого собственного вектора
    /// @param profile_index Индекс точки в расчетной сетке
    /// @param eigen_index Индекс собственного значения
    /// @param u Вектор переменных в точке
    virtual var_type GetRightEigenVector(
        size_t profile_index, size_t eigen_index, const var_type& u) const = 0;
    /// @brief Где-то используется, но в pde_solvers нет (9.05.2025)
    virtual double get_wave_strength(
        size_t profile_index, size_t eigen_index, const var_type& u) const = 0;
    /// @brief Вычисление правой части системы ОДУ
    /// Реализация метода из базового класса ode_t. Вычисляет правую часть как:
    /// Ux = A^-1 * b, где A - матрица коэффициентов, b - вектор источников
    virtual right_party_type ode_right_party(
        size_t grid_index, const var_type& point_vector) const override
    {
        // A * Ux = b => Ux = A^-1 * b

        equation_coeffs_type Ainv = getEquationsCoeffsInv(grid_index, point_vector);
        right_party_type b = getSourceTerm(grid_index, point_vector);
        right_party_type result = Ainv * b;
        return result;
    }
};

}