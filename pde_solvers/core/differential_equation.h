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
    typedef typename fixed_system_types<Dimension>::var_type var_type;
    typedef typename fixed_system_types<Dimension>::right_party_type right_party_type;
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
    using typename ode_t<Dimension>::equation_coeffs_type;
    using typename ode_t<Dimension>::right_party_type;
    using typename ode_t<Dimension>::var_type;
public:

    virtual equation_coeffs_type getEquationsCoeffs(
        size_t grid_index, const var_type& point_vector) const = 0;

    virtual equation_coeffs_type getEquationsCoeffsInv(
        size_t grid_index, const var_type& point_vector) const = 0;
    virtual right_party_type getSourceTerm(
        size_t grid_index, const var_type& point_vector) const = 0;


    /// @brief Получение собственных чисел и соответствующих им ЛЕВЫХ собственных векторов
    /// \param curr
    /// \param index
    /// \return Список собственных чисел, список собственных векторов
    virtual std::pair<var_type, equation_coeffs_type> GetLeftEigens(
        size_t index, const var_type& u) const = 0;

    virtual std::pair<var_type, equation_coeffs_type> GetRightEigens(
        size_t index, const var_type& u) const = 0;

    virtual var_type GetRightEigenVector(
        size_t profile_index, size_t eigen_index, const var_type& u) const = 0;

    virtual double get_wave_strength(
        size_t profile_index, size_t eigen_index, const var_type& u) const = 0;

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