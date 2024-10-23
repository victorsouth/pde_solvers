#pragma once

namespace pde_solvers {
;

/// @brief Уравнение адвекции постоянного объемного расхода Q по длине трубы
/// Реализация без привязки к трубе
class pipe_advection_pde_t : public pde_t<1>
{
public:
    using pde_t<1>::equation_coeffs_type;
    using pde_t<1>::right_party_type;
    using pde_t<1>::var_type;
protected:
    /// @brief Объемный расход, полагается постоянным по длине трубы
    const double volumetric_flow;
    /// @brief Площадь сечения трубы
    const double pipe_area;
    /// @brief Расчетная сетка
    const vector<double>& pipe_profile_coordinates;
public:
    pipe_advection_pde_t(double pipe_area, double volumetric_flow,
        const vector<double>& pipe_profile_coordinates)
        : volumetric_flow(volumetric_flow)
        , pipe_area(pipe_area)
        , pipe_profile_coordinates(pipe_profile_coordinates)
    {}

    /// @brief Возвращает известную уравнению сетку
    virtual const vector<double>& get_grid() const override {
        return pipe_profile_coordinates;
    }

    /// @brief Левая часть уравнения давекции
    /// @param grid_index Индекс узла расчетной сетки
    /// @param point_vector Вектор значений переменных
    /// @return Коэффициенты уравнения
    virtual equation_coeffs_type getEquationsCoeffs(
        size_t grid_index, const var_type& point_vector) const override
    {
        double v = volumetric_flow / pipe_area;
        return v;
    }

    virtual equation_coeffs_type getEquationsCoeffsInv(
        size_t grid_index, const var_type& point_vector) const override
    {
        double v = volumetric_flow / pipe_area;
        return 1 / v;
    }

    /// @brief Получение собственных чисел и соответствующих им ЛЕВЫХ собственных векторов
    /// \return Список собственных чисел, список собственных векторов
    virtual pair<var_type, equation_coeffs_type> GetLeftEigens(
        size_t grid_index, const var_type& point_vector) const override
    {
        double v = getEquationsCoeffs(grid_index, point_vector);
        return std::make_pair(v, v);
    }

    virtual pair<var_type, equation_coeffs_type> GetRightEigens(
        size_t grid_index, const var_type& point_vector) const override
    {
        double v = getEquationsCoeffs(grid_index, point_vector);
        return std::make_pair(v, v);
    }

    /// @brief Правая часть уравнения адвекция - нулевая
    virtual right_party_type getSourceTerm(
        size_t grid_index, const var_type& point_var) const override
    {
        return 0.0;
    }
    virtual var_type GetRightEigenVector(
        size_t profile_index, size_t eigen_index, const var_type& u) const override
    {
        throw std::logic_error("not implemented");
    }

    virtual double get_wave_strength(
        size_t profile_index, size_t eigen_index, const var_type& u) const override
    {
        throw std::logic_error("not implemented");
    }
};

/// @brief Уравление адвекции (транспортное уравнение) на основе объемного расхода
class PipeQAdvection : public pde_t<1>
{
public:
    using pde_t<1>::equation_coeffs_type;
    using pde_t<1>::right_party_type;
    using pde_t<1>::var_type;
protected:
    /// @brief Труба
    const pipe_properties_t& pipe;
    /// @brief Объемный расход
    const vector<double>& Q;
public:
    PipeQAdvection(const pipe_properties_t& pipe,
        const vector<double>& vol_flow)
        : pipe(pipe)
        , Q(vol_flow)
    {}

    const pipe_properties_t& get_pipe() const
    {
        return pipe;
    }

    /// @brief Возвращает известную уравнению сетку
    virtual const vector<double>& get_grid() const override {
        return pipe.profile.coordinates;
    }

    /// @brief Левая часть
    /// @param index 
    /// @return 
    virtual equation_coeffs_type getEquationsCoeffs(
        size_t grid_index, const var_type& point_vector) const override
    {
        double S_0 = pipe.wall.getArea();
        double v = Q[grid_index] / S_0;
        return v;
    }

    virtual equation_coeffs_type getEquationsCoeffsInv(
        size_t grid_index, const var_type& point_vector) const override
    {
        double S_0 = pipe.wall.getArea();
        double v = Q[grid_index] / S_0;
        return 1 / v;
    }

    /// @brief Получение собственных чисел и соответствующих им ЛЕВЫХ собственных векторов
    /// \return Список собственных чисел, список собственных векторов
    virtual pair<var_type, equation_coeffs_type> GetLeftEigens(
        size_t grid_index, const var_type& point_vector) const override
    {
        double v = getEquationsCoeffs(grid_index, point_vector);
        return std::make_pair(v, v);
    }

    virtual pair<var_type, equation_coeffs_type> GetRightEigens(
        size_t grid_index, const var_type& point_vector) const override
    {
        double v = getEquationsCoeffs(grid_index, point_vector);
        return std::make_pair(v, v);
    }

    /// @brief Правая часть уравнения адвекция - нулевая
    virtual right_party_type getSourceTerm(
        size_t grid_index, const var_type& point_var) const override
    {
        return 0.0;
    }
    virtual var_type GetRightEigenVector(
        size_t profile_index, size_t eigen_index, const var_type& u) const override
    {
        throw std::logic_error("not implemented");
    }

    virtual double get_wave_strength(
        size_t profile_index, size_t eigen_index, const var_type& u) const override
    {
        throw std::logic_error("not implemented");
    }
};

}