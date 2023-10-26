#pragma once

/// @brief Уравление адвекции (транспортное уравнение) на основе объемного расхода
class PipeQAdvection : public pde_t<1>
{
public:
    using pde_t<1>::equation_coeffs_type;
    using pde_t<1>::right_party_type;
    using pde_t<1>::var_type;
protected:
    /// @brief Труба
    const PipeProperties& pipe;
    /// @brief Объемный расход
    const vector<double>& Q;
public:
    PipeQAdvection(const PipeProperties& pipe,
        const vector<double>& vol_flow)
        : pipe(pipe)
        , Q(vol_flow)
    {}

    const PipeProperties& get_pipe() const
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
