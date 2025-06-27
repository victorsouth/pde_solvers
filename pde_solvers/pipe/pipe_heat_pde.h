#pragma once

namespace pde_solvers {
;


/// @brief Уравнение притока тепла через формулу Ньютона q = -Kt(T - T_outer)
/// Беспартийный расчет
/// по документу "Полный вывод НЕизотермических уравнений..."
class PipeHeatInflowConstArea : public pde_t<1> {
public:
    using pde_t<1>::equation_coeffs_type;
    using pde_t<1>::right_party_type;
    using pde_t<1>::var_type;
protected:
    const pipe_noniso_properties_t& pipe;
    const oil_parameters_t& oil;

    /// @brief массовый расход (в ячейках или точках??)
    const std::vector<double>& mass_flow;

public:

    PipeHeatInflowConstArea(const pipe_noniso_properties_t& pipe, const oil_parameters_t& oil,
        const std::vector<double>& mass_flow)
        : pipe(pipe)
        , oil(oil)
        , mass_flow(mass_flow)
    {}

    /// @brief Возвращает известную уравнению сетку
    virtual const std::vector<double>& get_grid() const override {
        return pipe.profile.coordinates;
    }

    /// @brief Левая часть
    /// @param index 
    /// @return 
    virtual equation_coeffs_type getEquationsCoeffs(
        size_t grid_index, const var_type& point_vector) const override
    {
        double S_0 = pipe.wall.getArea();
        double density = oil.density.nominal_density;
        double v = mass_flow[grid_index] / (S_0 * density);
        return v;
    }

    virtual equation_coeffs_type getEquationsCoeffsInv(
        size_t grid_index, const var_type& point_vector) const override
    {
        return 1 / getEquationsCoeffs(grid_index, point_vector);
    }

    /// @brief Получение собственных чисел и соответствующих им ЛЕВЫХ собственных векторов
    /// \return Список собственных чисел, список собственных векторов
    virtual std::pair<var_type, equation_coeffs_type> GetLeftEigens(
        size_t grid_index, const var_type& point_vector) const override
    {
        double v = getEquationsCoeffs(grid_index, point_vector);
        return std::make_pair(v, v);
    }

    virtual std::pair<var_type, equation_coeffs_type> GetRightEigens(
        size_t grid_index, const var_type& point_vector) const override
    {
        double v = getEquationsCoeffs(grid_index, point_vector);
        return std::make_pair(v, v);
    }

    /// @brief Правая часть уравнения притока тепла
    /// по документу "Полный вывод НЕизотермических уравнений...", формула (8)
    virtual right_party_type getSourceTerm(
        size_t grid_index, const var_type& point_var) const override
    {
        double T = point_var;

        double d = pipe.wall.diameter;
        double S_0 = pipe.wall.getArea();
        double density = oil.density.nominal_density;
        double v = mass_flow[grid_index] / (S_0 * density);

        double Re = v * pipe.wall.diameter / oil.viscosity(T);
        double lambda = hydraulic_resistance_shifrinson(Re, pipe.wall.relativeRoughness());

        double n_inner = -lambda / 2 * v * v * v / d;
        double q = pipe.heat.get_heat_flow_newton(T);

        double Cp = oil.heat.HeatCapacity;

        double s = (M_PI * d * q - density * S_0 * n_inner)
            / (density * S_0 * Cp);

        return s;
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



/// @brief Уравление притока тепла, учитывающее теплообмен с динамическим грунтом
/// Беспартийный расчет
/// по документу "Полный вывод НЕизотермических уравнений..."
class PipeSoilHeatInflowConstArea : public pde_t<1>
{
public:
    using pde_t<1>::equation_coeffs_type;
    using pde_t<1>::right_party_type;
    using pde_t<1>::var_type;
protected:
    /// @brief Труба
    const zoned_pipe_properties& pipe;
    /// @brief Нефть (считается одной и той же по всей трубе)
    const oil_parameters_t& oil;
    /// @brief Массовый расход (в ячейках или точках??)
    const std::vector<double>& mass_flow;
    /// @brief Средняя температура грунта (если null, то используется ambient_temperature)
    const std::vector<double>* temperature_soil{ nullptr };
public:
    /// @brief Конструктор для модель с тепловой динамикой грунта
    PipeSoilHeatInflowConstArea(const zoned_pipe_properties& pipe, const oil_parameters_t& oil,
        const std::vector<double>& mass_flow, const std::vector<double>& temperature_soil)
        : pipe(pipe)
        , oil(oil)
        , mass_flow(mass_flow)
        , temperature_soil(&temperature_soil)
    {}
    /// @brief Конструктор для модели Шухова
    PipeSoilHeatInflowConstArea(const zoned_pipe_properties& pipe, const oil_parameters_t& oil,
        const std::vector<double>& mass_flow)
        : pipe(pipe)
        , oil(oil)
        , mass_flow(mass_flow)
        , temperature_soil(nullptr)
    {}

    /// @brief Возвращает известную уравнению сетку
    virtual const std::vector<double>& get_grid() const override {
        return pipe.profile.coordinates;
    }

    /// @brief Левая часть
    /// @param index 
    /// @return 
    virtual equation_coeffs_type getEquationsCoeffs(
        size_t grid_index, const var_type& point_vector) const override
    {
        double S_0 = pipe.wall.getArea();
        double density = oil.density.nominal_density;
        double v = mass_flow[grid_index] / (S_0 * density);
        return v;
    }

    virtual equation_coeffs_type getEquationsCoeffsInv(
        size_t grid_index, const var_type& point_vector) const override
    {
        return 1 / getEquationsCoeffs(grid_index, point_vector);
    }

    /// @brief Получение собственных чисел и соответствующих им ЛЕВЫХ собственных векторов
    /// \return Список собственных чисел, список собственных векторов
    virtual std::pair<var_type, equation_coeffs_type> GetLeftEigens(
        size_t grid_index, const var_type& point_vector) const override
    {
        double v = getEquationsCoeffs(grid_index, point_vector);
        return std::make_pair(v, v);
    }

    virtual std::pair<var_type, equation_coeffs_type> GetRightEigens(
        size_t grid_index, const var_type& point_vector) const override
    {
        double v = getEquationsCoeffs(grid_index, point_vector);
        return std::make_pair(v, v);
    }




    /// @brief Правая часть уравнения притока тепла
    /// по документу "Полный вывод НЕизотермических уравнений...", формула (8)
    virtual right_party_type getSourceTerm(
        size_t grid_index, const var_type& point_var) const override
    {
        double T = point_var;

        double d = pipe.wall.diameter;
        double S_0 = pipe.wall.getArea();
        double density = oil.density.nominal_density;
        double v = mass_flow[grid_index] / (S_0 * density);

        double Re = v * pipe.wall.diameter / oil.viscosity(T);
        double lambda = pipe.resistance_function(Re);

        double n_inner = -lambda / 2 * v * v * v / d;

        //const auto& zone = pipe.get_heat_zone(grid_index);
        const auto& model = pipe.get_heat_zone_coeff(grid_index);

        double q;
        if (temperature_soil != nullptr) {
            double Tsoil = (*temperature_soil)[grid_index];
            q = model.get_heat_flow(T, Tsoil, model.temperature_ambient);
        }
        else {
            q = model.get_heat_flow(T, model.temperature_ambient);
        }

        double Cp = oil.heat.HeatCapacity;

        double s = (M_PI * d * q - density * S_0 * n_inner)
            / (density * S_0 * Cp);

        return s;
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


/// @brief Уравление притока тепла, учитывающее теплообмен с динамическим грунтом,
/// партийность жидкости
/// по документу "Полный вывод НЕизотермических уравнений..."
class PipeSoilHeatInflowSortedConstArea : public pde_t<1>
{
public:
    using pde_t<1>::equation_coeffs_type;
    using pde_t<1>::right_party_type;
    using pde_t<1>::var_type;
protected:
    /// @brief Профиль свойств жидкости
    const fluid_properties_profile_t& oil;
    /// @brief Труба
    const zoned_pipe_properties& pipe;
    /// @brief Объемный расход
    const std::vector<double>& Q;
    /// @brief Средняя температура грунта (если null, то используется ambient_temperature)
    const std::vector<double>* temperature_soil{ nullptr };
public:
    /// @brief Конструктор для партийной модели с тепловой динамикой грунта
    PipeSoilHeatInflowSortedConstArea(const zoned_pipe_properties& pipe,
        const fluid_properties_profile_t& oil,
        const std::vector<double>& vol_flow, const std::vector<double>& temperature_soil)
        : pipe(pipe)
        , oil(oil)
        , Q(vol_flow)
        , temperature_soil(&temperature_soil)
    {}
    /// @brief Конструктор для модели Шухова
    PipeSoilHeatInflowSortedConstArea(const zoned_pipe_properties& pipe,
        const fluid_properties_profile_t& oil,
        const std::vector<double>& vol_flow)
        : pipe(pipe)
        , oil(oil)
        , Q(vol_flow)
        , temperature_soil(nullptr)
    {}

    /// @brief Возвращает известную уравнению сетку
    virtual const std::vector<double>& get_grid() const override {
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
    virtual std::pair<var_type, equation_coeffs_type> GetLeftEigens(
        size_t grid_index, const var_type& point_vector) const override
    {
        double v = getEquationsCoeffs(grid_index, point_vector);
        return std::make_pair(v, v);
    }

    virtual std::pair<var_type, equation_coeffs_type> GetRightEigens(
        size_t grid_index, const var_type& point_vector) const override
    {
        double v = getEquationsCoeffs(grid_index, point_vector);
        return std::make_pair(v, v);
    }

    /// @brief Правая часть уравнения притока тепла
    /// по документу "Полный вывод НЕизотермических уравнений...", формула (8)
    virtual right_party_type getSourceTerm(
        size_t grid_index, const var_type& point_var) const override
    {
        double T = point_var;

        double da = pipe.adaptation.diameter;
        double d = pipe.wall.diameter;
        double S_0 = pipe.wall.getArea();
        d *= da;
        S_0 *= da * da; // квадрат
        double density = oil.nominal_density[grid_index];
        double v = Q[grid_index] / S_0;

        double Re = v * d / oil.get_viscosity(grid_index, T);
        double lambda = pipe.resistance_function(Re);
        lambda *= pipe.adaptation.friction;

        double n_inner = -lambda / 2 * v * v * v / d;

        //auto& zone = pipe.get_heat_zone(grid_index);
        auto& model = pipe.get_heat_zone_coeff(grid_index);

        double q;
        if (temperature_soil != nullptr) {
            double Tsoil = (*temperature_soil)[grid_index];
            q = model.get_heat_flow(T, Tsoil, model.temperature_ambient);
        }
        else {
            q = model.get_heat_flow(T, model.temperature_ambient);
        }

        double Cp = oil.heat_capacity; // Можно использовать формулу Крэго.
        Cp *= pipe.adaptation.heat_capacity;

        double s_heat = M_PI * d * q;
        double s_friction = -density * S_0 * n_inner;

        double s = (M_PI * d * q - density * S_0 * n_inner)
            / (density * S_0 * Cp);

        return s;
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