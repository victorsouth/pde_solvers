#pragma once

namespace pde_solvers {

/// @brief Уравнение трубы для давления и расхода
/// "Схема Годунова для уравнения с постоянным сечением.docx"
class PipeModelPGConstArea : public pde_t<2>
{
public:
    using pde_t<2>::equation_coeffs_type;
    using pde_t<2>::right_party_type;
    using pde_t<2>::var_type;
protected:
    pipe_properties_t pipe;
    oil_parameters_t oil;

public:
    PipeModelPGConstArea(const pipe_properties_t& pipe, const oil_parameters_t& oil)
        : pipe(pipe)
        , oil(oil)
    {

    }

    /// @brief Возвращает известную уравнению сетку
    virtual const std::vector<double>& get_grid() const override {
        return pipe.profile.coordinates;
    }

    /// @brief Получение матрицы коэффициентов системы уравнений (Row-major)
    /// \param curr
    /// \param index
    /// \return Row-major матрица, массив вектор-строк
    virtual equation_coeffs_type getEquationsCoeffs(
        size_t grid_index, const var_type& point_vector) const override
    {
        double S_0 = pipe.wall.getArea();
        double c = pipe.getSoundVelocity(oil);

        equation_coeffs_type A; // Row-major матрица, массив вектор-строк
        A[0] = { 0, pow(c, 2) / S_0 };
        A[1] = { S_0, 0 };
        return A;
    }
    /// @brief Обратная матрица коэффициентов системы уравнений
    /// @param curr 
    /// @param index 
    /// @return 
    virtual equation_coeffs_type getEquationsCoeffsInv(
        size_t grid_index, const var_type& point_vector) const override
    {
        double S_0 = pipe.wall.getArea();
        double c = pipe.getSoundVelocity(oil);

        std::array<std::array<double, 2>, 2> Ainv;
        Ainv[0] = { 0, 1 / S_0 };
        Ainv[1] = { S_0 / pow(c, 2), 0 };

        return Ainv;
    }

    /// @brief Получение вектора правой части системы уравнений
    /// \param curr
    /// \param index
    /// \return
    virtual var_type getSourceTerm(size_t grid_index, const var_type& point_vector) const override
    {
        double p = point_vector[0];
        double G = point_vector[1];
        double rho = oil.density();
        double S_0 = pipe.wall.getArea();
        double v = G / (rho * S_0);
        double Re = v * pipe.wall.diameter / oil.viscosity();
        double lambda = pipe.resistance_function(Re);
        double tau_w = lambda / 8 * rho * v * abs(v);
        double s1 = -M_PI * pipe.wall.diameter * tau_w;

        var_type s = { 0, s1 };
        return s;
    }

    /// @brief Получение собственных чисел и соответствующих им собственных векторов
    /// \param curr
    /// \param index
    /// \return Список собственных чисел, список собственных векторов
    virtual std::pair<var_type, equation_coeffs_type> GetLeftEigens(
        size_t profile_index, const var_type& u) const override
    {
        std::pair<var_type, equation_coeffs_type> result;

        auto& values = result.first;
        auto& vectors = result.second;

        double pressure = u[0];

        /// По файлу "2021-04-12 Характеристическая форма.xmcd"
        double S_0 = pipe.wall.getArea();
        double c = pipe.getSoundVelocity(oil);

        values = {
            -c,
            c,
        };

        vectors[0] = {
            -S_0 / c,
            1
        };

        vectors[1] = {
            S_0 / c,
            1
        };

        return result;
    }

    /// @brief Получение собственных чисел и соответствующих им собственных векторов
    /// \param curr
    /// \param index
    /// \return Список собственных чисел, список собственных векторов
    virtual std::pair<var_type, equation_coeffs_type> GetRightEigens(
        size_t index, const var_type& u) const override
    {
        std::pair<var_type, equation_coeffs_type> result;

        auto& values = result.first;
        auto& vectors = result.second;

        double pressure = u[0];

        /// По файлу "2021-04-12 Характеристическая форма.xmcd"
        double S_0 = pipe.wall.getArea();
        double c = pipe.getSoundVelocity(oil);

        values = {
            -c,
            c ,
        };

        vectors[0] = {
            -c / S_0,
            1
        };
        vectors[1] = {
            c / S_0,
            1
        };


        return result;
    }

    virtual double get_wave_strength(
        size_t profile_index, size_t eigen_index, const var_type& u) const override
    {
        double p = u[0];
        double G = u[1];

        double S_0 = pipe.wall.getArea();
        double c = pipe.getSoundVelocity(oil);
        if (eigen_index == 0) {
            return 0.5 * G - 0.5 * p * S_0 / c;
        }
        else {
            return 0.5 * G + 0.5 * p * S_0 / c;
        }

    }

    virtual var_type GetRightEigenVector(
        size_t profile_index, size_t eigen_index, const var_type& u) const override
    {
        double pressure = u[0];

        /// По файлу "2021-04-12 Характеристическая форма.xmcd"
        double S_0 = pipe.wall.getArea();
        double c = pipe.getSoundVelocity(oil);

        if (eigen_index == 0) {
            return { -c / S_0, 1 };
        }
        else {
            return { c / S_0, 1 };
        }
    }

    /// @brief Расчет потоков
    /// \param curr
    /// \param index
    /// \return
    var_type getFlux(size_t index, const var_type& u) const
    {
        double P = u[0];
        double G = u[1];
        double c = pipe.getSoundVelocity(oil);
        double S_0 = pipe.wall.getArea();

        var_type F = {
            c * c * G / S_0, // поток для переменной давления
            S_0 * P // поток импульса
        };

        return F;
    }

    /// @brief Расчет задачи Римана на границе
    /// "Схема Годунова для уравнения с постоянным сечением.docx"
    /// @param point_index Индекс границы
    /// @param U_b Значение параметра на граничной ячейке на старом слое
    /// @param boundary_eq Уравнение на границе
    /// @return Решение задачи Римана на границе
    std::array<double, 2> riemann_problem_boundary(
        size_t point_index,
        const std::array<double, 2>& U_b,
        const std::pair<std::array<double, 2>, double>& boundary_eq) const
    {
        std::array<std::array<double, 2>, 2> A;
        std::array<double, 2> b;

        bool is_left_boundary = point_index == 0;

        if (is_left_boundary) {
            // Из уравнения обобщенных инвариантов Римана остается уравнение 
            // через характеристику с _правым_ уклоном (см. документ, там есть рисунок)
            // lambda > 0, K_R
            // Выбор уравнения с правым уклоном за счет хардкода eigen_index = 1
            std::array<double, 2> K_R = GetRightEigenVector(point_index, 1, U_b);
            A[0] = { 1 / K_R[0], -1 / K_R[1] };
            b[0] = U_b[0] / K_R[0] - U_b[1] / K_R[1];

            A[1] = boundary_eq.first;
            b[1] = boundary_eq.second;
            std::array<double, 2> U = solve_linear_system(A, b);
            return U;
        }
        else {
            // Из уравнения обобщенных инвариантов Римана остается уравнение 
            // через характеристику с _левым_ уклоном (см. документ, там есть рисунок)
            // lambda < 0, K_L
            // Выбор уравнения с правым уклоном за счет хардкода eigen_index = 1
            std::array<double, 2> K_L = GetRightEigenVector(point_index, 0, U_b);
            A[0] = { 1 / K_L[0], -1 / K_L[1] };
            b[0] = U_b[0] / K_L[0] - U_b[1] / K_L[1];

            A[1] = boundary_eq.first;
            b[1] = boundary_eq.second;

            std::array<double, 2> U = solve_linear_system(A, b);
            return U;
        }


    }



    /// @brief Расчет задачи Римана для внутренних границ
    /// @param point_index Индекс границы
    /// @param U_L Значение вектора параметров в ячейке слева от границы
    /// @param U_R Значение вектора параметров в ячейке справа от границы
    /// @return Решение задачи Римана на внутренней границе
    std::array<double, 2> riemann_problem_inner(
        size_t point_index,
        const std::array<double, 2>& U_L,
        const std::array<double, 2>& U_R) const
    {
        std::array<std::array<double, 2>, 2> A;
        std::array<double, 2> b;

        // lambda < 0, K_L
        std::array<double, 2> K_L = GetRightEigenVector(point_index, 0, U_L);
        A[0] = { 1 / K_L[0], -1 / K_L[1] };
        b[0] = U_L[0] / K_L[0] - U_L[1] / K_L[1];

        // lambda > 0, K_R
        std::array<double, 2> K_R = GetRightEigenVector(point_index, 1, U_R);
        A[1] = { 1 / K_R[0], -1 / K_R[1] };
        b[1] = U_R[0] / K_R[0] - U_R[1] / K_R[1];


        std::array<double, 2> U = solve_linear_system(A, b);

        return U;
    }

    static std::pair<std::array<double, 2>, double> const_pressure_equation(double pressure) {
        return { {1, 0}, pressure };
    }
    static std::pair<std::array<double, 2>, double> const_mass_flow_equation(double mass_flow) {
        return { {0, 1}, mass_flow };
    }

};

/// @brief Модель трубы с постоянным сечением, для которой задается внешняя температура
/// Величина температуры влияет на вязкость и в итоге на гидравлические потери
/// Влияние температуры на плотность не учитывается
/// Использование: неизотермические расчеты
class PipeModelPGConstAreaNonIsothermal : public PipeModelPGConstArea
{
    using pde_t<2>::equation_coeffs_type;
    using pde_t<2>::right_party_type;
    using pde_t<2>::var_type;
private:
    /// @brief Температура посчитанная неким внешним алгоритимом. В данной модели предполагается заданной
    const std::vector<double>& temperature;
public:
    PipeModelPGConstAreaNonIsothermal(const pipe_properties_t& pipe, 
        const oil_parameters_t& oil, const std::vector<double>& temperature)
        : PipeModelPGConstArea(pipe, oil)
        , temperature(temperature)
    {

    }
    virtual var_type getSourceTerm(size_t grid_index, const var_type& point_vector) const override
    {
        double p = point_vector[0];
        double G = point_vector[1];
        double rho = oil.density();
        double S_0 = pipe.wall.getArea();
        double v = G / (rho * S_0);
        double Re = v * pipe.wall.diameter / oil.viscosity(temperature[grid_index]);

        double lambda = pipe.resistance_function(Re);
        //double lambda = hydraulic_resistance_shifrinson(Re, pipe.wall.relativeRoughness());
        double tau_w = lambda / 8 * rho * v * abs(v);
        double s1 = -M_PI * pipe.wall.diameter * tau_w;

        var_type s = { 0, s1 };
        return s;
    }
};


/// @brief Уравнение трубы для давления и объемного расхода
/// Не учитываются эффекты сжимаемости и объемного расширения
/// Учитывается партийность, неизотермичность
/// В источниковый член перенесена конвекция импульса, вызванная сменой плотности
/// Описание в документе "Уравнения для PQ"
class PipeModelPQConstAreaSortedNonisothermal : public pde_t<2>
{
    using pde_t<2>::equation_coeffs_type;
    using pde_t<2>::right_party_type;
    using pde_t<2>::var_type;
protected:
    /// @brief Параметры трубы
    const pipe_properties_t& pipe;
    /// @brief Профиль свойств жидкости
    const fluid_properties_profile_t& oil;
    /// @brief Профиль температуры
    const std::vector<double>& temperature;

public:
    PipeModelPQConstAreaSortedNonisothermal(
        const pipe_properties_t& pipe, const fluid_properties_profile_t& oil, 
        const std::vector<double>& temperature)
        : pipe(pipe)
        , oil(oil)
        , temperature(temperature)
    {
    }

    /// @brief Возвращает известную уравнению сетку
    virtual const std::vector<double>& get_grid() const override {
        return pipe.profile.coordinates;
    }

    /// @brief Получение матрицы коэффициентов системы уравнений (Row-major)
    /// \param curr
    /// \param index
    /// \return Row-major матрица, массив вектор-строк
    virtual equation_coeffs_type getEquationsCoeffs(
        size_t grid_index, const var_type& point_vector) const override
    {
        double S_0 = pipe.wall.getArea();
        double density = oil.nominal_density[grid_index];

        double beta_S = pipe.wall.getCompressionRatio();
        double beta_rho = oil.get_compression_ratio();

        equation_coeffs_type A; // Row-major матрица, массив вектор-строк
        A[0] = { 0, 1 / (S_0 * (beta_S + beta_rho)) };
        A[1] = { S_0 / density, 0 };
        return A;
    }
    /// @brief Обратная матрица коэффициентов системы уравнений
    /// @param curr 
    /// @param index 
    /// @return 
    virtual equation_coeffs_type getEquationsCoeffsInv(
        size_t grid_index, const var_type& point_vector) const override
    {
        double S_0 = pipe.wall.getArea();
        double density = oil.nominal_density[grid_index];

        double beta_S = pipe.wall.getCompressionRatio();
        double beta_rho = oil.get_compression_ratio();

        std::array<std::array<double, 2>, 2> Ainv;
        Ainv[0] = { 0, density / S_0 };
        Ainv[1] = { S_0 * (beta_S + beta_rho), 0 };

        return Ainv;
    }



    /// @brief Получение вектора правой части системы уравнений
    virtual var_type getSourceTerm(size_t grid_index, const var_type& point_vector) const override
    {
        double p = point_vector[0];
        double Q = point_vector[1];
        double rho = oil.nominal_density[grid_index];

        double d = pipe.wall.diameter;
        double S_0 = pipe.wall.getArea();
        double da = pipe.adaptation.diameter;
        d *= da;
        S_0 *= da * da; // квадрат

        double v = Q / S_0;

        double T = temperature[grid_index];
        double Re = v * d / oil.get_viscosity(grid_index, T);
        double lambda = pipe.resistance_function(Re);
        lambda *= pipe.adaptation.friction;
        double tau_w = lambda / 8 * rho * v * abs(v);

        double height_gradient; // dz/dx
        double density_gradient; // d(\rho)/dx
        const double* density_profile = &oil.nominal_density[grid_index];
        const double* height_profile = &pipe.profile.heights[grid_index];
        const double* grid = &(get_grid()[grid_index]);
        if (grid_index == 0) {
            density_gradient = (density_profile[+1] - density_profile[0]) / (grid[+1] - grid[0]);
            height_gradient = (height_profile[+1] - height_profile[0]) / (grid[+1] - grid[0]);
        }
        else if (grid_index == pipe.profile.get_point_count() - 1) {
            density_gradient = (density_profile[0] - density_profile[-1]) / (grid[0] - grid[-1]);
            height_gradient = (height_profile[0] - height_profile[-1]) / (grid[0] - grid[-1]);
        }
        else {
            density_gradient = (density_profile[+1] - density_profile[-1]) / (grid[+1] - grid[-1]);
            height_gradient = (height_profile[+1] - height_profile[-1]) / (grid[+1] - grid[-1]);
        }

        double s1 =
            2 * S_0 * v * abs(v) / rho * density_gradient
            - M_PI * d * tau_w / rho
            - M_G * S_0 * height_gradient;

        var_type s = { 0, s1 };
        return s;
    }
    virtual std::pair<var_type, equation_coeffs_type> GetLeftEigens(
        size_t index, const var_type& u) const {
        throw std::logic_error("not impl");
    }

    virtual std::pair<var_type, equation_coeffs_type> GetRightEigens(
        size_t index, const var_type& u) const {
        throw std::logic_error("not impl");
    }

    virtual var_type GetRightEigenVector(
        size_t profile_index, size_t eigen_index, const var_type& u) const {
        throw std::logic_error("not impl");
    }

    virtual double get_wave_strength(
        size_t profile_index, size_t eigen_index, const var_type& u) const {
        throw std::logic_error("not impl");
    }

};

/// @brief Уравнение трубы для задачи PQ с учетом движения партий
/// Учитывается, что параметры партий могут задавать в точках, и в ячейках
class isothermal_pipe_PQ_parties_t : public ode_t<1>
{
public:
    using ode_t<1>::equation_coeffs_type;
    using ode_t<1>::right_party_type;
    using ode_t<1>::var_type;
protected:
    const std::vector<double>& rho_profile;
    const std::vector<double>& nu_profile;
    const pipe_properties_t& pipe;
    const double flow;
    const int solver_direction;
    size_t flag_for_points;
public:
    /// @brief Констуктор уравнения трубы
    /// @param pipe Ссылка на сущность трубы
    /// @param oil Ссылка на сущность нефти
    /// @param flow Объемный расход
    /// @param solver_direction Направление расчета по Эйлеру, должно обязательно совпадать с параметром солвера Эйлера
    isothermal_pipe_PQ_parties_t(const pipe_properties_t& pipe, const std::vector<double>& rho_profile, const std::vector<double>& nu_profile, double flow,
        int solver_direction, size_t flag_for_points = 0)
        : pipe(pipe)
        , rho_profile(rho_profile)
        , nu_profile(nu_profile)
        , flow(flow)
        , solver_direction(solver_direction)
        , flag_for_points(flag_for_points)
    {}

    /// @brief Возвращает известную уравнению сетку
    virtual const std::vector<double>& get_grid() const override {
        return pipe.profile.coordinates;
    }

    /// @brief Возвращает значение правой части ДУ
    /// @param grid_index Обсчитываемый индекс расчетной сетки
    /// @param point_vector Начальные условия
    /// @return Значение правой части ДУ в точке point_vector
    virtual right_party_type ode_right_party(
        size_t grid_index, const var_type& point_vector) const override
    {

        /// Обработка индекса в случае расчетов на границах трубы
        /// Чтобы не выйти за массив высот, будем считать dz/dx в соседней точке
        size_t rheo_index = grid_index;

        if (pipe.profile.get_point_count() == rho_profile.size())
        {
            // Случай расчета партий в точках (например для метода характеристик)
            if (solver_direction == +1)
                rheo_index += 1;
            else
                rheo_index -= 1;
        }
        else
        {
            // Случай расчета партий в ячейках (например для quickest ultimate) 
            rheo_index = solver_direction == +1
                ? grid_index
                : grid_index - 1;
        }
        double rho = rho_profile[rheo_index];
        double S_0 = pipe.wall.getArea();
        double v = flow / (S_0);
        double Re = v * pipe.wall.diameter / nu_profile[rheo_index];
        double lambda = pipe.resistance_function(Re);
        double tau_w = lambda / 8 * rho * v * abs(v);
        double height_derivative = pipe.profile.get_height_derivative(grid_index, solver_direction);
        double result = -4 * tau_w / pipe.wall.diameter - rho * M_G * height_derivative;
        return result;
    }
};

/// @brief Уравнение сохранение импульса с учетом беспартийной перекачки, но с учетом неизотермичности
class nonisothermal_pipe_PQ_noparties_t : public ode_t<1>
{
public:
    using ode_t<1>::equation_coeffs_type;
    using ode_t<1>::right_party_type;
    using ode_t<1>::var_type;
protected:
    /// @brief Профиль температур
    const std::vector<double>& temperature_profile;
    /// @brief Параметра нефти
    const oil_parameters_t& oil;
    /// @brief Параметры трубы
    const pipe_properties_t& pipe;
    /// @brief Объемный расход
    const double flow;
    /// @brief Направление расчета
    const int solver_direction;
    size_t flag_for_points;
public:
    /// @brief Констуктор уравнения трубы
    /// @param pipe Ссылка на сущность трубы
    /// @param oil Ссылка на сущность нефти
    /// @param flow Объемный расход
    /// @param solver_direction Направление расчета по Эйлеру, должно обязательно совпадать с параметром солвера Эйлера
    nonisothermal_pipe_PQ_noparties_t(const pipe_properties_t& pipe, const oil_parameters_t& oil, const std::vector<double>& temperature_profile, double flow,
        int solver_direction, size_t flag_for_points = 0)
        : pipe(pipe)
        , oil(oil)
        , temperature_profile(temperature_profile)
        , flow(flow)
        , solver_direction(solver_direction)
        , flag_for_points(flag_for_points)
    {
    }

    /// @brief Возвращает известную уравнению сетку
    virtual const std::vector<double>& get_grid() const override {
        return pipe.profile.coordinates;
    }

    /// @brief Возвращает значение правой части ДУ
    /// @param grid_index Обсчитываемый индекс расчетной сетки
    /// @param point_vector Начальные условия
    /// @return Значение правой части ДУ в точке point_vector
    virtual right_party_type ode_right_party(
        size_t grid_index, const var_type& point_vector) const override
    {

        /// Обработка индекса в случае расчетов на границах трубы
        /// Чтобы не выйти за массив высот, будем считать dz/dx в соседней точке
        size_t temp_index = grid_index;

        if (pipe.profile.get_point_count() == temperature_profile.size())
        {
            // Случай расчета партий в точках (например для метода характеристик)
            if (solver_direction == +1)
                temp_index += 1;
            else
                temp_index -= 1;
        }
        else
        {
            // Случай расчета партий в ячейках (например для quickest ultimate) 
            temp_index = solver_direction == +1
                ? grid_index
                : grid_index - 1;
        }
        double rho = oil.density.nominal_density;
        double S_0 = pipe.wall.getArea();
        double v = flow / (S_0);
        double Re = v * pipe.wall.diameter / oil.viscosity(temperature_profile[temp_index]);
        //double check = oil.viscosity(temperature_profile[reo_index]);
        double lambda = pipe.resistance_function(Re);
        double tau_w = lambda / 8 * rho * v * abs(v);
        double height_derivative = pipe.profile.get_height_derivative(grid_index, solver_direction);
        double result = -4 * tau_w / pipe.wall.diameter - rho * M_G * height_derivative;
        return result;
    }
};

}