#pragma once


namespace pde_solvers {

/// @brief Свойства конденсатопровода
/// Расширяет базовые свойства трубы добавлением кинематической вязкости
struct iso_nonbarotropic_pipe_properties_t : public pipe_properties_t {
    /// @brief Кинематическая вязкость, м²/с
    double kinematic_viscosity{1e-7};

    /// @brief Конструктор по умолчанию
    iso_nonbarotropic_pipe_properties_t() = default;

    /// @brief Объект со значениями по умолчанию
    static iso_nonbarotropic_pipe_properties_t default_values() {
        iso_nonbarotropic_pipe_properties_t result;
        double length = 5000;
        double dx = 200;
        result.profile.coordinates =
            pde_solvers::pipe_profile_uniform::generate_uniform_grid(0.0, length, dx);
        result.wall.diameter = 1;       
        return result;
    }

    /// @brief Конструктор из JSON данных
    /// @param json_data Данные о трубе в формате JSON
    iso_nonbarotropic_pipe_properties_t(const pde_solvers::pipe_json_data& json_data)
    {
        *this = iso_nonbarotropic_pipe_properties_t::default_values();
        profile.coordinates = { json_data.x_start, json_data.x_end };
        wall.diameter = json_data.diameter;
    }

    /// @brief Создает равномерный профиль трубы с заданным шагом
    /// TODO: Метод одинаков для всех профилей. Убрать дублирование
    /// @param desired_dx Желаемый шаг сетки, м
    void make_uniform_profile(double desired_dx) {
        profile.coordinates = pde_solvers::pipe_profile_uniform::generate_uniform_grid(
            profile.coordinates.front(), profile.get_length(), desired_dx
        );

        // TODO: заменить на фактический профиль
        profile.heights = std::vector<double>(profile.get_point_count());
    }
};

/// @brief Уравнение сохранения импульса для конденсатопровода с учетом движения партий
/// Используется для задач PQ (задан расход, найти профиль давления) и PP (заданы давления, найти расход методом Ньютона)
/// Учитывает гидравлические потери на трение и влияние гравитации
class iso_nonbaro_impulse_equation_t : public ode_t<1>
{
public:
    /// @brief Тип коэффициентов уравнения
    using ode_t<1>::equation_coeffs_type;
    /// @brief Тип правой части дифференциального уравнения
    using ode_t<1>::right_party_type;
    /// @brief Тип переменной дифференциального уравнения
    using ode_t<1>::var_type;
protected:
    /// @brief Профиль плотности вдоль трубы, кг/м³
    const std::vector<double>& density_profile;
    /// @brief Свойства конденсатопровода
    const iso_nonbarotropic_pipe_properties_t& pipe;
    /// @brief Объемный расход, м³/с
    const double flow;
    /// @brief Направление расчета по Эйлеру (+1 - от входа к выходу, -1 - от выхода ко входу)
    const int solver_direction;
public:
    /// @brief Констуктор уравнения трубы
    /// @param pipe Ссылка на сущность трубы
    /// @param rho_profile Профиль плотности
    /// @param flow Объемный расход
    /// @param solver_direction Направление расчета по Эйлеру, должно обязательно совпадать с параметром солвера Эйлера
    iso_nonbaro_impulse_equation_t(const iso_nonbarotropic_pipe_properties_t& pipe,
        const std::vector<double>& density_profile, 
        double flow,
        int solver_direction)
        : pipe(pipe)
        , density_profile(density_profile)
        , flow(flow)
        , solver_direction(solver_direction)
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
        size_t rheo_index = grid_index;

        if (pipe.profile.get_point_count() == density_profile.size())
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
        double rho = density_profile[rheo_index];
        double S_0 = pipe.wall.getArea();
        double v = flow / (S_0);
        double Re = v * pipe.wall.diameter / pipe.kinematic_viscosity;
        double lambda = pipe.resistance_function(Re);
        double tau_w = lambda / 8 * rho * v * std::abs(v);
        double height_derivative = pipe.profile.get_height_derivative(grid_index, solver_direction);
        double result = -4 * tau_w / pipe.wall.diameter - rho * M_G * height_derivative;
        return result;
    }
};

enum class improver_efficiency_formula_t
{
    polynomial = 1,
    /// @brief Белоусов 1986 формула (7.9) 
    virk = 2,
    /// @brief Егоров 2013 формула (2) Белоусов 1986 формула (7.13) 
    virk_modified = 3,
    /// @brief Кулик 2025 опытная формула
    potential = 4
};

/// @brief Параметры противотурбулентной присадки
struct improver_properties_t
{
    /// @brief коэффициенты модели присадки, меняют смысл в зависимости от efficiency_formula
    std::vector<double> coeffs;
    /// @brief формула расчета эффективности присадки
    improver_efficiency_formula_t efficiency_formula;
    
    /// @brief расчет эффективности присадки с учетом заданного типа efficiency_formula
    double drag_reduction_effectiveness(double concentration) const
    {
        switch (efficiency_formula)
        {
        case improver_efficiency_formula_t::polynomial:
        {
            if (coeffs.empty())
            {
                throw std::runtime_error("dimension of coeffs[] must be not empty");
            }
            double dr = 0;
            for (std::size_t i = 0; i < coeffs.size(); i++)
            {
                dr += coeffs[i] * std::pow(concentration, i);
            }

            return dr;
        }
        // Белоусов 1986 формула (7.9) 
        case improver_efficiency_formula_t::virk:
        {
            if (coeffs.size() != 2)
            {
                throw std::runtime_error("dimension of coeffs[] must be 2");
            }
            double dr_max = coeffs[0];
            double alpha = coeffs[1];

            double dr = dr_max * (alpha * concentration) / (1 + alpha * concentration);

            return dr;
        }
        //Егоров 2013 формула (2) Белоусов 1986 формула (7.13)     
        case improver_efficiency_formula_t::virk_modified:
        {
            if (coeffs.size() != 2)
            {
                throw std::runtime_error("dimension of coeffs[] must be 2");
            }
            double A = coeffs[0];
            double B = coeffs[1];

            double dr = concentration / (A + B * concentration);

            return dr;
        }
        //Кулик 2025 опытная формула
        case improver_efficiency_formula_t::potential:
        {
            if (coeffs.size() != 2)
            {
                throw std::runtime_error("dimension of coeffs[] must be 2");
            }
            double a = coeffs[0];
            double b = coeffs[1];

            double dr = a * std::pow(concentration, b);

            return dr;
        }
        default:
            throw std::runtime_error("Not implemented");
        }
    }
};

/// @brief Свойства конденсатопровода
/// Расширяет базовые свойства трубы добавлением кинематической вязкости и свойствами присадки
struct iso_nonbaro_improver_pipe_properties_t : public pipe_properties_t {
    /// @brief Кинематическая вязкость, м²/с
    double kinematic_viscosity{std::numeric_limits<double>::quiet_NaN()};
    /// @brief Параметры присадки
    improver_properties_t improver;

    /// @brief Конструктор по умолчанию
    iso_nonbaro_improver_pipe_properties_t() = default;

    /// @brief Объект со значениями по умолчанию
    static iso_nonbaro_improver_pipe_properties_t default_values() {
        throw std::runtime_error("not implemented");
    }

    /// @brief Конструктор из JSON данных
    /// @param json_data Данные о трубе в формате JSON
    iso_nonbaro_improver_pipe_properties_t(const pde_solvers::pipe_json_data& json_data)
    {
        throw std::runtime_error("not implemented");
    }

    /// @brief Создает равномерный профиль трубы с заданным шагом
    /// TODO: Метод одинаков для всех профилей. Убрать дублирование
    /// @param desired_dx Желаемый шаг сетки, м
    void make_uniform_profile(double desired_dx) {
        profile.coordinates = pde_solvers::pipe_profile_uniform::generate_uniform_grid(
            profile.coordinates.front(), profile.get_length(), desired_dx
        );

        // TODO: заменить на фактический профиль
        profile.heights = std::vector<double>(profile.get_point_count());
    }
};
/// @brief Уравнение сохранения импульса для конденсатопровода с учетом движения партий
/// Используется для задач PQ (задан расход, найти профиль давления) и PP (заданы давления, найти расход методом Ньютона)
/// Учитывает гидравлические потери на трение и влияние гравитации
class iso_nonbaro_improver_impulse_equation_t : public ode_t<1>
{
public:
    /// @brief Тип коэффициентов уравнения
    using ode_t<1>::equation_coeffs_type;
    /// @brief Тип правой части дифференциального уравнения
    using ode_t<1>::right_party_type;
    /// @brief Тип переменной дифференциального уравнения
    using ode_t<1>::var_type;
protected:
    /// @brief Профиль плотности вдоль трубы, кг/м³
    const std::vector<double>& density_profile;
    /// @brief Профиль концентрации присадки вдоль трубы
    const std::vector<double>& improver_concentration_profile;
    /// @brief Свойства конденсатопровода
    const iso_nonbaro_improver_pipe_properties_t& pipe;
    /// @brief Объемный расход, м³/с
    const double flow;
    /// @brief Направление расчета по Эйлеру (+1 - от входа к выходу, -1 - от выхода ко входу)
    const int solver_direction;
public:
    /// @brief Констуктор уравнения трубы
    /// @param pipe Ссылка на сущность трубы
    /// @param rho_profile Профиль плотности
    /// @param flow Объемный расход
    /// @param solver_direction Направление расчета по Эйлеру, должно обязательно совпадать с параметром солвера Эйлера
    iso_nonbaro_improver_impulse_equation_t(const iso_nonbaro_improver_pipe_properties_t & pipe,
        const std::vector<double>& density_profile,
        const std::vector<double>& improver_concentration_profile,
        double flow,
        int solver_direction)
        : pipe(pipe)
        , density_profile(density_profile)
        , improver_concentration_profile(improver_concentration_profile)
        , flow(flow)
        , solver_direction(solver_direction)
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
        size_t rheo_index = grid_index;

        if (pipe.profile.get_point_count() == density_profile.size())
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
        double rho = density_profile[rheo_index];
        double S_0 = pipe.wall.getArea();
        double v = flow / (S_0);
        double Re = v * pipe.wall.diameter / pipe.kinematic_viscosity;
        double improver_concentration = improver_concentration_profile[rheo_index];
        double DR = pipe.improver.drag_reduction_effectiveness(improver_concentration);
        double lambda = pipe.resistance_function(Re) * (1 - DR);
        double tau_w = lambda / 8 * rho * v * std::abs(v);
        double height_derivative = pipe.profile.get_height_derivative(grid_index, solver_direction);
        double result = -4 * tau_w / pipe.wall.diameter - rho * M_G * height_derivative;
        return result;
    }
};

}

