#pragma once


namespace pde_solvers {

/// @brief Свойства конденсатопровода
/// Расширяет базовые свойства трубы добавлением кинематической вязкости
struct condensate_pipe_properties_t : public pipe_properties_t {
    /// @brief Кинематическая вязкость, м²/с
    double kinematic_viscosity{1e-7};
};

/// @brief Уравнение сохранения импульса для конденсатопровода с учетом движения партий
/// Используется для задач PQ (задан расход, найти профиль давления) и PP (заданы давления, найти расход методом Ньютона)
/// Учитывает гидравлические потери на трение и влияние гравитации
class condensate_pipe_PQ_parties_t : public ode_t<1>
{
public:
    using ode_t<1>::equation_coeffs_type;
    using ode_t<1>::right_party_type;
    using ode_t<1>::var_type;
protected:
    /// @brief Профиль плотности вдоль трубы, кг/м³
    const std::vector<double>& rho_profile;
    /// @brief Свойства конденсатопровода
    const condensate_pipe_properties_t& pipe;
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
    condensate_pipe_PQ_parties_t(const condensate_pipe_properties_t& pipe,
        const std::vector<double>& rho_profile, 
        double flow,
        int solver_direction)
        : pipe(pipe)
        , rho_profile(rho_profile)
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
        double Re = v * pipe.wall.diameter / pipe.kinematic_viscosity;
        double lambda = pipe.resistance_function(Re);
        double tau_w = lambda / 8 * rho * v * std::abs(v);
        double height_derivative = pipe.profile.get_height_derivative(grid_index, solver_direction);
        double result = -4 * tau_w / pipe.wall.diameter - rho * M_G * height_derivative;
        return result;
    }
};

}

