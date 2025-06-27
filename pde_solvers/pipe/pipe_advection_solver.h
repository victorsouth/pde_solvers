#pragma once

namespace pde_solvers {

/// @brief Решатель транспортного уравнения методом характеристик, 
/// при этом считается, что скорость по длине трубопровода постоянна,
/// а число Куранта всегда равно единице
class advection_moc_solver
{
public:
    /// @brief Конструктор транспортного солвера
    /// @param pipe Параметры трубопровода 
    /// @param vol_flow Объёмный расход
    /// @param prev Предыдущий слой
    /// @param next Новый слой
    advection_moc_solver(const pipe_properties_t& pipe, double vol_flow,
        std::vector<double>& prev, std::vector<double>& next)
        : pipe{ pipe }
        , volumetric_flow{ vol_flow }
        , prev{ prev }
        , next{ next }
    {}

    /// @brief Расчёт нового слоя
    /// @param dt Шаг моделирования
    /// @param par_in Значение параметра среды, втекающей в начало трубопровода
    /// @param par_out Значение параметра среды, втекающей в конец трубопровода при обратном течении
    void step(const double dt,const double par_in, const double par_out)
    {
        double p = interpolation_offset(dt);

        int direction = get_eigen_value() > 0 ? 1 : -1;
        size_t start_index = direction > 0 ? 1 : (next.size()) - 2;
        size_t end_index = direction < 0 ? -1 : next.size();
        next[start_index - direction] = direction > 0 ? par_in : par_out;
        for (size_t index = start_index; index != end_index; index += direction)
        {
            next[index] = prev[index - direction] * p + prev[index] * (1 - p);
        }
    }

    /// @brief Расчёт шага по времени, при котором Курант равен единице (Cr = 1)
    double prepare_step() const
    {
        const std::vector<double>& grid = pipe.profile.coordinates;

        double max_eigen = get_eigen_value();
        double dx = grid[1] - grid[0];
        double courant_step = dx / max_eigen;

        return abs(courant_step);
    }

protected:

    /// @brief модель трубопровода
    const pipe_properties_t& pipe;
    /// @brief Объемный расход
    const double volumetric_flow;
    /// @brief Предыдущий слой
    std::vector<double>& prev;
    /// @brief Новый слой
    std::vector<double>& next;


    /// @brief Расчёт собственного значения
    double get_eigen_value() const
    {
        double S = pipe.wall.getArea();
        return volumetric_flow / S;
    }

    /// @brief Расчёт смещения в случае Cr < 1
    /// @param dt Шаг моделирования
    double interpolation_offset(const double dt) const
    {
        return dt / prepare_step();
    }
};


}