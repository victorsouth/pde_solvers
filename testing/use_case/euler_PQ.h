#pragma once

/// @brief Уравнение трубы для задачи PQ
class Pipe_model_for_PQ_t : public ode_t<1>
{
public:
    using ode_t<1>::equation_coeffs_type;
    using ode_t<1>::right_party_type;
    using ode_t<1>::var_type;
protected:
    pipe_properties_t& pipe;
    oil_parameters_t& oil;
    double flow;

public:
    /// @brief Констуктор уравнения трубы
    /// @param pipe Ссылка на сущность трубы
    /// @param oil Ссылка на сущность нефти
    /// @param flow Объемный расход
    Pipe_model_for_PQ_t(pipe_properties_t& pipe, oil_parameters_t& oil, double flow)
        : pipe(pipe)
        , oil(oil)
        , flow(flow)
    {

    }

    /// @brief Возвращает известную уравнению сетку
    virtual const vector<double>& get_grid() const override {
        return pipe.profile.coordinates;
    }

    /// @brief Возвращает значение правой части ДУ
    /// см. файл 2023-11-09 Реализация стационарных моделей с прицелом на квазистационар.docx
    /// @param grid_index Обсчитываемый индекс расчетной сетки
    /// @param point_vector Начальные условия
    /// @return Значение правой части ДУ в точке point_vector
    virtual right_party_type ode_right_party(
        size_t grid_index, const var_type& point_vector) const override
    {
        double rho = oil.density();
        double S_0 = pipe.wall.getArea();
        double v = flow / (S_0);
        double Re = v * pipe.wall.diameter / oil.viscosity();
        double lambda = pipe.resistance_function(Re, pipe.wall.relativeRoughness());
        double tau_w = lambda / 8 * rho * v * abs(v);

        /// Обработка индекса в случае расчетов на границах трубы
        /// Чтобы не выйти за массив высот, будем считать dz/dx в соседней точке
        grid_index = grid_index == 0 ? grid_index + 1 : grid_index;
        grid_index = grid_index == pipe.profile.heights.size() - 1 ? grid_index - 1 : grid_index;

        double height_derivative = (pipe.profile.heights[grid_index] - pipe.profile.heights[grid_index - 1]) /
            (pipe.profile.coordinates[grid_index] - pipe.profile.coordinates[grid_index - 1]);

        return { ((-4) / pipe.wall.diameter) * tau_w - rho * M_G * height_derivative };
    }

};

TEST(EulerPQ, UseCase)
{
    /// Создаем сущность трубы
    pipe_properties_t pipe;

    /// Задаем сетку трубы
    pipe.profile.coordinates = { 0, 1000, 2000 };

    /// Задаем высотные отметки и 
    pipe.profile.heights = vector<double>(pipe.profile.coordinates.size(), 0);

    /// Создаем буфер из двух слоев, каждый совпадает по размеру с pipe.profile.heights
    ring_buffer_t<vector<double>> buffer(2, pipe.profile.heights);

    /// Создаем сущность нефти
    oil_parameters_t oil;

    /// Задаем объемнй расход нефти, [м3/с]
    double Q = 0.8;

    /// Создаем расчетную модель трубы
    Pipe_model_for_PQ_t pipeModel(pipe, oil, Q);

    /// Получаем указатель на начало слоя в буфере
    profile_wrapper<double, 1> start_layer(buffer.current());

    ///Задаем начальное давление
    double Pout = 5e5;

    /// Модифицированный метод Эйлера для модели pipeModel,
    /// расчет ведется справа-налево относительно сетки,
    /// начальное условие Pout, 
    /// результаты расчета запишутся в слой, на который указывает start_layer
    solve_euler_corrector<1>(pipeModel, -1, Pout, &start_layer);
}