#pragma once

/// @brief Структура для хранения массы в нулевой момент времени, на входе и на выходе трубы
struct mass_balance_data_t {
    /// @brief Масса флюида в трубе на начальный момент времени
    double M0 = 0.0;
    /// @brief Масса флюида поступающего на вход трубы за время моделирования
    double M_in = 0.0;
    /// @brief Масса флюида выходящего из трубы за время моделирования
    double M_out = 0.0;
};

/// @brief класс для вычисления массы флюида в трубе
class mass_data_processor {
private:
    /// @brief Объем ячейки
    double V_cell;
    /// @brief Временной шаг
    double dt;
    /// @brief Для отслеживания массы в балансовом расчете
    mass_balance_data_t mass;
public:
    /// @brief Вычисляет объём одной расчётной ячейки на основе геометрических параметров трубы.
    /// @param pipe параметры трубы
    /// @param time_step шаг расчета
    mass_data_processor(const pde_solvers::simple_pipe_properties& pipe, 
        double time_step)
        : dt(time_step)
    {
        V_cell = (M_PI * pipe.diameter * pipe.diameter / 4.0) * pipe.dx;
    }
    /// @brief Вычисляет массу флюида в трубе на текущем слое
    /// @param layer параметры расчетного слоя
    /// @return
    double calc_pipe_mass(const pde_solvers::qsm_advection_layer& layer) const {
        double total_mass = 0;
        for (double ρ : layer.value)
            total_mass += ρ * V_cell;
        return total_mass;
    }
    /// @brief Текущий баланс
    /// @return
    const mass_balance_data_t& get_current_balance() const {
        return mass;
    }
    /// @brief Обрабатывает данные слоя на текущем временном шаге
    void process_data(size_t step_index,
        const pde_solvers::qsm_advection_layer& layer)
    {
        if (step_index == 0) {
            mass.M0 = calc_pipe_mass(layer);
        }
        else {
            double Q = layer.volumetric_flow;
            double ρ_in = layer.value_in;
            double ρ_out = layer.value_out;
            mass.M_in += ρ_in * Q * dt;
            mass.M_out += ρ_out * Q * dt;
        }
    }
};
/// @brief Исходные параметры трубы
/// @return
simple_pipe_properties prepare_conservation_pipe() {
    simple_pipe_properties pipe;
    pipe.length = 50e3;
    pipe.dx = 200.0;
    pipe.diameter = 0.7;
    return pipe;
}

///@brief Тест на проверку консервативности схемы QUICKEST Ultimate.
TEST(QuickestUltimate, HasConservativity) {
    // Физические параметры задачи
    const double ρ1 = 600.0;     ///< Плотность до переключения [кг/м³]
    const double ρ2 = 650.0;     ///< Плотность после переключения [кг/м³]
    const double Q = 0.5;        ///< Постоянный объёмный расход [м³/с]
    const double dt = 60;        ///< Временной шаг расчёта [с]
    //Временные параметры моделирования
    double T = 120 * 3600;       ///< Общее время моделирования [с]
    double t_sw = 6 * 3600;      ///< Время переключения плотности на входе [с]
    //Дискретизация по времени
    const int steps = static_cast<int>(T / dt);
   /* const int sw_step = static_cast<int>(t_sw / dt);   */
    // Подготовка временной сетки
    std::vector<double> times(steps + 1);
    for (int i = 0; i <= steps; ++i) {
        times[i] = i * dt;
    }
    // Подготовка граничных условий
    std::vector<pde_solvers::qsm_advection_task_boundaries_t> boundaries(steps + 1);
    for (int i = 0; i <= steps; ++i) {
        boundaries[i].volumetric_flow = Q;
        boundaries[i].value = (i * dt <= t_sw) ? ρ1 : ρ2;
    }

    //Инициализация геометрии трубы
    simple_pipe_properties sp = prepare_conservation_pipe();
    //Создание объектов для расчёта
    pde_solvers::qsm_advection_task_t task(sp);
    mass_data_processor processor(sp, dt);

    //Получение результатов расчёта
    auto& final_pipe_state = task.get_buffer().current();
    double M_final = processor.calc_pipe_mass(final_pipe_state);
    mass_balance_data_t M = processor.get_current_balance();
    double M_ref = std::max({
        std::abs(M.M0),
        std::abs(M.M_in),
        std::abs(M.M_out),
        1.0
        });
    double tol = std::max(1e-6, 1e-10 * M_ref);
    EXPECT_NEAR(M_final, M.M0 + M.M_in - M.M_out, tol);
}
