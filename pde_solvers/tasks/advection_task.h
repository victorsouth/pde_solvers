#pragma once

namespace pde_solvers {
;

/// @brief Проблемно-ориентированный слой для расчета адвекции одного параметра
struct qsm_advection_layer {
    /// @brief Объемный расход
    double volumetric_flow;
    /// @brief Значение параметра на входе (плотность)
    double value_in;
    /// @brief Значение параметра на выходе (плотность)
    double value_out;
    /// @brief Профиль параметра (ячейки)
    std::vector<double> value;
    /// @brief Профиль вспомогательных расчетов для метода конечных объемов
    quickest_ultimate_fv_solver_traits<1>::specific_layer specific;
    /// @brief Инициализация профилей
    /// @param point_count Количество точек
    qsm_advection_layer(size_t point_count)
        : value(point_count - 1)
        , specific(point_count)
    {
    }

    /// @brief Подготовка обертки слоя для расчета методом конечных объемов
    static quickest_ultimate_fv_wrapper<1> get_value_wrapper(qsm_advection_layer& layer)
    {
        return quickest_ultimate_fv_wrapper<1>(layer.value, layer.specific);
    }
};

/// @brief Структура, содержащая в себе краевые условия задачи адвекиции
struct qsm_advection_task_boundaries_t {
    /// @brief Объемный расход
    double volumetric_flow;
    /// @brief Значение параметра на входе/выходе
    double value;
    /// @brief Конструктор по умолчанию
    qsm_advection_task_boundaries_t() = default;


    /// @brief Создание структуры со значениями по умолчанию
    static qsm_advection_task_boundaries_t default_values() {
        qsm_advection_task_boundaries_t result;
        result.volumetric_flow = 0.2;
        result.value = 600;// значение плотности
        return result;
    }
};


/// @brief Расчетная задача (task) для расчета адвекции Quickest-Ultimate
class qsm_advection_task_t {
public:
    /// @brief Тип слоя
    using layer_type = qsm_advection_layer;
    /// @brief Тип буфера
    using buffer_type = ring_buffer_t<layer_type>;
    /// @brief Тип граничных условий
    using boundaries_type = qsm_advection_task_boundaries_t;
private:
    /// @brief Модель трубы
    pipe_properties_t pipe;
    /// @brief Буфер для двух слоев адвекции (текущий и предыдущий)
    buffer_type buffer;

public:
    /// @brief Конструктор
    /// @param pipe Модель трубопровода
    qsm_advection_task_t(const simple_pipe_properties& pipe)
        : pipe(pipe_properties_t::build_simple_pipe(pipe))
        , buffer(2 /*количество слоев*/, pipe.get_point_count())
    {

    }

    /// @brief Начальный стационарный расчёт. 
    /// Ставим по всей трубе параметр из initial_conditions
    void solve(const qsm_advection_task_boundaries_t& initial_conditions)
    {
        auto& current = buffer.current();
        for (double& value : current.value) {
            value = initial_conditions.value;
        }
        current.volumetric_flow = initial_conditions.volumetric_flow;
        current.value_in = current.value_out = initial_conditions.value;
    }
    

public:
    /// @brief Рассчёт шага моделирования, включающий в себя расчёт шага движения партии
    /// Функция делает сдвиг буфера (advance) так, что buffer.current после вызова содержит свежерассчитанный слой
    /// @param dt временной шаг моделирования
    /// @param boundaries Краевые условия
    void step(double dt, const qsm_advection_task_boundaries_t& boundaries) {
        // считаем партии с помощью QUICKEST-ULTIMATE
        size_t n = pipe.profile.get_point_count();

        // Используем move-семантику для избежания лишних копирований
        std::vector<double> Q_profile;
        Q_profile.reserve(n);
        Q_profile.assign(n, boundaries.volumetric_flow); /// задаем по трубе новый расход из временного ряда

        // Создаем модель переноса
        PipeQAdvection advection_model(pipe, std::move(Q_profile));

        // Сдвигаем буфер
        buffer.advance(+1);
        auto value_wrapper = buffer.get_buffer_wrapper(&qsm_advection_layer::get_value_wrapper);

        // Создаем и запускаем решатель
        quickest_ultimate_fv_solver solver_rho(advection_model, value_wrapper);
        solver_rho.step(dt, boundaries.value, boundaries.value);

        // Обновляем текущий слой
        auto& current = buffer.current();
        current.volumetric_flow = boundaries.volumetric_flow;

        // Устанавливаем граничные значения в зависимости от направления потока
        if (current.volumetric_flow >= 0) {
            current.value_in = boundaries.value;
            current.value_out = current.value.back();
        }
        else {
            current.value_out = boundaries.value;
            current.value_in = current.value.front();
        }
    }

    /// @brief Возвращает ссылку на буфер
    auto& get_buffer()
    {
        return buffer;
    }
};

}
