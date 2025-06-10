#pragma once
#include <string>
#include <vector>
#include <Eigen/Sparse>
#include "pipe/pipe_hydraulic_struct.h"
#include <fixed/fixed_optimizer.h>
//#include "../timeseries/timeseries_helpers.h"

namespace pde_solvers {
;

/// @brief Структура для хранения текущих значений параметров идентификации
struct ident_parameters_isothermal_qsm {
    /// @brief Параметр для идентификации по диаметру
    double diameter_adaptation{ 1 };
    /// @brief Параметр для идентификации по коэффициенту гидравлического сопротивления
    double friction_adaptation{ 1 };

    /// @brief Применение параметров идентификации в модели трубопровода
    /// @param nominal_pipe Номинальная (исходная) модель трубопровода
    /// @param pipe_to_ident Модель трубопровода, изменённая с учётом параметров идентификации
    void set_adaptation(const pipe_properties_t& nominal_pipe, pipe_properties_t* pipe_to_ident) const
    {
        pipe_to_ident->wall.diameter = diameter_adaptation * nominal_pipe.wall.diameter;
        pipe_to_ident->wall.resistance_function_adaptation = friction_adaptation;
    }

};
/// @brief Структура для выбора параметра идентификации
struct ident_isothermal_qsm_pipe_settings {
    /// @brief Флаг идентификации по диаметру
    bool ident_diameter{ false };
    /// @brief Флаг идентификации по коэффициенту гидравлического сопротивления
    bool ident_friction{ false };

    /// @brief Проверка на то, чтобы был выбран только один параметр, 
    /// так как проводить идентификацию сразу по двум этим параметрам не имеет смысла - они связаны между собой
    void check_parameters() const {
        // если ничего не задано или наоборот, все задано, кидаем исключение
        if (!(ident_diameter ^ ident_friction))
            throw std::runtime_error("Identification parameters are incorrectly selected");
    }

    /// @brief Распаковываем переданный в residuals VectorXd в зависимости от того, какой параметр был выбран
    /// @param packed_ident_parameters Переданный в функцию расчёта невязок набор аргументов
    /// @return Распакованный набор параметров идентификации
    ident_parameters_isothermal_qsm unpack_ident_parameters(const Eigen::VectorXd& packed_ident_parameters) const {
        
        ident_parameters_isothermal_qsm result;
        // Распаковываем набор аргументов
        size_t index = 0;
        if (ident_diameter) {
            result.diameter_adaptation = packed_ident_parameters(index++);
        }
        if (ident_friction) {
            result.friction_adaptation = packed_ident_parameters(index++);
        }
        return result;
    }
};

/// @brief Класс для идентификации изотермической квазистационарной модели трубопровода
/// на данных с реального ЛУ по одному из двух параметров:
/// 1. диаметру
/// 2. коэффициенту гидравлического сопротивления
class ident_isothermal_qsm_pipe_parameters_t : public fixed_least_squares_function_t
{
    /// @brief Настройка, определяющая по какому из параметров будет проводиться идентификация
    const ident_isothermal_qsm_pipe_settings settings;
    /// @brief Номинальная (исходная) модель трубопровода
    const pipe_properties_t pipe_nominal;
    /// @brief Модель трубопровода с учётом текущего значения параметра идентификации
    pipe_properties_t pipe_to_ident;

    /// @brief Предпосчитанная временная сетка для моделирования работы ЛУ
    const std::vector<double>& times;
    /// @brief Интерполированные значения временных рядов СДКУ 
    /// (давление, расход, плотность и вязкость в начале ЛУ) - краевые условия
    const std::vector<std::vector<double>>& control_data;
    /// @brief Интерполированные значения временного рядка СДКУ (давление в конце ЛУ) - эталонные значения
    const std::vector<double>& etalon_pressure;
public:
    /// @brief Конструктор класса для идентификации
    /// @param settings Настройка, определяющая выбор параметра идентификации
    /// @param pipe Исходная модель трубопровода
    /// @param times Предпосчитанная временная сетка
    /// @param control_data Предпосчитанные краевые условия
    /// @param etalon_pressure Предпосчитанные эталонные значения
    ident_isothermal_qsm_pipe_parameters_t(
        const ident_isothermal_qsm_pipe_settings& settings,
        const pipe_properties_t& pipe, const std::vector<double>& times, const std::vector<std::vector<double>>& control_data, const std::vector<double>& etalon_pressure)
        : settings(settings)
        , pipe_nominal{ pipe }
        , pipe_to_ident{ pipe }
        , times{ times }
        , control_data{ control_data }
        , etalon_pressure{ etalon_pressure }
    {
        // Проверяем коректность выбора параметра идентификации
        settings.check_parameters();
    }
protected:
    /// @brief Расчёт невязок - проводится моделирование работы ЛУ и считается отклонение расчётного и фактического давления в конце ЛУ 
    /// @param ident_parameters Текущие значения параметров адаптации
    /// @return Вектор невязок - расхождения расчётного и фактического давления на выхое ЛУ
    std::vector<double> calc_vector_residuals(const ident_parameters_isothermal_qsm& ident_parameters)
    {
        // Сущность для сбора расчётных данных давления в конце ЛУ
        isothermal_qsm_batch_Pout_collector_t collector(times);
        // Применяем текущие параметры идентификации на модели трубопровода
        ident_parameters.set_adaptation(pipe_nominal, &pipe_to_ident);

        // Проводим гидравлический изотермический квазистационарный расчёт
        isothermal_quasistatic_PQ_task_t<quickest_ultimate_fv_solver> task(pipe_to_ident);
        isothermal_quasistatic_batch<quickest_ultimate_fv_solver, isothermal_qsm_batch_Pout_collector_t::layer_type>(
            task,
            times,
            control_data,
            &collector
        );

        const vector<double>& calc_pressure = collector.get_pressure_out_calculated();
        vector<double> simulation_result(times.size());
        // Считаем расхождение расчёта и факта
        std::transform(calc_pressure.begin(), calc_pressure.end(), etalon_pressure.begin(), simulation_result.begin(),
            [](double etalon, double calc) { return etalon - calc; });


        return simulation_result;
    }
public:
    /// @brief Инициализация расчёта вектора невязок
    /// @param parameters Текущие значения параметров идентификации
    /// @return Вектор невязок - расхождения расчётного и фактического давления на выхое ЛУ
    virtual Eigen::VectorXd residuals(const Eigen::VectorXd& parameters) override {

        ident_parameters_isothermal_qsm ident_parameters = settings.unpack_ident_parameters(parameters);
        //ident_parameters.diameter_adaptation = parameters(0);


        std::vector<double> simulation_result = calc_vector_residuals(ident_parameters);

        Eigen::Map<Eigen::VectorXd> result(simulation_result.data(), simulation_result.size());

        return result;
    }

    /// @brief Запуск алгоритма идентификации
    /// @param result Сущность, хранящая в себе результаты выполнения алгоритма идентификации
    /// @param analysis Сущность, хранящая в себе данные для аналитики хода идентификации
    /// @return Оптимальное значение параметра идентификации
    double ident(fixed_optimizer_result_t* result = nullptr, fixed_optimizer_result_analysis_t* analysis = nullptr) {
        // Настроечные параметры процесса идентификации
        fixed_optimizer_parameters_t parameters;
        if (analysis != nullptr) {
            parameters.analysis.objective_function_history = true;
            parameters.analysis.steps = true;
            parameters.analysis.argument_history = true;
        }

        fixed_optimizer_result_t result_local;
        if (result == nullptr) {
            result = &result_local;
        }

        // Задание начального значения параметра идентификации 
        Eigen::VectorXd initial_ident_parameter(1);
        initial_ident_parameter(0) = 1;
        // Запуск алгоритма идентификации 
        fixed_optimize_gauss_newton::optimize(*this, initial_ident_parameter, parameters, result, analysis);

        if (result->result_code == numerical_result_code_t::Converged) {
            return result->argument(0);
        }
        else {
            return std::numeric_limits<double>::quiet_NaN();
        }
    }
};

}
