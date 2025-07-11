﻿#pragma once

namespace pde_solvers {
;

/// @brief Структура для хранения текущих значений параметров идентификации
struct ident_parameters_nonisothermal_qsm {
    /// @brief Параметр для идентификации по коэффициенту теплообмена
    double htc{ 1 };

    /// @brief Применение параметров идентификации в модели трубопровода
    /// @param nominal_pipe Номинальная (исходная) модель трубопровода
    /// @param pipe_to_ident Модель трубопровода, изменённая с учётом параметров идентификации
    void set_parameter(const pipe_noniso_properties_t& nominal_pipe, pipe_noniso_properties_t* pipe_to_ident) const
    {
        pipe_to_ident->heat.ambient_heat_transfer = htc;
    }

};

/// @brief Структура для выбора параметра идентификации
struct ident_nonisothermal_qsm_pipe_settings {
    /// @brief Флаг идентификации по коэффициенту теплообмена
    bool ident_htc{ false };

    /// @brief Проверка на то, чтобы параметр идентификации был задан
    void check_parameters() const {
        // если не задан, кидаем исключение
        if (!(ident_htc))
            throw std::runtime_error("Identification parameters are incorrectly selected");
    }

    /// @brief Распаковываем переданный в residuals VectorXd в зависимости от того, какой параметр был выбран
    /// @param packed_ident_parameters Переданный в функцию расчёта невязок набор аргументов
    /// @return Распакованный набор параметров идентификации
    ident_parameters_nonisothermal_qsm unpack_ident_parameters(const Eigen::VectorXd& packed_ident_parameters) const {

        ident_parameters_nonisothermal_qsm result;
        // Распаковываем набор аргументов
        size_t index = 0;
        if (ident_htc) {
            result.htc = packed_ident_parameters(index++);
        }
        return result;
    }
};

/// @brief Класс для идентификации неизотермической квазистационарной модели трубопровода
/// на данных с реального ЛУ по коэффициенту теплообмена
class ident_nonisothermal_qsm_pipe_parameters_t : public fixed_least_squares_function_t
{
    /// @brief Настройка, определяющая по какому из параметров будет проводиться идентификация
    const ident_nonisothermal_qsm_pipe_settings settings;
    /// @brief Номинальная (исходная) модель трубопровода
    const pipe_noniso_properties_t pipe_nominal;
    /// @brief Модель трубопровода с учётом текущего значения параметра идентификации
    pipe_noniso_properties_t pipe_to_ident;
    
    /// @brief Параметры нефти
    const oil_parameters_t oil;
    /// @brief Предпосчитанная временная сетка для моделирования работы ЛУ
    const std::vector<double>& times;
    /// @brief Интерполированные значения временных рядов СДКУ 
    /// (темп, расход, плотность и вязкость в начале ЛУ) - краевые условия
    const std::vector<std::vector<double>>& control_data;
    /// @brief Интерполированные значения временного рядка СДКУ (темп в конце ЛУ) - эталонные значения
    const std::vector<double>& etalon_temperature;
    /// @brief Режим разделения задачи на Адвекци, Шухова и Шухова с Адвекцией
    noniso_qsm_model_type model_type;

public:
    /// @brief Конструктор класса для идентификации
    /// @param settings Настройка, определяющая выбор параметра идентификации
    /// @param pipe Исходная модель трубопровода
    /// @param times Предпосчитанная временная сетка
    /// @param control_data Предпосчитанные краевые условия
    /// @param etalon_temperature Предпосчитанные эталонные значения
    ident_nonisothermal_qsm_pipe_parameters_t(
        const ident_nonisothermal_qsm_pipe_settings& settings,
        const pipe_noniso_properties_t& pipe, 
        const oil_parameters_t& oil, 
        const std::vector<double>& times, 
        const std::vector<std::vector<double>>& control_data, 
        const std::vector<double>& etalon_temperature, 
        noniso_qsm_model_type model_type)
        : settings(settings)
        , pipe_nominal{ pipe }
        , pipe_to_ident{ pipe }
        , oil{ oil }
        , times{ times }
        , control_data{ control_data }
        , etalon_temperature{ etalon_temperature }
        , model_type(model_type)
    {
        // Проверяем коректность выбора параметра идентификации
        settings.check_parameters();
    }
protected:
    /// @brief Расчёт невязок - проводится моделирование работы ЛУ и считается отклонение расчётного и фактического давления в конце ЛУ 
    /// @param ident_parameters Текущие значения параметров адаптации
    /// @return Вектор невязок - расхождения расчётного и фактического давления на выхое ЛУ
    std::vector<double> calc_vector_residuals(const ident_parameters_nonisothermal_qsm& ident_parameters)
    {
        // Сущность для сбора расчётных данных температуры в конце ЛУ
        nonisothermal_qsm_batch_Tout_collector_t collector(times);
        // Применяем текущие параметры идентификации на модели трубопровода
        ident_parameters.set_parameter(pipe_nominal, &pipe_to_ident);

        // Проводим гидравлический изотермический квазистационарный расчёт
        qsm_noniso_T_task_t task(pipe_to_ident, oil, model_type);
        quasistatic_batch(task, times, control_data, &collector);
        
        const std::vector<double>& calc_temp = collector.get_temp_out_calculated();
        std::vector<double> simulation_result(times.size());
        // Считаем расхождение расчёта и факта
        std::transform(calc_temp.begin(), calc_temp.end(), etalon_temperature.begin(), simulation_result.begin(),
            [](double etalon, double calc) { return etalon - calc; });


        return simulation_result;
    }
public:
    /// @brief Инициализация расчёта вектора невязок
    /// @param parameters Текущие значения параметров идентификации
    /// @return Вектор невязок - расхождения расчётного и фактического давления на выхое ЛУ
    virtual Eigen::VectorXd residuals(const Eigen::VectorXd& parameters) override {

        ident_parameters_nonisothermal_qsm ident_parameters = settings.unpack_ident_parameters(parameters);
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
