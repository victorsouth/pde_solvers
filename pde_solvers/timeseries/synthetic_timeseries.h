#pragma once

#include <iostream>
#include <random>
#include <algorithm>
#include <iomanip>
#include "timeseries_helpers.h" 

using std::vector;
using std::pair;
using std::string;
using std::time_t;

/// @brief Исходны данные и настроечные параметры
struct timeseries_generator_settings {
    /// @brief Время начала моделирования, с
    std::time_t start_time;
    /// @brief Время моделирования, с
    std::time_t duration;
    /// @brief Минимальное значение размаха шага, с
    std::time_t sample_time_min;
    /// @brief Максимальное значение размаха шага, с
    std::time_t sample_time_max;
    /// @brief Относительное минимальное отклонение значения параметров, доли
    double value_relative_decrement;
    /// @brief Относительное максимальное отклонение значения параметров, доли
    double value_relative_increment;
    /// @brief Настроечные параметры по умолчанию 
    static timeseries_generator_settings default_settings() {
        timeseries_generator_settings result;
        result.start_time = std::time(nullptr);
        result.duration = 250000;
        result.sample_time_min = 200;
        result.sample_time_max = 400;
        result.value_relative_decrement = 0.0002;
        result.value_relative_increment = 0.0002;
        return result;
    }
};
/// @brief Класс для генерации синтетических временных рядов
class synthetic_time_series_generator {
    // Типы синонимов для улучшения читаемости кода
    using TimeVector = vector<time_t>;
    using ParamVector = vector<double>;
    using ParamPair = pair<TimeVector, ParamVector>;
    std::random_device rd; // Генератор случайных чисел
    std::mt19937 gen; // Генератор псевдослучайных чисел
public:
    /// @brief Конструктор класса, в котором синтезируются врмененные ряды, похожие на реальные данные (в зависимости от того, какие настройки)
    /// @param initial_values Исходные данные для создания временных рядов
    /// @param settings Настройки генератора временных рядов
    synthetic_time_series_generator(const vector<pair<string, double>> initial_values, const timeseries_generator_settings& settings)
        : initial_values(initial_values),
        settings(settings) {
        std::uniform_real_distribution<double> timeDis(settings.sample_time_min, settings.sample_time_max);

        for (const auto& param : initial_values) {
            TimeVector timeValues;
            ParamVector paramValues;

            std::uniform_real_distribution<double> normalDis(param.second * (1 - settings.value_relative_increment), param.second * (1 + settings.value_relative_decrement));

            time_t timeStep = timeDis(gen);
            for (time_t time = settings.start_time; time <= settings.start_time + settings.duration; time += timeStep) {
                timeValues.push_back(time); // Создается ряд времени для определенного параметра
                timeStep = timeDis(gen);
            }

            std::transform(timeValues.begin(), timeValues.end(), std::back_inserter(paramValues),
                [&](time_t) { return normalDis(gen); });

            data.push_back({ timeValues, paramValues });
        }
    }
    /// @brief Применение скачка к временному ряду (резкое изменение выбираемого параметра на определенное значение)
    /// @param jump_time Время, когда происходит скачок
    /// @param jump_value Значение скачка
    /// @param paramName Имя параметра, к которому применяется скачок
    void apply_jump(time_t jump_time, double jump_value, const string& paramName) {
        for (size_t i = 0; i < initial_values.size(); ++i) {
            if (initial_values[i].first == paramName) {
                auto it = std::lower_bound(data[i].first.begin(), data[i].first.end(), settings.start_time + jump_time);
                size_t position = std::distance(data[i].first.begin(), it);
                // Также используется normalDis для имитации показаний датчиков
                std::uniform_real_distribution<double> normalDis(jump_value * (1 - settings.value_relative_increment), jump_value * (1 + settings.value_relative_increment));
                for (size_t j = position; j < data[i].first.size(); ++j) {
                    double value = normalDis(gen);
                    data[i].second[j] = value;
                }
            }
        }
    }
    /// @brief Получение сгенерированных данных
    /// @return Вектор временных рядов
    const vector<ParamPair>& get_data() const {
        return data;
    }

private:
    /// @brief Исходные данные, обязательно должны присутствовать два опциональных параметра
    const vector<pair<string, double>> initial_values;
    /// @brief Настройки генератора
    timeseries_generator_settings settings;
    /// @brief Данные временных рядов
    vector<ParamPair> data;
};