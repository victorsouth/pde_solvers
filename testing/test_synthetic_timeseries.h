#pragma once
#include <pde_solvers/timeseries.h>

/// @brief Проверка 
TEST(SyntheticTimeSeries, CheckGeneration)
{
    // Задаем настроечные параметры по умолчанию
    timeseries_generator_settings settings = timeseries_generator_settings::default_settings();
    vector<pair<string, double>> timeseries_initial_values = {
        { "Q", 0.3 }, // Задаю расход, м3/с
        { "p_in", 5e6}, // Задаю давление, Па
        { "rho_in", 850 }, // Задаем плотность, кг/м3
        { "visc_in", 17e-6}, // Задаем вязкость, м2/сек
    };
    synthetic_time_series_generator data_generator(timeseries_initial_values, settings);
    // Получаем временные ряды из генератора
    const auto& data = data_generator.get_data();
    // Записываем временной ряд в вектор временных рядов
    vector_timeseries_t params(data);
    // Задаём интересующий нас момент времени
    time_t test_time = static_cast<time_t>(std::time(nullptr) + 200000);
    // Интерополируем значения параметров в заданный момент времени
    vector<double> values_in_test_time = params(test_time);
    // Проверяем корректность генерирования данных
    ASSERT_NEAR(0.3, values_in_test_time[0], settings.value_relative_increment);
}

/// @brief Проверка 
TEST(SyntheticTimeSeries, CheckJump)
{
    // Задаем настроечные параметры по умолчанию
    timeseries_generator_settings settings = timeseries_generator_settings::default_settings();
    vector<pair<string, double>> timeseries_initial_values = {
        { "Q", 0.3 }, // Задаю расход, м3/с
        { "p_in", 5e6}, // Задаю давление, Па
        { "rho_in", 850 }, // Задаем плотность, кг/м3
        { "visc_in", 17e-6}, // Задаем вязкость, м2/сек
    };
    synthetic_time_series_generator data_generator(timeseries_initial_values, settings);
    // Вводим изменения в наши временные ряды
    const time_t jump_time_Q = 150000;
    // На сколько изменится значение нашего расхода
    const double jump_value_Q = -0.1;
    // Применяем скачок
    data_generator.apply_jump(jump_time_Q, jump_value_Q, "Q");
    // Получаем временные ряды из генератора
    const auto& data = data_generator.get_data();
    // Записываем временной ряд в вектор временных рядов
    vector_timeseries_t params(data);
    // Задаём интересующий нас момент времени
    time_t test_time = static_cast<time_t>(std::time(nullptr) + 200000);
    // Интерополируем значения параметров в заданный момент времени
    vector<double> values_in_test_time = params(test_time);
    // Проверяем корректность генерирования данных
    ASSERT_NEAR(0.2, values_in_test_time[0], settings.value_relative_increment);
}

/// @brief Базовый пример генерации синтетических временных рядов
TEST(SyntheticTimeSeries, UseCase_PrepareTimeSeries)
{
    // Задаем настроечные параметры по умолчанию
    timeseries_generator_settings settings = timeseries_generator_settings::default_settings();
    vector<pair<string, double>> timeseries_initial_values = {
        { "Q", 0.3 }, // Задаю расход, м3/с
        { "p_in", 5e6}, // Задаю давление, Па
        { "rho_in", 850 }, // Задаем плотность, кг/м3
        { "visc_in", 17e-6}, // Задаем вязкость, м2/сек
    };
    synthetic_time_series_generator data_generator(timeseries_initial_values, settings);
    // Вводим изменения в наши временные ряды
    const time_t jump_time_Q = 150000;
    // На сколько изменится значение нашего расхода
    const double jump_value_Q = -0.1;
    // Применяем скачок
    data_generator.apply_jump(jump_time_Q, jump_value_Q, "Q");
    // Получаем временные ряды из генератора
    const auto& data = data_generator.get_data();
    // Записываем временной ряд в вектор временных рядов
    vector_timeseries_t params(data);
    // Задаём интересующий нас момент времени
    time_t test_time = static_cast<time_t>(std::time(nullptr) + 200000);
    // Интерополируем значения параметров в заданный момент времени
    vector<double> values_in_test_time = params(test_time);
}