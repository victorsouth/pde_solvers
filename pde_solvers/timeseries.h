#pragma once

#include <fstream>
#include <chrono>
#include <iostream>
#include <random>
#include <algorithm>
#include <iomanip>
#include <ctime>

#include "timeseries/timeseries_helpers.h"
#include "timeseries/csv_readers.h" 
#include "timeseries/vector_timeseries.h" 
#include "timeseries/synthetic_timeseries.h" 


/// @brief Для моделирования работы ЛУ считываем значения параметров с реального трубопровода
/// и на основе этих данных предпосчитываем интерполированные значения 
/// краевых условий (плотность, вязкость, расход и давление в начале ЛУ) 
/// и эталонного параметра (давление на выходе трубопровода)
/// для временной сетки с шагом равным одной минуте
/// @param path_to_real_data Путь к данным с реального трубопровода
/// @return Кортёж, состоящий из 3-х элементов
/// 1. Временная сетка
/// 2. Предпосчитанные краевые условия
/// 3. Предпосчитанные эталонные значения
inline  std::tuple<std::vector<double>, std::vector<std::vector<double>>, std::vector<double>>
prepare_timeseries_data(
    double step,
    const std::string& start_period,
    const std::string& end_period,
    const std::vector<std::pair<std::string, std::string>>& parameters)
{
    using namespace std;
    // Временные ряды краевых условий
    std::vector<std::pair<std::vector<time_t>, std::vector<double>>> control_tag_data;
    // Временные ряды эталонных данных
    std::vector<std::pair<std::vector<time_t>, std::vector<double>>> etalon_tag_data;
    // Задаём период
    // Считываем временные ряды параметров
    csv_multiple_tag_reader tags(parameters);

    control_tag_data = tags.read_csvs(start_period, end_period); //чувствителен к отсутствию явного формата секунд для некоторых csv файлов (по опыту, когда секунды :00)
    etalon_tag_data = { control_tag_data.back() };
    control_tag_data.pop_back();

    // Помещаем временные ряды в вектор
    vector_timeseries_t control_parameters_time_series(control_tag_data);
    vector_timeseries_t etalon_parameters_time_series(etalon_tag_data);


    // Определяем начало периода
    time_t start_period_time = std::max(control_parameters_time_series.get_start_date(), etalon_parameters_time_series.get_start_date());
    // Определяем конец периода
    time_t end_period_time = std::min(control_parameters_time_series.get_end_date(), etalon_parameters_time_series.get_end_date());
    // Определяем продолжительность периода
    time_t duration = (end_period_time - start_period_time);

    // Считаем количество точек в сетке
    size_t dots_count = static_cast<size_t>(ceil(duration / step) + 0.00001);

    std::vector<double>  times = std::vector<double>(dots_count);
    std::vector<double>  times_vector = std::vector<double>(dots_count);
    std::vector<std::vector<double>> control_data = std::vector<std::vector<double>>(dots_count);
    std::vector<double> etalon_temp = std::vector<double>(dots_count);

    for (size_t i = 0; i < dots_count; i++)
    {
        // Заполняем временную сетку
        times[i] = step * i;
        // Определяем момент времени
        time_t t = start_period_time + static_cast<time_t>(times[i] + 0.5);
        times_vector[i] = t;  // Передача временных меток 
        // Получаем интерполированные значения краевых условий и эталонных значений
        control_data[i] = control_parameters_time_series(t);
        etalon_temp[i] = etalon_parameters_time_series(t).front();

    };

    return std::make_tuple(std::move(times_vector), std::move(control_data), std::move(etalon_temp));
}

