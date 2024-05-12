#pragma once

#include <chrono>
using std::vector;
using std::pair;


/// @brief Векторный верменной ряд
class vector_timeseries_t {
private:
    vector<pair<vector<time_t>, vector<double>>> data;
    time_t start_date;
    time_t end_date;
    mutable vector<size_t> left_bound;

public:
    /// @brief Получение количества значений определённого параметра
    /// @param numb Номер параметра 
    size_t get_elements_count(size_t numb) const
    {
        return data[numb].first.size();
    }
    /// @brief Конструктор
    /// @param data Вектор временных рядов, каждый элемент которого
    /// представляет собой пару, в которой первый элемент это временная сетка,
    /// а второй - вектор значений параметров в соответствующие моменты времени
    vector_timeseries_t(const vector<pair<vector<time_t>, vector<double>>>& data)
        : data(data)
    {
        if (data.empty())
            return;

        std::tie(start_date, end_date) = get_timeseries_period(data);

        left_bound = vector<size_t>(data.size(), 0);

    };
    /// @brief Получение времени начала периода 
    time_t get_start_date() const {
        return start_date;
    };
    /// @brief Расчёт Длительности периода
    time_t get_duration() const {
        return end_date - start_date;
    };
    /// @brief Рассчитывает астрономическое время по модельному
    time_t get_astronomic_time(double model_time) const {
        return static_cast<time_t>(model_time + start_date + 0.5);
    };
    /// @brief Возвращает последнее время временных рядов
    time_t get_end_date() const {
        return end_date;
    }

    /// @brief Возвращает интерполированные значения 
    /// временных рядов в момент времени t
    /// @param t Момент времени
    /// @return Интерполированные значения
    vector<double> operator()(time_t t) const
    {
        for (size_t i = 0; i < data.size(); ++i) {
            const auto& times = data[i].first;
            if (times[left_bound[i]] > t) {
                // если один раз сдвинулись правее некоторой точки, то обратно вернуться уже не разрешаем
                throw std::logic_error("wrong time value");
            }
        }

        vector<double> result(data.size());
        for (size_t i = 0; i < data.size(); ++i) {
            const auto& times = data[i].first;
            const auto& values = data[i].second;
            // it - указывает на элемент либо равный t, либо больший 
            auto it = std::lower_bound(times.begin() + left_bound[i], times.end(), t);
            if (it != times.end()) {
                size_t k = it - times.begin();
                // запоминаем левую границу я
                left_bound[i] = k == 0 ? 0 : k - 1;
                if (times[k] == t) {
                    // точное число есть в данных
                    result[i] = values[k];
                }
                else {
                    // точное число режит между values[k-1] и values[k]
                    if (k == 0) {
                        result[i] = std::numeric_limits<double>::quiet_NaN();
                    }
                    else {
                        time_t t_prev = times[k - 1];
                        time_t t_next = times[k];

                        // Если t=t_prev, будет 0, если t=t_next, будет 1
                        double alpha = 1.0 * (t - t_prev) / (t_next - t_prev);

                        double v_prev = values[k - 1];
                        double v_next = values[k];

                        result[i] = (1 - alpha) * v_prev + alpha * v_next;
                    }
                }

            }
            else {
                result[i] = std::numeric_limits<double>::quiet_NaN();
            }
        };


        return result;
    }

private:
    /// @brief Определение начала и конца периода
    /// @param data Временные ряды параметров
    /// @return Начало и конец периода
    static pair<time_t, time_t> get_timeseries_period(const vector<pair<vector<time_t>, vector<double>>>& data)
    {
        time_t start_date = std::numeric_limits<time_t>::min();
        time_t end_date = std::numeric_limits<time_t>::max();;
        for (size_t i = 0; i < data.size(); ++i) {
            if (!data[i].first.empty()) {
                start_date = std::max(start_date, data[i].first.front());
                end_date = std::min(end_date, data[i].first.back());
            }
        }

        if (start_date > end_date) {
            throw std::logic_error("wrong data period");
        }
        return std::make_pair(start_date, end_date);
    }
};