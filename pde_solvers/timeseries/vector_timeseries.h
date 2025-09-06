#pragma once


using std::vector;
using std::pair;


/// @brief Способы интерполяции временных рядов
/// Step - Ступенчатая интерполяция
/// Linear - линейная интерполяция
enum class InterplationMethod {
    Step, Linear
};


/// @brief Векторный верменной ряд
template <typename DataType = double>
class vector_timeseries_t {
private:
    /// @brief Метод интерполяции
    InterplationMethod interpolation_method;
    /// @brief Самое позднее начальное время среди временных рядов
    time_t start_date{ std::numeric_limits<time_t>::min() };
    /// @brief Самое ранее конечное время среди временных рядов
    time_t end_date{ std::numeric_limits<time_t>::max() };
    /// @brief Начальные точки индексов временных рядов, создающие левую границу при поиске во времени
    mutable std::vector<size_t> left_bound;
    /// @brief Исходные временные ряды
    std::vector<std::pair<std::vector<time_t>, std::vector<DataType>>> data;

public:
    /// @brief Геттер для хранящихся данных
    const std::vector<std::pair<std::vector<time_t>, std::vector<DataType>>>& get_data() const {
        return data;
    }
    /// @brief Проверяет, что по какой-то причине никаких данных нет
    bool empty() const {
        return data.empty();
    }
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
    vector_timeseries_t(const std::vector<std::pair<std::vector<time_t>, std::vector<DataType>>>& data, InterplationMethod interpolation_method = InterplationMethod::Linear)
        : data(data)
        , interpolation_method(interpolation_method)
    {
        if (data.empty())
            return;
            
        std::tie(start_date, end_date) = get_timeseries_period(data);

        left_bound = std::vector<size_t>(data.size(), 0);

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
    std::vector<DataType> operator()(time_t t) const
    {
        for (size_t i = 0; i < data.size(); ++i) {
            const auto& times = data[i].first;
            if (times[left_bound[i]] > t) {
                // если один раз сдвинулись правее некоторой точки, то обратно вернуться уже не разрешаем
                throw std::logic_error("wrong time value");
            }
        }

        std::vector<DataType> result(data.size());
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
                    // точное число лежит между values[k-1] и values[k]
                    if (k == 0) {
                        result[i] = std::numeric_limits<DataType>::quiet_NaN();
                    }
                    else if (interpolation_method == InterplationMethod::Linear){
                        time_t t_prev = times[k - 1];
                        time_t t_next = times[k];

                        // Если t=t_prev, будет 0, если t=t_next, будет 1
                        double alpha = 1.0 * (t - t_prev) / (t_next - t_prev);

                        DataType v_prev = values[k - 1];
                        DataType v_next = values[k];

                        result[i] = (1 - alpha) * v_prev + alpha * v_next;
                    }
                    else if (interpolation_method == InterplationMethod::Step) {
                        result[i] = values[k - 1];
                    }
                    else {
                        throw std::logic_error("wrong interpolation method");
                    }
                }

            }
            else {
                result[i] = std::numeric_limits<DataType>::quiet_NaN();
            }
        };


        return result;
    }

private:
    /// @brief Определение начала и конца периода
    /// @param data Временные ряды параметров
    /// @return Начало и конец периода
    static std::pair<time_t, time_t> get_timeseries_period(const std::vector<std::pair<std::vector<time_t>, std::vector<DataType>>>& data)
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