﻿#pragma once

#include <chrono>

/// @brief Перевод UNIX времени в строку формата dd:mm:yyyy HH:MM:SS
/// @param t время UNIX
/// @return строку формата dd:mm:yyyy HH:MM:SS
inline std::string UnixToString(time_t t) {

    struct tm tm;
    gmtime_s(&tm, &t);
    char buffer[200];
    strftime(buffer, sizeof(buffer), "%d.%m.%Y %H:%M:%S", &tm);
    return std::string(buffer);
}

/// @brief Перевод строки формата dd:mm:yyyy HH:MM:SS в переменную UNIX времени
/// @param source строка формата dd:mm:yyyy HH:MM:SS
/// @return время UNIX
inline std::time_t StringToUnix(const std::string& source)
{

    using namespace std::chrono;
    using namespace std::literals::string_literals;
    auto in = std::istringstream(source);
    auto tp = sys_seconds{};
    in >> parse("%d.%m.%Y %H:%M:%S"s, tp);
    time_t result = system_clock::to_time_t(tp);
    //string s = UnixToString(result); // для проверки
    return result;
}

/// @brief Перевод из строкого типа в тип double
/// @param str Переменная строкового типа
/// @param delim разделитель целой и дробной частей
/// @return Переменная типа double
inline double str2double(const std::string& str, char delim = '.')
{
    if (delim == '.') {
        std::stringstream ss;
        ss.str(str);
        double result;
        ss >> result;
        return result;
    }
    else {
        size_t i = str.find(delim);
        if (i == std::string::npos) {
            std::stringstream ss;
            ss.str(str);
            double result;
            ss >> result;
            return result;
        }
        else {
            std::string str_ = str;
            str_[i] = '.';
            std::stringstream ss;
            ss.str(str_);
            double result;
            ss >> result;
            return result;
        }

    }
}

/// @brief Деление строки на подстроки по знаку
/// @param str_to_split Изначальная строка
/// @param delimeter Делитель
/// @return Вектор подстрок
inline vector<std::string> split_str(std::string& str_to_split, char delimeter)
{
    std::stringstream str_stream(str_to_split);
    std::vector<std::string> split_output;
    std::string str;
    while (getline(str_stream, str, delimeter))
    {
        split_output.push_back(str);
    }

    return split_output;
}

/// @brief Перевод единиц измерения
class dimension_converter
{
public:
    /// @brief Пересчёт единиц измерений
    /// @param value Текущее значение
    /// @param c Коэффициенты перевода
    /// @return Значение параметра в заданных единицах измерения
    static double convert_dimension(double value, std::pair<double, double> c = std::make_pair(1.0, 0.0))
    {
        return value / c.first - c.second;
    };

    /// @brief Перевод единиц измерения
    /// @param value Текущее значение
    /// @param dimension Инструкция перевода
    /// @return Значение в требуемых единицах измерения
    double convert(double value, const std::string& dimension) const
    {
        auto it = units.find(dimension);
        if (it != units.end()) {
            value = convert_dimension(value, it->second);
        }

        return value;
    }

private:
    /// @brief Коэффициенты для перевода единиц
    const map<std::string, std::pair<double, double>> units{
        {"m3/s", {1.0, 0}},
        {"m3/h-m3/s", {3600, 0.0}},
        {"m3/h", {1 / 3600, 0.0}},
        {"K-C", {1.0, -KELVIN_OFFSET}},
        {"K", {1.0, 0.0}},
        {"kgf/cm2(e)", {1 / TECHNICAL_ATMOSPHERE, -ATMOSPHERIC_PRESSURE / TECHNICAL_ATMOSPHERE}},
        {"kgf/cm2", {1 / TECHNICAL_ATMOSPHERE, -ATMOSPHERIC_PRESSURE / TECHNICAL_ATMOSPHERE}},
        {"kgf/cm2(a)", {1 / TECHNICAL_ATMOSPHERE, 0.0}},
        {"MPa-Pa", {1e-6, 0.0}},
        {"MPa-kPa", {1e-3, 0.0}},
        {"mm^2/s-m^2/s", {1e+6, 0.0}},
    };
};