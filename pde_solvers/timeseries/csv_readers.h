#pragma once

#include <fstream>
#include "timeseries_helpers.h" 


using std::pair;
using std::string;
using std::vector;


/// @brief Чтение исторических данных по одному тегу
class csv_tag_reader
{
public:
    /// @brief Чтение исторических данных
    /// @param input_stream Входной поток
    /// @param dimension Инструкция для перевода единиц измерения
    /// @param time_begin Начало периода
    /// @param time_end Конец периода
    /// @return Временной ряд
    static pair<vector<time_t>, vector<double>> read_from_stream(std::istream& input_stream, const string& dimension, time_t time_begin = std::numeric_limits<time_t>::min(),
        time_t time_end = std::numeric_limits<time_t>::max())
    {
        vector<time_t> t;
        vector<double> x;

        string line_from_file;
        vector<string> split_line_to_file;
        while (getline(input_stream, line_from_file))
        {
            split_line_to_file = split_str(line_from_file, ';');
            const auto& date = split_line_to_file[0];
            time_t ut = StringToUnix(date);

            if (ut < time_begin) {
                continue;
            }
            if (ut > time_end) {
                break;
            }

            double value = str2double(split_line_to_file[1], ',');

            converter_dimension converter;
            value = converter.convert(value, dimension);

            t.emplace_back(std::move(ut));
            x.emplace_back(std::move(value));
        }

        return std::make_pair(t, x);
    }
private:
    /// @brief Чтение одного файла 
    /// @param filename Название файла
    /// @param dimension Инструкция перевода 
    /// единиц измерения
    /// @param time_begin Начало периода
    /// @param time_end Конец периода
    /// @return Временной ряд
    static pair<vector<time_t>, vector<double>> read_from_file(const string& filename, const string& dimension, time_t time_begin = std::numeric_limits<time_t>::min(),
        time_t time_end = std::numeric_limits<time_t>::max())
    {
        std::ifstream infile(filename);
        if (!infile)
            throw std::runtime_error("file is not exist");

        return read_from_stream(infile, dimension, time_begin, time_end);
    }
public:
    /// @brief Конструктор
    /// @param data Название тега и инструкция перевода единиц измерения
    csv_tag_reader(const pair<string, string>& data)
        :tagname_dim{ data }
    {

    };

    /// @brief Чтение параметра для заданного периода
    /// @param unixtime_from Начало периода задаётся строкой формата dd:mm:yyyy HH:MM:SS
    /// @param unixtime_to Конец периода задаётся строкой формата dd:mm:yyyy HH:MM:SS
    /// @return возвращает временной ряд в формате 
    /// для хранения в vector_timeseries_t
    pair<vector<time_t>, vector<double>> read_csv(const string& unixtime_from, const string& unixtime_to) const
    {
        time_t start_period = StringToUnix(unixtime_from);
        time_t end_period = StringToUnix(unixtime_to);

        return read_csv(start_period, end_period);
    }

    /// @brief Чтение параметров для заданного периода
    /// @param start_period Начало периода задаётся типом данных time_t
    /// @param end_period Конец периода задаётся типом данных time_t
    /// @return возвращает временной ряд в формате 
    /// для хранения в vector_timeseries_t
    pair<vector<time_t>, vector<double>> read_csv(
        time_t start_period = std::numeric_limits<time_t>::min(),
        time_t end_period = std::numeric_limits<time_t>::max()) const
    {
        pair<vector<time_t>, vector<double>> data;

        string extension = ".csv";

        data = read_from_file(tagname_dim.first + extension, tagname_dim.second + extension, start_period, end_period);

        return data;
    };

private:
    /// @brief Название тега и инструкция перевода единиц измерения
    const pair<string, string> tagname_dim;
};

/// @brief Чтение параметров из исторических данных для нескольких тегов
class csv_multiple_tag_reader
{
public:

    /// @brief Конструктор
    /// @param data Вектор пар, в которых первый элемент название тега, 
    /// а второй - инструкция для перевода единиц измерения
    csv_multiple_tag_reader(const vector<pair<string, string>>& data)
        :filename_dim{ data }
    {

    };

    /// @brief Чтение параметров для заданного периода
    /// @param unixtime_from Начало периода задаётся строкой формата dd:mm:yyyy HH:MM:SS
    /// @param unixtime_to Конец периода задаётся строкой формата dd:mm:yyyy HH:MM:SS
    /// @return возвращает временной ряд в формате 
    /// для хранения в vector_timeseries_t
    vector<pair<vector<time_t>, vector<double>>> read_csvs(const string& unixtime_from, const string& unixtime_to) const
    {
        time_t start_period = StringToUnix(unixtime_from);
        time_t end_period = StringToUnix(unixtime_to);

        return read_csvs(start_period, end_period);
    }

    /// @brief Чтение параметров для заданного периода
    /// @param start_period Начало периода задаётся типом данных time_t
    /// @param end_period Конец периода задаётся типом данных time_t
    /// @return возвращает временной ряд в формате 
    /// для хранения в vector_timeseries_t
    vector<pair<vector<time_t>, vector<double>>> read_csvs(
        time_t start_period = std::numeric_limits<time_t>::min(),
        time_t end_period = std::numeric_limits<time_t>::max()) const
    {
        vector<pair<vector<time_t>, vector<double>>> data;

        for (size_t i = 0; i < filename_dim.size(); i++)
        {
            csv_tag_reader tag_reader(filename_dim[i]);
            pair<vector<time_t>, vector<double>> data_from_single_csv =
                tag_reader.read_csv(start_period, end_period);
            data.emplace_back(std::move(data_from_single_csv));
        }

        return data;
    };

private:
    /// @brief Вектор пар, из которых первый элемент название тега,
    /// а второй - инструкция для перевода единиц измерения
    const vector<pair<string, string>> filename_dim;
};