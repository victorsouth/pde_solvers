#pragma once

#include <unordered_map>
#include <filesystem>


template <typename DataTypeOriginal, typename DataTypeConverted>
inline std::string convert_csv_values(const std::string& filename,
    const std::unordered_map<DataTypeOriginal, DataTypeConverted>& conversion_map,
    const DataTypeConverted& exception_value = std::numeric_limits<DataTypeConverted>::max())
{
    using std::filesystem::path;
    
    // Проверка существования исходного файла
    std::ifstream infile(filename + ".csv");
    if (!infile)
        throw std::runtime_error("file is not exist");

    // Формирование имени конвертированного файла в директории с исходным файлом
    std::string converted_filename = filename + "_converted";
    
    if (std::filesystem::exists(converted_filename)) {
        // Временной ряд уже был сконвертирован ранее
        return converted_filename;
    }

    // Чтение значений из исходного файла
    std::pair<std::vector<time_t>, std::vector<DataTypeOriginal>>
        original_timeseries = csv_tag_reader::read_from_stream<DataTypeOriginal>(
            infile, "" /* Безразмерный параметр */);
    
    // Конвертация значений
    const std::vector<std::time_t>& timestamps = original_timeseries.first;
    std::vector<DataTypeConverted> converted_values;

    for (const DataTypeOriginal& original : original_timeseries.second) {
        DataTypeConverted converted_value = conversion_map.at(original);
        if (converted_value == exception_value)
            throw std::runtime_error("exception_value found. Check original csv");
        converted_values.emplace_back(converted_value);
    }

    // Запись значений
    write_csv_tag_file<DataTypeConverted>(converted_filename, timestamps, converted_values);

    return converted_filename;
}
