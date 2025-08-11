#pragma once



/// @brief Статус величины, привязанной на расчетную схему
enum class measurement_status_t {
    measured,
    inconfident_calculation,
    confident_calculation,
    unused
};

/// @brief Параметры измерения (без указания физической величины)
struct measurement_untyped_data_t {
    /// @brief Внешний идентификатор измерения на стороне БД
    long long id{ -1 };
    /// @brief Значение измерения
    double value{ std::numeric_limits<double>::quiet_NaN() };
    /// @brief Достоверность измерения (пока используется, если измерение берется из расчета)
    measurement_status_t status{ measurement_status_t::measured };
    /// @brief Проверка наличия измерения одного параметра
    bool has_value() const
    {
        return std::isfinite(value);
    }
};

/// @brief Типизированные значение временного ряда
template <typename MeasurementEnumClass>
struct measurement_data_t : public measurement_untyped_data_t {
    /// @brief Тип измерения
    MeasurementEnumClass type{ MeasurementEnumClass(0) };
};


/// @brief Информация о теге временного ряда (имя, размерность)
struct tag_info_t {
    /// @brief Имя тега
    std::string name;
    /// @brief Размерность тега
    std::string dimension;
};

/// @brief Чтение исторических данных по одному тегу
class csv_tag_reader {
public:
    /// @brief Чтение исторических данных
    /// @param input_stream Входной поток
    /// @param dimension Инструкция для перевода единиц измерения
    /// @param time_begin Начало периода
    /// @param time_end Конец периода
    /// @return Временной ряд в формате пары векторов: [Метки времени; Значения параметра]
    static std::pair<std::vector<time_t>, std::vector<double>> read_from_stream(std::istream& input_stream, const std::string& dimension, time_t time_begin = std::numeric_limits<time_t>::min(),
        time_t time_end = std::numeric_limits<time_t>::max())
    {
        std::vector<time_t> t;
        std::vector<double> x;

        std::string line_from_file;
        std::vector<std::string> split_line_to_file;
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

            dimension_converter converter;
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
    /// @return Временной ряд в формате пары векторов: [Метки времени; Значения параметра]
    static std::pair<std::vector<time_t>, std::vector<double>> read_from_file(const std::string& filename, const std::string& dimension, time_t time_begin = std::numeric_limits<time_t>::min(),
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
    csv_tag_reader(const std::pair<std::string, std::string>& data)
        : csv_tag_reader(data.first, data.second)
    { }
    /// @brief Конструтор
    /// @param tag Имя тега
    /// @param dimension Инструкция для перевода единиц измерения
    csv_tag_reader(const std::string& tag, const std::string& dimension)
        : tagname{ tag }, dim{ dimension }
    { }
    /// @brief Конструктор
    csv_tag_reader(const tag_info_t& tag)
        : csv_tag_reader(tag.name, tag.dimension)
    { }

    /// @brief Чтение параметра для заданного периода
    /// @param unixtime_from Начало периода задаётся строкой формата dd:mm:yyyy HH:MM:SS
    /// @param unixtime_to Конец периода задаётся строкой формата dd:mm:yyyy HH:MM:SS
    /// @return возвращает временной ряд в формате 
    /// для хранения в vector_timeseries_t
    std::pair<std::vector<time_t>, std::vector<double>> read_csv(const std::string& unixtime_from, const std::string& unixtime_to) const
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
    std::pair<std::vector<time_t>, std::vector<double>> read_csv(
        time_t start_period = std::numeric_limits<time_t>::min(),
        time_t end_period = std::numeric_limits<time_t>::max()) const
    {
        std::pair<std::vector<time_t>, std::vector<double>> data;

        std::string extension = ".csv";

        data = read_from_file(tagname + extension, dim, start_period, end_period);

        return data;
    };

private:
    /// @brief Название тега и инструкция перевода единиц измерения
    const std::string tagname;
    const std::string dim;
};

/// @brief Чтение параметров из исторических данных для нескольких тегов
class csv_multiple_tag_reader
{
private:
    /// @brief Вектор пар, из которых первый элемент название тега,
    /// а второй - инструкция для перевода единиц измерения
    std::vector<tag_info_t> filename_dim;

public:

    /// @brief Конструктор
    /// @param data Вектор пар, в которых первый элемент название тега, 
    /// а второй - инструкция для перевода единиц измерения
    csv_multiple_tag_reader(const std::vector<std::pair<std::string, std::string>>& tag_list)
    {
        for (const std::pair<std::string, std::string>& tag : tag_list)
        {
            tag_info_t tag_info;
            tag_info.name = tag.first;
            tag_info.dimension = tag.second;
            filename_dim.push_back(tag_info);
        }
    }

    /// @brief Конструктор по списку тегов
    csv_multiple_tag_reader(const std::vector<tag_info_t>& tag_list)
        : filename_dim(tag_list)
    { }

    /// @brief Конструктор по списку тегов
    template <typename TagInfo>
    csv_multiple_tag_reader(const std::vector<TagInfo>& tag_list, const std::string& path = "")
    {
        for (const TagInfo& tag : tag_list)
        {
            tag_info_t tag_info;
            tag_info.name = path + "/" + tag.name;
            tag_info.dimension = tag.dimension;
            filename_dim.push_back(tag_info);
        }
    }
    /// @brief Чтение параметров для заданного периода
    /// @param unixtime_from Начало периода задаётся строкой формата dd:mm:yyyy HH:MM:SS
    /// @param unixtime_to Конец периода задаётся строкой формата dd:mm:yyyy HH:MM:SS
    /// @return возвращает временной ряд в формате 
    /// для хранения в vector_timeseries_t
    std::vector<std::pair<std::vector<time_t>, std::vector<double>>> read_csvs(const std::string& unixtime_from, const std::string& unixtime_to) const
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
    std::vector<std::pair<std::vector<time_t>, std::vector<double>>> read_csvs(
        time_t start_period = std::numeric_limits<time_t>::min(),
        time_t end_period = std::numeric_limits<time_t>::max()) const
    {
        std::vector<std::pair<std::vector<time_t>, std::vector<double>>> data(filename_dim.size());

#pragma omp parallel for schedule(dynamic)
        for (int i = 0; i < static_cast<int>(filename_dim.size()); i++)
        {
            csv_tag_reader tag_reader(filename_dim[i]);
            std::pair<std::vector<time_t>, std::vector<double>> data_from_single_csv =
                tag_reader.read_csv(start_period, end_period);
            data[i] = std::move(data_from_single_csv);
        }

        return data;
    };

};


/// @brief Записывает в формате, пригодном для csv_reader
inline void write_csv_tag_file(const std::string& filename,
    const std::vector<std::time_t>& astro_times,
    const std::vector<double>& values)
{
    std::ofstream file(filename + ".csv");
    for (std::size_t index = 0; index < astro_times.size(); ++index) {
        file << UnixToString(astro_times[index]) << ";" << values[index] << std::endl;
    }

}


/// @brief Запись временного ряда со статусом в файл
/// @param filename Имя .csv файла для записи
/// @param astro_times Временные метки в астрономическом времени
/// @param timeseries Вектор значений
inline void write_csv_tag_file(const std::string& filename,
    const std::vector<std::time_t>& astro_times,
    const std::vector<measurement_untyped_data_t>& timeseries)
{
    std::ofstream file(filename + ".csv");
    for (std::size_t index = 0; index < astro_times.size(); ++index) {
        double value = timeseries[index].value;
        int status = (int)timeseries[index].status;

        file << UnixToString(
            astro_times[index]) << ";" <<
            value << ";" <<
            status << std::endl;
    }

}


/// @brief Записывает в формате, пригодном для csv_reader
/// @param filename Имя файла
/// @param astro_times Моменты времени
/// @param getter Геттер для получения Значения 
/// (передается индекс моментов времени и само значение времени)
template <typename ValueGetter>
inline void write_csv_tag_file(const std::string& filename,
    const std::vector<std::time_t>& astro_times,
    ValueGetter& getter)
{
    std::ofstream file(filename + ".csv");
    for (std::size_t index = 0; index < astro_times.size(); ++index) {
        double value = getter(index, astro_times[index]);

        file << UnixToString(astro_times[index]) << ";" << value << std::endl;
    }

}
