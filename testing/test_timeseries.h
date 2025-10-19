#pragma once
#include <pde_solvers/timeseries.h>

/// @brief Проверяет функционал конвертации значений временных рядов в csv формате 
TEST(ConvertCsvValues, HandelsIntToIntConversion)
{
    // Исходные данные
    const std::string test_timestamp = "10.08.2021 08:30:50";
    const int original_value = 10;
    const int original_mapped = 2;

    // Формируем исходных csv файл
    std::string test_folder = prepare_test_folder();
    std::string original_filename = test_folder + "original";
   
    write_csv_tag_file<int>(original_filename,
        { StringToUnix(test_timestamp) }, // временные метки
        { original_value } // значения
    ); 

    // Формируем сопоставление
    const std::unordered_map<int, int> conversion_map = {
        {original_value, original_mapped}
    };
    
    // Конвертация файла
    std::string new_filename = convert_csv_values(original_filename, conversion_map);

    // Зачитываем преобразованный файл
    auto reader = csv_tag_reader(new_filename, "" /* Безразмерный параметр */);
    auto [new_times, new_values] = reader.read_csv<int>();

    ASSERT_EQ(UnixToString(new_times.front()), test_timestamp);
    ASSERT_EQ(new_values.front(), conversion_map.at(original_value));
}

/// @brief  Проверяет способность чтения данных в формате с плавающей точкой из потока 
TEST(CsvRead, ReadStreamDouble)
{
    // Записываем в поток данные параметра
    std::stringstream ss;
    ss << "10.08.2021 08:30:50;5" << std::endl;
    ss << "10.08.2021 08:40:50;6" << std::endl;
    ss << "10.08.2021 08:50:50;6.2" << std::endl;

    // Читаем данных из потока
    auto [time, values] = csv_tag_reader::read_from_stream(ss, "MPa");

    // Проверяем корректность считанных данных
    ASSERT_EQ(UnixToString(time[1]), "10.08.2021 08:40:50");
    ASSERT_NEAR(6200000.0, values[2], 1);
}

/// @brief  Проверяет способность чтения целочисленных данных из потока 
TEST(CsvRead, ReadStreamInt)
{
    // Записываем в поток данные параметра
    std::stringstream ss;
    ss << "10.08.2021 08:30:50;1" << std::endl;
    ss << "10.08.2021 08:40:50;-1" << std::endl;
    ss << "10.08.2021 08:50:50;2" << std::endl;

    // Читаем данных из потока
    auto [time, values] = csv_tag_reader::read_from_stream<int>(ss, "");

    // Проверяем корректность считанных данных
    ASSERT_EQ(UnixToString(time[1]), "10.08.2021 08:40:50");
    ASSERT_EQ(-1, values[1]);
}

/// @brief  Проверяет способность чтения данных в формате bool из потока 
TEST(CsvRead, ReadStreamBool)
{
    // Записываем в поток данные параметра
    std::stringstream ss;
    ss << "10.08.2021 08:30:50;1" << std::endl;
    ss << "10.08.2021 08:40:50;0" << std::endl;

    // Читаем данных из потока
    auto [time, values] = csv_tag_reader::read_from_stream<bool>(ss, "");

    // Проверяем корректность считанных данных
    ASSERT_EQ(UnixToString(time[1]), "10.08.2021 08:40:50");
    ASSERT_TRUE(values[0]);
    ASSERT_FALSE(values[1]);
}


/// @brief  Проверяет наличие сключения при попытке обработать как bool значение, отличное от 0 или 1 
TEST(CsvRead, ReadStreamBoolThrowsExceptionIfNotZeroOrOne)
{
    // Записываем в поток данные параметра
    std::stringstream ss;
    ss << "10.08.2021 08:40:50;2" << std::endl;

    EXPECT_ANY_THROW(csv_tag_reader::read_from_stream<bool>(ss, ""));
}
/// @brief  Проверяет наличие исключения при попытке чтения данных неподдерживаемого типа
TEST(CsvRead, ReadStreamThrowsExceptionOnUnsupportedType)
{
    // Тип данных, неподдерживаемый в csv_tag_reader::read_from_stream
    struct test_type {};
    
    // Записываем в поток любые данные
    std::stringstream ss;
    ss << "10.08.2021 08:40:50;abc" << std::endl;

    // Проверяем наличие исключения при попытке обработать неизвестный тип данных
    EXPECT_ANY_THROW(csv_tag_reader::read_from_stream<test_type>(ss, ""));
}

/// @brief Проверяет способность чтения данных в формате с плавающей точкой из потока за заданный период
TEST(CsvRead, ReadStreamWithPeriodDouble)
{
    // Записываем в поток данные параметра
    std::stringstream ss;
    ss << "10.08.2021 08:30:50;5" << std::endl;
    ss << "10.08.2021 08:40:50;6" << std::endl;
    ss << "10.08.2021 08:50:50;6.2" << std::endl;
    ss << "10.08.2021 09:55:50;5.5" << std::endl;

    // Читаем данные из потока для заданного периода 
    auto [time, values] = csv_tag_reader::read_from_stream(ss, "MPa", 
        StringToUnix("10.08.2021 08:35:50;5"), StringToUnix("10.08.2021 09:30:50"));

    // Проверяем корректность считанных данных
    ASSERT_EQ(UnixToString(time[0]), "10.08.2021 08:40:50");
    ASSERT_EQ(time.size(), 2);
    ASSERT_NEAR(6200000, values[1], 1);
}

/// @brief Проверка функции интерполяции временных рядов 
TEST(VectorTimeseries, InterpolateTimeseries)
{
    // Записываем в поток данные параметра
    std::stringstream pres;
    pres << "10.08.2021 08:30:50;5" << std::endl;
    pres << "10.08.2021 08:40:50;6" << std::endl;
    pres << "10.08.2021 08:50:50;6.2" << std::endl;
    pres << "10.08.2021 09:55:50;5.5" << std::endl;

    // Читаем данные из потока
    auto press = csv_tag_reader::read_from_stream(pres, "MPa");

    // Записываем временной ряд в вектор временных рядов
    vector_timeseries_t timeseries({ press });

    // Интерполируем значение параметра в заданный момент времени
    time_t test_time = StringToUnix("10.08.2021 08:45:50");
    std::vector<double> inter_values = timeseries(test_time);

    // Проверяем соответствие полученного значения
    ASSERT_NEAR(inter_values[0], 6100000, 10);

}

/// @brief Проверка функционала по выбору способа интерполяции
TEST(VectorTimeseries, CheckInterpolationMethods)
{
    // Записываем в поток данные параметра
    std::stringstream pres;
    pres << "10.08.2021 08:30:50;5" << std::endl;
    pres << "10.08.2021 08:40:50;6" << std::endl;
    pres << "10.08.2021 08:50:50;6.2" << std::endl;
    pres << "10.08.2021 09:55:50;5.5" << std::endl;

    // Читаем данные из потока для заданного периода
    auto press = csv_tag_reader::read_from_stream(pres, "MPa");

    // Заданный момент времени
    time_t test_time = StringToUnix("10.08.2021 08:45:50");

    // Записываем временной ряд в вектор временных рядов и указываем линейный способ интерполяции
    vector_timeseries_t timeseries_linear({ press }, InterplationMethod::Linear);

    // Интерполируем значение параметра в заданный момент времени
    std::vector<double> inter_values = timeseries_linear(test_time);

    // Проверяем соответствие полученного значения при линейной интерполяции
    ASSERT_NEAR(inter_values[0], 6100000, 10);

    // Записываем временной ряд в вектор временных рядов и указываем ступенчатый способ интерполяции
    vector_timeseries_t timeseries_step({ press }, InterplationMethod::Step);

    // Интерполируем значение параметра в заданный момент времени
    inter_values = timeseries_step(test_time);

    // Проверяем соответствие полученного значения при ступенчатой интерполяции
    ASSERT_NEAR(inter_values[0], 6000000, 10);
}

/// @brief Проверка обработки случая, когда момент времени
/// смещается левее точки, в которой уже происходила интерполяция
TEST(VectorTimeseries, CheckWrongTime)
{
    // Записываем в поток данные параметра
    std::stringstream pres;
    pres << "10.08.2021 08:30:50;5" << std::endl;
    pres << "10.08.2021 08:40:50;6" << std::endl;
    pres << "10.08.2021 08:50:50;6.2" << std::endl;
    pres << "10.08.2021 09:55:50;5.5" << std::endl;

    // Читаем данные из потока для заданного периода
    auto press = csv_tag_reader::read_from_stream(pres, "MPa");

    // Записываем временной ряд в вектор временных рядов
    vector_timeseries_t timeseries({ press });

    // Интерполируем значение параметра в заданный момент времени
    time_t test_time = StringToUnix("10.08.2021 08:45:50");
    std::vector<double> inter_values = timeseries(test_time);

    // Выбираем момент времени левее точки предыдущей интерполяции
    time_t wrong_time = StringToUnix("10.08.2021 08:35:50");

    // Проверяем наличие исключения
    ASSERT_ANY_THROW(timeseries(wrong_time));
}

/// @brief Пример использование библиотеки timeseries.h 
TEST(Timeseries, DISABLED_UseCase)
{
    using namespace std::string_literals;

    // Записываем пути к историческим данным
    std::string folder = "data/";
    std::vector<std::pair<std::string, std::string>>parameters =
    {
        { folder + "rho_in", "kg/m3"s },
        { folder + "visc_in", "mm^2/s-m^2/s"s },
        { folder + "p_in", "MPa"s },
        { folder + "Q_in", "m3/h-m3/s"s },
        { folder + "p_out", "MPa"s },
        { folder + "Q_out", "m3/h-m3/s"s }
    };

    // Задаём период
    std::string start_period = "01.08.2021 00:00:00";
    std::string end_period = "01.09.2021 00:00:00";

    // Считываем временные ряды параметров
    csv_multiple_tag_reader tags(parameters);
    auto tag_data = tags.read_csvs(start_period, end_period);

    // Помещаем временные ряды в вектор
    vector_timeseries_t params(tag_data);

    // Задаём интересующий нас момент времени
    time_t test_time = StringToUnix("10.08.2021 08:53:50");

    // Интерополируем значения параметров в заданный момент времени
    std::vector<double> values_in_test_time = params(test_time);
}
