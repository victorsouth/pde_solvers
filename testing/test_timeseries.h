#pragma once
#include <pde_solvers/timeseries.h>
/// @brief  Проверяет способность чтения данных из потока
TEST(CsvRead, ReadStream)
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

/// @brief Проверяет способность чтения данных из потока за заданный период
TEST(CsvRead, ReadStreamWithPeriod)
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
