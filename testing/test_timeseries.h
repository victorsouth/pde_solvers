#pragma once

TEST(CsvRead, ReadStream)
{
    stringstream ss;
    ss << "10.08.2021 08:30:50;5" << std::endl;
    ss << "10.08.2021 08:40:50;6" << std::endl;
    ss << "10.08.2021 08:50:50;6.2" << std::endl;

    auto [time, values] = csv_tag_reader::read_from_stream(ss, "MPa-kPa");

    ASSERT_EQ(UnixToString(time[1]), "10.08.2021 08:40:50");
    ASSERT_NEAR(6200.0, values[2], 1);
}

TEST(CsvRead, ReadStreamWithPeriod)
{
    stringstream ss;
    ss << "10.08.2021 08:30:50;5" << std::endl;
    ss << "10.08.2021 08:40:50;6" << std::endl;
    ss << "10.08.2021 08:50:50;6.2" << std::endl;
    ss << "10.08.2021 09:55:50;5.5" << std::endl;

    auto [time, values] = csv_tag_reader::read_from_stream(ss, "MPa-kPa", 
        StringToUnix("10.08.2021 08:35:50;5"), StringToUnix("10.08.2021 09:30:50"));

    ASSERT_EQ(UnixToString(time[0]), "10.08.2021 08:40:50");
    ASSERT_EQ(time.size(), 2);
    ASSERT_NEAR(6200.0, values[1], 1);
}

TEST(VectorTimeseries, Interpolation)
{
    stringstream pres;
    pres << "10.08.2021 08:30:50;5" << std::endl;
    pres << "10.08.2021 08:40:50;6" << std::endl;
    pres << "10.08.2021 08:50:50;6.2" << std::endl;
    pres << "10.08.2021 09:55:50;5.5" << std::endl;

    auto press = csv_tag_reader::read_from_stream(pres, "MPa-kPa");

    vector_timeseries_t timeseries({ press });

    time_t test_time = StringToUnix("10.08.2021 08:45:50");
    vector<double> inter_values = timeseries(test_time);

    ASSERT_NEAR(inter_values[0], 6100, 10);

    time_t wrong_time = StringToUnix("10.08.2021 08:35:50");

    ASSERT_ANY_THROW(timeseries(wrong_time));

}


TEST(Timeseries, UseCase)
{
    string folder = "data/";
    vector<pair<string, string>>parameters =
    {
        { folder + "rho_in", "kg/m3" },
        { folder + "visc_in", "mm^2/s-m^2/s"s },
        { folder + "p_in", "MPa-kPa"s },
        { folder + "Q_in", "m3/h-m3/s"s },
        { folder + "p_out", "MPa-kPa"s },
        { folder + "Q_out", "m3/h-m3/s"s }
    };

    string start_period = "01.08.2021 00:00:00";
    string end_period = "01.09.2021 00:00:00";

    csv_multiple_tag_reader tags(parameters);

    auto tag_data = tags.read_csvs(start_period, end_period);

    vector_timeseries_t params(tag_data);

    time_t test_time = StringToUnix("10.08.2021 08:53:50");

    vector<double> values_in_test_time = params(test_time);
}