#pragma once

#include <iostream>
#include <random>
#include <algorithm>
#include <iomanip>
#include "timeseries_helpers.h" 

using std::vector;
using std::pair;
using std::string;
using std::time_t;


// ���� ��������� ��� ��������� ���������� ����
using TimeVector = vector<time_t>;
using ParamVector = vector<double>;
using ParamPair = pair<TimeVector, ParamVector>;

/// @brief ������� ������ � ����������� ���������
struct timeseries_generator_settings {
    /// @brief ����� ������ �������������, �
    std::time_t start_time;
    /// @brief ����� �������������, �
    std::time_t duration;
    /// @brief ����������� �������� ������� ����, �
    std::time_t sample_time_min;
    /// @brief ������������ �������� ������� ����, �
    std::time_t sample_time_max;
    /// @brief ������������� ����������� ���������� �������� ����������, ����
    double value_relative_decrement;
    /// @brief ������������� ������������ ���������� �������� ����������, ����
    double value_relative_increment;
    /// @brief ����������� ��������� �� ��������� 
    static timeseries_generator_settings default_settings() {
        timeseries_generator_settings result;
        result.start_time = std::time(nullptr);
        result.duration = 250000;
        result.sample_time_min = 200;
        result.sample_time_max = 400;
        result.value_relative_decrement = 0.0002;
        result.value_relative_increment = 0.0002;
        return result;
    }

};
/// @brief ����� ��� ��������� ������������� ��������� �����
class synthetic_time_series_generator {
public:
    /// @brief ����������� ������
    /// @param settings ��������� ���������� ��������� �����
    synthetic_time_series_generator(const vector<pair<string, double>> initial_values, const timeseries_generator_settings& settings)
        : initial_values_(initial_values),
        settings_(settings),
        gen(rd()) {
        std::uniform_real_distribution<double> timeDis(settings_.sample_time_min, settings_.sample_time_max);

        for (const auto& param : initial_values_) {
            TimeVector timeValues;
            ParamVector paramValues;

            std::uniform_real_distribution<double> normalDis(param.second * (1 - settings_.value_relative_increment), param.second * (1 + settings_.value_relative_decrement));

            time_t timeStep = timeDis(gen);
            for (time_t time = settings_.start_time; time <= settings_.start_time + settings_.duration; time += timeStep) {
                timeValues.push_back(time);
                timeStep = timeDis(gen);
            }

            std::transform(timeValues.begin(), timeValues.end(), std::back_inserter(paramValues),
                [&](time_t) { return normalDis(gen); });

            data.push_back({ timeValues, paramValues });
        }
    }
    /// @brief ���������� ������ � ���������� ����
    /// @param jump_time �����, ����� ���������� ������
    /// @param jump_value �������� ������
    /// @param paramName ��� ���������, � �������� ����������� ������
    void apply_jump(time_t jump_time, double jump_value, const string& paramName) {
        for (size_t i = 0; i < initial_values_.size(); ++i) {
            if (initial_values_[i].first == paramName) {
                auto it = std::lower_bound(data[i].first.begin(), data[i].first.end(), settings_.start_time + jump_time);
                size_t position = std::distance(data[i].first.begin(), it);
                std::uniform_real_distribution<double> normalDis(jump_value * (1 - settings_.value_relative_increment), jump_value * (1 + settings_.value_relative_increment));
                for (size_t j = position; j < data[i].first.size(); ++j) {
                    double value = normalDis(gen);
                    data[i].second[j] = value;
                }
            }
        }
    }
    /// @brief ��������� ��������������� ������
    /// @return ������ ��������� �����
    const vector<ParamPair>& get_data() const {
        return data;
    }

private:
    /// @brief �������� ������, ����������� ������ �������������� ��� ������������ ���������
    const vector<pair<string, double>> initial_values_;
    /// @brief ��������� ����������
    timeseries_generator_settings settings_;
    /// @brief ������ ��������� �����
    vector<ParamPair> data;
    /// @brief ��������� ��������� �����
    std::random_device rd;
    /// @brief ��������� ��������������� �����
    std::mt19937 gen;
};