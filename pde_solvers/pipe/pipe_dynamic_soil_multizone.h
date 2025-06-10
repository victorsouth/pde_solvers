#pragma once

namespace pde_solvers {
;

/*
    Структуры данных и расчеты для однозонной модели.
    Реализована сложная система с подготовкой коэффициентов:
      исходные параметры, параметры эквивалентной модели, коэффициенты модели 
    Это нужно для идентификации. Коэффициенты модели зависят от эквивалентных параметров и коэффициентов идентификации   

*/

enum class HeatModelVer { Legacy, V2 };


/// @brief Параметры идентификации тепловой зоны эквивалентной модели 
/// для многозонной модели
struct heat_zone_adaptation_t {
    /// @brief Адаптация теплопроводности слоя изоляции
    double conductivity_isolation{ 1.0 };
    /// @brief Адаптация теплопроводности слоя грунта
    double conductivity_soil{ 1.0 };
    /// @brief Адаптация удельной объемная теплоемкость грунта
    double volume_heat_capacity_soil{ 1.0 };
    /// @brief Адаптация температуры термически невозмущенного грунта
    double temperature_ambient_soil{ 1.0 };
};


/// @brief Коэффициенты двухлойной тепловой модели, зависящие от эквивалентных параметров
struct equivalent_heat_zone_coefficients_t {
    double al1; // alpha1
    double al2; // alpha2
    double b;

    double A1;
    double A2;
    double A3;
    double A4;

    double A;
    double B;
    double C;
    double F;

    /// @brief Температура термически невозмущенного грунта
    /// Поскольку коэффициенты в этой структуре вычисляются с учетом параметров адаптации, 
    /// здесь температура грунта тоже с учетом идентификации
    double temperature_ambient{ KELVIN_OFFSET - 5 };

    /// @brief Время релаксации T грунта
    double get_relaxation_time() const {
        return -1.0 / F;
    }
    /// @brief Коэффициент усиления влияния T жидкости на T грунта
    double get_gain() const {
        return -C / F;
    }
    /// @brief Тепловой поток от термически ВОЗМУЩЕННОГО грунта в трубопровод
    /// @param T_oil Температура жидкости в трубопроводе
    /// @param T_soil Температура термически возмущенного грунта
    /// @param ambient_temperature Температура термически невозмущенного грунта
    double get_heat_flow(double T_oil, double T_soil, double ambient_temperature) const {
        double DT = T_oil - ambient_temperature;
        double DTsr = T_soil - ambient_temperature;
        double q = (A + B * C) * DT + F * B * DTsr;
        return q;
        double q_st = A * DT; // для сравнения
    }
    /// @brief Тепловой поток от термически НЕВОЗМУЩЕННОГО грунта в трубопровод
    /// @param T_oil 
    /// @param ambient_temperature 
    /// @return 
    double get_heat_flow(double T_oil, double ambient_temperature) const {
        double DT = T_oil - ambient_temperature;
        double q = A * DT;
        return q;
    }

    /// @brief Расчет средней темперутры грунта после расчет профиля температуры жидкости по Шухову
    /// @param T_oil Температура жидкости, рассчитанной по Шухову
    /// @param ambient_temperature Температура наружного грунта
    double get_soil_temperature_static(double T_oil, double ambient_temperature) const {
        //Расчет коэффициента связи средней температуры грунта и температуры флюида
        //double K1 = lambda * (pow(al2, 2) - 2 * log(al2) - 1) / (2 * b * log(al1) * log(al2) * (pow(al2, 2) - 1));
        double K = -C / F;

        // Приращение температуры жидкости над Тнар
        double DT1 = T_oil - ambient_temperature;
        double result = K * DT1 + ambient_temperature;

        return result;
    }

    /// @brief Расчет температуры внешнего слоя изоляции
    /// @param T_oil Температура жидкости
    /// @param T_soil Средняя температура грунта
    /// @param ambient_temperature Температура наружного грунта
    double get_isolation_outer_temperature(double T_oil, double T_soil, double ambient_temperature) const {
        double DT1 = T_oil - ambient_temperature;
        double DTsr = T_soil - ambient_temperature;
        double result = (A3 + A4 * C) * DT1 + A4 * F * DTsr + ambient_temperature;
        return result;
    }
    /// @brief Расчет средней температуры грунта на новом слое
    /// @param dt Временной шаг
    /// @param T_oil Температура жидкости в трубе
    /// @param T_soil_old Старая температура грунта
    /// @param ambient_temperature Температура наружного грунта
    /// @return Средня температура грунта на новом слое
    double get_next_soil_temperature(double dt, double T_oil, double T_soil_old, double ambient_temperature) const
    {
        double DT1 = T_oil - ambient_temperature;
        double DTsr = T_soil_old - ambient_temperature;
        double dT = dt * (C * DT1 + F * DTsr);
        double T_soil_new = T_soil_old + dT;
        return T_soil_new;
    }
};



/// @brief Теплофизические и геометрические параметры 
/// эквивалентной двухслойной тепловой модели
struct equivalent_heat_zone_parameters_t {
    /// @brief Координата начала зоны. 
    size_t coordinate_begin;
    size_t coordinate_end;

    /// @brief Внутренний радиус трубопровода
    double r1{ std::numeric_limits<double>::quiet_NaN() };
    /// @brief Внешний радиус эквивалентного слоя изоляции
    double r1eq{ std::numeric_limits<double>::quiet_NaN() };
    /// @brief Внешний радиус эквивалентного слоя грунта
    double r2eq{ std::numeric_limits<double>::quiet_NaN() };
    /// @brief Теплопроводность слоя изоляции
    double conductivity_isolation{ std::numeric_limits<double>::quiet_NaN() };
    /// @brief Теплофизические параметры грунта
    thermophysical_properties_t soil;
    /// @brief Температура термически невозмущенного грунта
    double _temperature_ambient{ KELVIN_OFFSET - 5 };

    /// @brief Коэффициенты для расчетных формул эквивалентных моделей
    equivalent_heat_zone_coefficients_t get_coefficients(
        const heat_zone_adaptation_t& ident) const
    {
        equivalent_heat_zone_coefficients_t result;

        result.temperature_ambient = _temperature_ambient + 100 * (ident.temperature_ambient_soil - 1);

        double soilHeatCapacityVolume = soil.density * soil.heat_capacity * ident.volume_heat_capacity_soil;
        double soilConductivity = soil.conductivity * ident.conductivity_soil;
        double isolationCoductivity = conductivity_isolation * ident.conductivity_isolation;

        // Коэффициент температуропроводности грунта, отнесенный к квадрату внутреннего радиуса трубопровода
        double aSqr = soilConductivity / (soilHeatCapacityVolume * pow(r1, 2));
        double a = sqrt(aSqr);

        // Отношение радиусов
        double al1 = result.al1
            = r1eq / r1; // alpha1
        double al2 = result.al2
            = r2eq / r1eq; // alpha2
        double b = result.b
            = isolationCoductivity / log(al1)
            + soilConductivity / log(al2);


        // Коэффициенты уравнения связи температур T1 и Tsr
        double& A1 = result.A1 = (pow(al2, 2) - pow(al2, 2) * log(al2) - log(al2) - 1) / (8 * aSqr * log(al2));
        double& A2 = result.A2 = (pow(al2, 2) - 2 * log(al2) - 1) / (2 * log(al2) * (pow(al2, 2) - 1));
        double& A3 = result.A3 = isolationCoductivity / (b * log(al1));
        double& A4 = result.A4 = soilConductivity * (2 * log(al2) - pow(al2, 2) + 1);

        // Коэффициенты для выражения теплового потока

        result.A = -2 * M_PI * soilConductivity * isolationCoductivity / (b * log(al1) * log(al2));
        result.B = M_PI * soilConductivity * (2 * log(al2) - pow(al2, 2) + 1) *
            (b * log(al2) - soilConductivity) / (2 * aSqr * log(al2));
        result.C = -A3 * A2 / (A1 + A4 * A2);
        result.F = 1 / (A1 + A4 * A2);

        return result;
    }
};


/// @brief Первичные теплофизические параметры (для зоны трубопровода с одинаковыми свойствами)
struct primary_heat_zone_parameters_t {
    /// @brief Координата начала зоны. 
    /// Конец зоны определеяется либо следующей зоной, либо концом трубопровода
    size_t coordinate_begin;
    /// @brief Толщина стенки, слоев изоляции и защитного слоя, м
    /// стальная стенка трубы, первый слой изоляции, второй слой изоляции, защитный слой
    std::array<double, 4> isolation_thickness{ 0.010, 0.003, 0.050, 0.005 };
    /// @brief Теплопроводность материалов трубопровода, Вт*м-1*К-1
    /// сталь, первый слой изоляции, второй слой изоляции, защитный слой
    std::array<double, 4> isolation_conductivity{ 40, 0.2, 0.035, 0.3 };
    /// @brief Глубина залегания - толщина слоя грунта, м
    double depth = 1.5;
    /// @brief Теплофизические параметры грунта
    thermophysical_properties_t soil;
    /// @brief Температура термически невозмущенного грунта
    double temperature_ambient{ KELVIN_OFFSET - 5 };

    /// @brief Общая толщина физического слоя изоляции
    double get_total_isolation_thickness() const {
        return std::accumulate(isolation_thickness.begin(), isolation_thickness.end(), 0.0);
    }

    /// @brief Расчет коэффициента теплопроводности \lambda всей изоляции
    /// @param pipe_inner_diameter 
    /// @return 
    double get_total_isolation_conductivity(double pipe_inner_diameter) const {
        double Kt_sum = get_total_isolation_HTC(pipe_inner_diameter);

        double r_inner = pipe_inner_diameter / 2;
        double r_outer = r_inner + get_total_isolation_thickness();

        double result = Kt_sum * r_inner * log(r_outer / r_inner);
        return result;
    }

    /// @brief Коэффициент теплообмена всей изоляции, приведенный к внутреннему диаметру трубы
    /// Вклад теплопроводности жидкости исключен ввиду его малости
    /// @param pipe_inner_diameter Внутренний диаметр трубопровода
    double get_total_isolation_HTC(double pipe_inner_diameter) const {
        // Вклад теплопроводности жидкости исключен ввиду его малости
        // считался по формуле: invKt1 = 1 / (a1 * 2 * r1)

        // invKt1 - величина, обратная коэффициенту теплообмена всей изоляции (temp1), 
        // приведенная к внутреннему диаметру трубы
        double invKt1 = 0;
        double r_prev = pipe_inner_diameter / 2;
        for (size_t i = 0; i < isolation_thickness.size(); i++)
        {
            double r_next = r_prev + isolation_thickness[i];
            invKt1 += log(r_next / r_prev) / (2 * isolation_conductivity[i]);
            r_prev = r_next;
        }
        double Kt1 = 1 / (pipe_inner_diameter * invKt1);
        return Kt1;
    }

    /// @brief Коэффициент теплообмена грунта, приведенный к в внешнему диаметру слоя изоляции
    /// @param pipe_inner_diameter Внутренний диаметр трубопровода
    double get_soil_HTC(double pipe_inner_diameter) const {

        double r1 = pipe_inner_diameter / 2; // внутренний радиус трубы
        double rn = r1 + get_total_isolation_thickness(); // внешний радиус слоя изоляции
        double dn = 2 * rn; // внешний диаметр слоя изоляции

        double delta_h = 0.1; // толщина дополнительного (фиктивного)
        // слоя изоляции для учета ухода тепла со свободной поверхности грунта в атмосферу, м

        double h0 = 2 * (depth - rn) + delta_h;

        double invKt2 = log(h0 / dn +
            sqrt(pow((h0) / dn, 2) - 1))
            / (2 * soil.conductivity);

        // invKt2 - величина, обратная коэффициенту теплообмена грунта,
        // приведенная в внешнему диаметру слоя изоляции
        double Kt2 = 1 / (dn * invKt2);
        return Kt2;
    }

    /// @brief Расчет коэффициентов теплообмена эквивалентных слоев Kt1, Kt2
    /// Ошибка в расчете диаметров dsn, оставлено для совместимости
    pair<double, double> get_equivalent_heat_transfers_legacy(double pipe_inner_diameter) const
    {
        // Коэффициент внутренней теплоотдачи. Значение прибито, мало влияет
        double a1 = 257;
        // Внутренний радиус труб ()
        double r1 = pipe_inner_diameter / 2;

        // Коэффициент теплопередачи эквивалентного слоя изоляции
        double temp1 = 1 / (a1 * 2 * r1);
        double temp0 = temp1;

        double ds = 2 * r1;
        for (size_t i = 0; i < isolation_thickness.size(); i++)
        {
            double dsn = ds + isolation_thickness[i];
            temp1 += log(dsn / ds) / (2 * isolation_conductivity[i]);
            ds = dsn;
        }
        double Kt1 = 1 / (2 * r1 * temp1);

        double temp2 = log((2 * (depth - ds / 2) + 0.1) / ds +
            sqrt(pow((2 * (depth - ds / 2) + 0.1) / ds, 2) - 1))
            / (2 * soil.conductivity);

        double Kt2 = 1 / (ds * temp2);

        return std::make_pair(Kt1, Kt2);
    }
    /// @brief Вычисляет радиусы эквивалентных слоев изоляции (r1) и грунта (r2)
    /// @param Kt1 
    /// @param Kt2 
    /// Коэффициент теплопроводности эквивалентного слоя изоляции. 
    /// По умолчанию не задается, а принимается равным теплопроводности реального слоя изоляции
    /// Если задан, то имеет смысл параметра идентификации
    /// @return 
    pair<double, double> get_equivalent_radiuses_legacy(double Kt1, double Kt2,
        double pipe_inner_diameter) const
    {
        // Коэффициенты экв. слоя изоляции, берем самый толстый слой
        double isolatonConductivity = isolation_conductivity[2];

        // Эквивалентные радиусы
        double r1 = pipe_inner_diameter / 2;
        double r1eq = r1 * exp(isolatonConductivity / (Kt1 * r1));
        double r2eq = r1eq * exp(soil.conductivity / (Kt2 * r1eq));

        return std::make_pair(r1eq, r2eq);
    }

    /// @brief Расчет параметров эквивалентной модели тепловой динамики в системе трубопровод-грунт
    /// По сравнению с legacy, 
    ///  a) убрана ошибка в расчете Kt1 для слоев изоляции
    ///  b) убран учет теплообмена в жидкости
    ///  с) изменена логика расчета r1eq
    ///  d) как следствие, расчет r2eq тоже отличается
    /// @param fluid Параметры перекачиваемой жидкости
    /// @return Параметров эквивалентной модели 
    equivalent_heat_zone_parameters_t get_heat_eqivalent_model_alt(double pipe_inner_diameter) const
    {
        equivalent_heat_zone_parameters_t result;
        result.r1 = pipe_inner_diameter / 2;

        // эквивалентирование слоя изоляции цилиндрической стенкой за счет расчета теплопроводности стенки
        result.r1eq = result.r1 + get_total_isolation_thickness();
        result.conductivity_isolation = get_total_isolation_conductivity(pipe_inner_diameter);

        // эквивалентирование грунта цилиндрической стенкой за счет подбора радиуса стенки
        double Kt2 = get_soil_HTC(pipe_inner_diameter);
        result.r2eq = result.r1eq * exp(soil.conductivity / (Kt2 * result.r1eq));
        result.soil = soil;
        result._temperature_ambient = temperature_ambient;
        return result;
    }

    /// @brief Расчет эквивалентной модели для совместимости с первой версией модели
    /// @param pipe_inner_diameter 
    /// @return 
    equivalent_heat_zone_parameters_t get_heat_eqivalent_model_lurie(double pipe_inner_diameter) const
    {
        equivalent_heat_zone_parameters_t result;
        result.r1 = pipe_inner_diameter / 2;

        result.conductivity_isolation = isolation_conductivity[2];

        auto [Kt1, Kt2] = get_equivalent_heat_transfers_legacy(pipe_inner_diameter);
        std::tie(result.r1eq, result.r2eq)
            = get_equivalent_radiuses_legacy(Kt1, Kt2, pipe_inner_diameter);

        result.soil = soil;
        result._temperature_ambient = temperature_ambient;

        return result;
    }

    equivalent_heat_zone_parameters_t get_heat_eqivalent_model(
        HeatModelVer model_version, double pipe_inner_diameter) const
    {
        if (model_version == HeatModelVer::Legacy) {
            return get_heat_eqivalent_model_lurie(pipe_inner_diameter);
        }
        else {
            return get_heat_eqivalent_model_alt(pipe_inner_diameter);
        }
    }

};


/// @brief Трубопровод с зонированными теплофизическими свойствами грунта
struct zoned_pipe_properties : public pipe_properties_t {
    /// @brief Первичные параметры зон (задаются пользователем априорно)
    vector<primary_heat_zone_parameters_t> heat_zones;
    /// @brief Параметры идентификации (задаются пользователем из редактора или в процессе идентификации)
    vector<heat_zone_adaptation_t> adaptation_zones;
    /// @brief Параметры эквивалентных зон 
    /// Вычисляются по первичным параметрам, потом больше не меняются
    /// По сути, функция от первичных зон
    vector<equivalent_heat_zone_parameters_t> eqheat_zones;
    /// @brief Коэффициенты модели для эквивалнтных зон 
    /// (вычисляются по параметрам эквивалентных зон и параметрам идентификации)
    vector<equivalent_heat_zone_coefficients_t> coeff_zones;
    vector<size_t> pressure_sensor_indices;
    HeatModelVer model_version{ HeatModelVer::V2 };

    void init_pressure_sensors(vector<double> coordinates) {
        std::sort(coordinates.begin(), coordinates.end());

        double dx = profile.coordinates[1] - profile.coordinates[0];
        pressure_sensor_indices = { 0 };
        for (double x : coordinates) {
            size_t index = static_cast<int>(x / dx + 0.5);
            pressure_sensor_indices.push_back(index);
        }
        pressure_sensor_indices.push_back(profile.coordinates.size() - 1);
    }

    /// @brief Инциализация эквивалентных моделей, коэффициентов модели с учетом параметров адаптации
    void init_zone_coeffs() {
        eqheat_zones = get_equivalent_zones();
        coeff_zones = get_coeff_zones(adaptation_zones);
    }



    const equivalent_heat_zone_coefficients_t& get_heat_zone_coeff(size_t index) const {
        auto iter = std::find_if(eqheat_zones.rbegin(), eqheat_zones.rend(),
            [index](const equivalent_heat_zone_parameters_t& zone)
            {
                return index >= zone.coordinate_begin;
            });

        if (iter == eqheat_zones.rend())
            throw std::logic_error("wrong zone index");

        size_t zone_index = iter - eqheat_zones.rbegin(); // индекс зоны с конца
        return *(coeff_zones.rbegin() + zone_index);
    }

    const equivalent_heat_zone_parameters_t& get_heat_zone(size_t index) const {

        auto iter = std::find_if(eqheat_zones.rbegin(), eqheat_zones.rend(),
            [index](const equivalent_heat_zone_parameters_t& zone)
            {
                return index >= zone.coordinate_begin;
            });

        if (iter == eqheat_zones.rend())
            throw std::logic_error("wrong zone index");

        return *iter;
    }

    /// @brief Возвращает параметры эквивалентных теплофизических моделей
    vector<equivalent_heat_zone_coefficients_t> get_coeff_zones(
        const vector<heat_zone_adaptation_t>& adapt_zones) const
    {
        vector<equivalent_heat_zone_coefficients_t> result;
        for (size_t index = 0; index < eqheat_zones.size(); ++index) {
            if (adapt_zones.size() == eqheat_zones.size()) {
                result.emplace_back(eqheat_zones[index].get_coefficients(adapt_zones[index]));
            }
            else {
                heat_zone_adaptation_t zone_ident;
                result.emplace_back(eqheat_zones[index].get_coefficients(zone_ident));
            }
        }
        return result;
    }

    /// @brief Параметры по умолчанию - одни единицы
    vector<heat_zone_adaptation_t> get_initial_ident_parameters() const {
        vector<heat_zone_adaptation_t> result(eqheat_zones.size());
        return result;
    }


    /// @brief Возвращает параметры эквивалентных теплофизических моделей
    vector<equivalent_heat_zone_parameters_t> get_equivalent_zones() const {
        vector<equivalent_heat_zone_parameters_t> result;

        //for (const auto& zone : heat_zones) 
        for (size_t index = 0; index < heat_zones.size(); ++index)
        {
            const auto& zone = heat_zones[index];
            auto eqzone = zone.get_heat_eqivalent_model(model_version, wall.diameter);

            eqzone.coordinate_begin = zone.coordinate_begin;
            if (index + 1 < heat_zones.size()) {
                eqzone.coordinate_end = heat_zones[index + 1].coordinate_begin;
            }
            else {
                eqzone.coordinate_end = profile.get_point_count() - 1;
            }

            result.emplace_back(eqzone);
        }

        return result;
    }

};


}