#pragma once

namespace pde_solvers {
;

/*
  Структуры данных и расчеты для однозонной модели.
   - Параметры идентификации однозонной модели
   - Исходные параметры Системы "труба, изоляция, грунт"
   - Эквивалентная модель Системы: эквивалентные радиусы слоев, коэффициенты моделей
   - Отнаследованная структура трубопровода, к которому добавлена тепловая модель
   - Вычисления по формулам статей Лурье-Чупраковой
*/

/// @brief Параметры идентификации однозонной эквивалентной модели
struct thermal_model_ident_parameters_t {
    /// @brief Теплопроводность слоя изоляции
    double conductivity_isolation{ std::numeric_limits<double>::quiet_NaN() };
    /// @brief Теплопроводность слоя грунта
    double conductivity_soil{ std::numeric_limits<double>::quiet_NaN() };
    /// @brief Удельная объемная теплоемкость грунта
    double volume_heat_capacity_soil{ std::numeric_limits<double>::quiet_NaN() };
};



/// @brief Параметры тепловой модели системы "труба, изоляция, грунт". 
/// Одна зона грунта
struct pipe_heat_model_t {
    /// @brief Толщина стенки, слоев изоляции и защитного слоя, м
    /// первый слой изоляции, второй слой изоляции, защитный слой
    std::array<double, 4> layerThickness{ 0.010, 0.003, 0.050, 0.005 };

    /// @brief Теплопроводность материалов трубопровода, Вт*м-1*К-1
    /// сталь, первый слой изоляции, второй слой изоляции, защитный слой
    std::array<double, 4> thermalConductivity{ 40, 0.2, 0.035, 0.3 };

    /// @brief Температура окружающей среды, К
    double ambientTemperature{ KELVIN_OFFSET - 5 };

    // @brief Глубина залегания - толщина слоя грунта, м
    double depthOfOccurrence = 1.5;
    /// @brief Теплофизические свойства грунта
    thermophysical_properties_t soil;
    /// @brief Коэффициент теплообмена с окружающей средой
    double ambient_heat_transfer{ 0 };

    /// @brief Простая зависимость теплового потока
    ///  q = -Kt(T - T_ground)
    /// Согласно Лурье 2017, стр. 57 называется формулой Ньютона
    /// @param temperature Температура потока в трубе
    /// @return Тепловой поток
    double get_heat_flow_newton(double temperature) const {
        double q = -ambient_heat_transfer * (temperature - ambientTemperature);
        return q;
    }
    /// @brief Возвращает начальные оценки параметров идентификации, рассчитанные авторским способом
    thermal_model_ident_parameters_t get_ident_parameters() const {
        thermal_model_ident_parameters_t result;
        result.conductivity_isolation = thermalConductivity[2];
        result.conductivity_soil = soil.conductivity;
        result.volume_heat_capacity_soil = soil.density * soil.heat_capacity;
        return result;
    }
    /// @brief Изменяет параметры идентификации
    /// @param ident_parameters 
    void set_ident_parameters(const thermal_model_ident_parameters_t& ident_parameters) {
        thermalConductivity[2] = ident_parameters.conductivity_isolation;
        soil.conductivity = ident_parameters.conductivity_soil;
        soil.heat_capacity = ident_parameters.volume_heat_capacity_soil / soil.density;
    }

};

/// @brief Параметры эквивалентной модели
struct equivalent_thermal_model_t {
    /// @brief Коэффициент теплообмена в эквивалентном слое изоляции
    double Kt1;
    /// @brief Коэффициент теплообмена в эквивалентном слое грунта
    double Kt2;
    /// @brief Внешний радиус эквивалентного слоя изоляции
    double r1eq;
    /// @brief Внешний радиус эквивалентного слоя грунта
    double r2eq;

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
    double get_heat_flow(double T_oil, double T_soil, double ambient_temperature) {
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
    double get_heat_flow(double T_oil, double ambient_temperature) {
        double DT = T_oil - ambient_temperature;
        double q = A * DT;
        return q;
    }
};


/// @brief Трубопровод с тепловой однозонной моделью
struct pipe_noniso_properties_t : public pipe_properties_t {
    /// @brief Тепловая модель
    pipe_heat_model_t heat;



    /// @brief Расчет коэффициентов теплообмена эквивалентных слоев Kt1, Kt2
    /// @param fluid Параметры флюида
    /// @return 
    std::pair<double, double> get_equivalent_heat_transfers(const oil_parameters_t& fluid) const
    {
        // Коэффициент внутренней теплоотдачи
        double a1 = fluid.heat.internalHeatTransferCoefficient;

        // Внутренний радиус труб ()
        //double r1 = (PipeP.wall.diameter - PipeP.heat.layerThickness[0]) / 2;
        double r1 = wall.diameter / 2;

        // Коэффициент теплопередачи эквивалентного слоя изоляции
        double temp1 = 1 / (a1 * 2 * r1);
        double temp0 = temp1;

        double ds = 2 * r1;
        for (size_t i = 0; i < heat.layerThickness.size(); i++)
        {
            double dsn = ds + heat.layerThickness[i];
            temp1 += log(dsn / ds) / (2 * heat.thermalConductivity[i]);
            ds = dsn;
        }
        double Kt1 = 1 / (2 * r1 * temp1);

        double temp2 = log((2 * (heat.depthOfOccurrence - ds / 2) + 0.1) / ds +
            sqrt(pow((2 * (heat.depthOfOccurrence - ds / 2) + 0.1) / ds, 2) - 1))
            / (2 * heat.soil.conductivity);

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
    std::pair<double, double> get_equivalent_radiuses(double Kt1, double Kt2) const
    {
        // Коэффициенты экв. слоя изоляции, берем самый толстый слой
        double isolatonConductivity = heat.thermalConductivity[2];

        // Эквивалентные радиусы
        double r1 = wall.diameter / 2;
        double r1eq, r2eq;
        r1eq = r1 * exp(isolatonConductivity / (Kt1 * r1));
        r2eq = r1eq * exp(heat.soil.conductivity / (Kt2 * r1eq));

        return std::make_pair(r1eq, r2eq);
    }

    /// @brief Расчет параметров эквивалентной модели тепловой динамики в системе трубопровод-грунт
    /// @param fluid Параметры перекачиваемой жидкости
    /// @return Параметров эквивалентной модели 
    equivalent_thermal_model_t get_heat_eqivalent_model(const oil_parameters_t& fluid,
        const thermal_model_ident_parameters_t& ident = thermal_model_ident_parameters_t()
    ) const
    {
        equivalent_thermal_model_t result;
        double r1 = wall.diameter / 2;

        std::tie(result.Kt1, result.Kt2) = get_equivalent_heat_transfers(fluid);
        std::tie(result.r1eq, result.r2eq) = get_equivalent_radiuses(result.Kt1, result.Kt2);

        double isolatonConductivity = heat.thermalConductivity[2];
        double soilConductivity = heat.soil.conductivity;
        double soilHeatCapacityVolume = heat.soil.density * heat.soil.heat_capacity;
        if (isfinite(ident.conductivity_isolation))
            isolatonConductivity = ident.conductivity_isolation;
        if (isfinite(ident.conductivity_soil))
            soilConductivity = ident.conductivity_soil;
        if (isfinite(ident.volume_heat_capacity_soil))
            soilHeatCapacityVolume = ident.volume_heat_capacity_soil;


        // Коэффициент температуропроводности грунта, отнесенный к квадрату внутреннего радиуса трубопровода
        double aSqr = soilConductivity / (soilHeatCapacityVolume * pow(r1, 2));
        double a = sqrt(aSqr);

        // Отношение радиусов
        double al1 = result.al1
            = result.r1eq / r1; // alpha1
        double al2 = result.al2
            = result.r2eq / result.r1eq; // alpha2
        double b = result.b
            = isolatonConductivity / log(al1)
            + soilConductivity / log(al2);


        // Коэффициенты уравнения связи температур T1 и Tsr
        double& A1 = result.A1;
        double& A2 = result.A2;
        double& A3 = result.A3;
        double& A4 = result.A4;
        A1 = (pow(al2, 2) - pow(al2, 2) * log(al2) - log(al2) - 1) / (8 * aSqr * log(al2));
        A2 = (pow(al2, 2) - 2 * log(al2) - 1) / (2 * log(al2) * (pow(al2, 2) - 1));
        A3 = isolatonConductivity / (b * log(al1));
        A4 = soilConductivity * (2 * log(al2) - pow(al2, 2) + 1);

        // Коэффициенты для выражения теплового потока

        result.A = -2 * M_PI * soilConductivity * isolatonConductivity / (b * log(al1) * log(al2));
        result.B = M_PI * soilConductivity * (2 * log(al2) - pow(al2, 2) + 1) *
            (b * log(al2) - soilConductivity) / (2 * aSqr * log(al2));
        result.C = -A3 * A2 / (A1 + A4 * A2);
        result.F = 1 / (A1 + A4 * A2);

        return result;
    }

};



/// @brief Расчет профиля средней темературы грунта по Лурье, Чупракова 2021, формула 14
inline void compute_soil_mean_temperature_distribution(
    const pipe_noniso_properties_t& pipe,
    const oil_parameters_t& oil,
    const std::vector<double>& temperature_oil,
    std::vector<double>* _result,
    const thermal_model_ident_parameters_t& ident = thermal_model_ident_parameters_t()
)
{
    std::vector<double>& temperature = *_result;

    double lambda = pipe.heat.thermalConductivity[2]; //Стоит разобраться, почему выбирают именно этот слой

    equivalent_thermal_model_t heat_dynamic_model = pipe.get_heat_eqivalent_model(oil, ident);

    double al1 = heat_dynamic_model.al1;
    double al2 = heat_dynamic_model.al2;
    double b = heat_dynamic_model.b;

    //Расчет коэффициента связи средней температуры грунта и температуры флюида
    double K = lambda * (pow(al2, 2) - 2 * log(al2) - 1) / (2 * b * log(al1) * log(al2) * (pow(al2, 2) - 1));
    double K1 = -heat_dynamic_model.C / heat_dynamic_model.F;

    // Расчет распределения Тср
    for (size_t i = 0; i < temperature.size(); i++)
    {
        // среднее превышение температурой гранута Tгр(r) наружной температуры Tнар = Tгр(r_2)
        double DT1 = temperature_oil[i] - pipe.heat.ambientTemperature;
        temperature[i] = K1 * DT1 + pipe.heat.ambientTemperature;
    }
}

/// @brief Пересчет профиля темературы флюида по Лурье, Чупракова 2021, стр. 6
inline void compute_new_fluid_temperature_distribution(
    //const PipeHeatInflowConstArea& heat_pde,
    const pipe_noniso_properties_t& pipe,
    const oil_parameters_t& oil,
    double mass_flow,
    double Tin,
    double dt,
    const std::vector<double>& _heat_flow,
    const std::vector<double>& temperature_oil_old,
    std::vector<double>* temperature_oil_new)
{
    const std::vector<double>& qT = _heat_flow;
    std::vector<double>& T1 = *temperature_oil_new;

    double d = pipe.wall.diameter;
    double S_0 = pipe.wall.getArea();
    double density = oil.density.nominal_density;
    double C = oil.heat.HeatCapacity;


    T1[0] = Tin;
    for (size_t i = 0; i < T1.size() - 1; i++)
    {
        //double s1 = heat_pde.getSourceTerm(i, temperature_oil_old[i]);

        double v = mass_flow / (S_0 * density);
        double Re = v * pipe.wall.diameter / oil.viscosity(temperature_oil_old[i]);
        double lambda = hydraulic_resistance_shifrinson(Re, pipe.wall.relativeRoughness());

        double q = 0.5 * (qT[i] + qT[i + 1]);
        double s = ((4 * q / (density * C * d)) + lambda * pow(v, 3) / (2 * C * d));
        double dT = dt * s;
        T1[i + 1] = temperature_oil_old[i] + dT;
    }
}

/// @brief Пересчет профиля темературы флюида. 
/// Попытка избавиться от артефакта метода
inline void compute_new_fluid_temperature_distribution2(
    //const PipeHeatInflowConstArea& heat_pde,
    const pipe_noniso_properties_t& pipe,
    const oil_parameters_t& oil,
    double mass_flow,
    double Tin,
    double dt,
    const std::vector<double>& _heat_flow,
    const std::vector<double>& temperature_oil_old,
    const std::vector<double>& Tsr,
    std::vector<double>* temperature_oil_new,
    const thermal_model_ident_parameters_t& ident
)
{
    const std::vector<double>& qT = _heat_flow;
    std::vector<double>& T1 = *temperature_oil_new;

    double d = pipe.wall.diameter;
    double S_0 = pipe.wall.getArea();
    double density = oil.density.nominal_density;
    double C_oil = oil.heat.HeatCapacity;

    equivalent_thermal_model_t heat_dynamic_model = pipe.get_heat_eqivalent_model(oil, ident);
    double A = heat_dynamic_model.A;
    double B = heat_dynamic_model.B;
    double C = heat_dynamic_model.C;
    double F = heat_dynamic_model.F;

    T1[0] = Tin;
    for (size_t i = 0; i < T1.size() - 1; i++)
    {
        auto right_party = [&](size_t i, double T) {
            double v = mass_flow / (S_0 * density);
            double Re = v * pipe.wall.diameter / oil.viscosity(T);
            double lambda = hydraulic_resistance_shifrinson(Re, pipe.wall.relativeRoughness());

            double DT1 = T - pipe.heat.ambientTemperature;
            double DTsr = Tsr[i] - pipe.heat.ambientTemperature;

            double q = (A + B * C) * DT1 + F * B * DTsr;

            double s = ((4 * q / (density * C_oil * d)) + lambda * pow(v, 3) / (2 * C_oil * d));
            return s;
        };

        //double s1 = heat_pde.getSourceTerm(i, temperature_oil_old[i]);
        //double q = 0.5 * (qT[i] + qT[i + 1]);
        //double s = ((4 * q / (density * C * d)) + lambda * pow(v, 3) / (2 * C * d));

        double rp1 = right_party(i, temperature_oil_old[i]);
        double Tnext_estimate = temperature_oil_old[i] + dt * rp1;
        double rp2 = right_party(i + 1, Tnext_estimate);

        T1[i + 1] = temperature_oil_old[i] + dt * 0.5 * (rp1 + rp2);
    }
}


/// @brief Пересчет профиля средней темературы грунта по Лурье, Чупракова 2021, стр. 6
inline void compute_new_soil_mean_temperature_distribution(
    const pipe_noniso_properties_t& pipe,
    const oil_parameters_t& oil,
    double dt,
    const std::vector<double>& _temperature,
    const std::vector<double>& _mean_dirt_temperature,
    std::vector<double>* _result,
    const thermal_model_ident_parameters_t& ident = thermal_model_ident_parameters_t()

)
{
    const std::vector<double>& T1 = _temperature;
    const std::vector<double>& Tsr_old = _mean_dirt_temperature;
    std::vector<double>& Tsr = *_result;

    equivalent_thermal_model_t heat_dynamic_model = pipe.get_heat_eqivalent_model(oil, ident);
    double C = heat_dynamic_model.C;
    double F = heat_dynamic_model.F;

    for (size_t i = 0; i < Tsr_old.size(); i++)
    {
        double DT1 = T1[i] - pipe.heat.ambientTemperature;
        double DTsr = Tsr_old[i] - pipe.heat.ambientTemperature;
        double dT = dt * (C * DT1 + F * DTsr);
        Tsr[i] = Tsr_old[i] + dT;
    }
}

/// @brief Расчет теплового потока по Лурье, Чупракова 2021, формула 13
inline void compute_heat_flow_distribution(
    const pipe_noniso_properties_t& pipe,
    const oil_parameters_t& oil,
    const std::vector<double>& temperature,
    const std::vector<double>& mean_amb_temperature,
    std::vector<double>* _qT,
    const thermal_model_ident_parameters_t& ident = thermal_model_ident_parameters_t()
)
{
    const std::vector<double>& T1 = temperature;
    const std::vector<double>& Tsr = mean_amb_temperature;
    std::vector<double>& qT = *_qT;

    equivalent_thermal_model_t heat_dynamic_model = pipe.get_heat_eqivalent_model(oil, ident);
    double A = heat_dynamic_model.A;
    double B = heat_dynamic_model.B;
    double C = heat_dynamic_model.C;
    double F = heat_dynamic_model.F;

    // Расчет распределения qT
    for (size_t i = 0; i < qT.size(); i++)
    {
        //double qt_furie = pipe.heat.get_heat_flow_newton(T1[i]); // для сравнения

        double DT1 = T1[i] - pipe.heat.ambientTemperature;
        double DTsr = Tsr[i] - pipe.heat.ambientTemperature;

        double qt_lurie_st = A * DT1;

        double qt_lurie = (A + B * C) * DT1 + F * B * DTsr;
        qT[i] = qt_lurie;
    }
}

/// @brief Расчет граничной температуры между изоляцией и грунтом
inline void compute_isolation_boundary_temperature(
    const pipe_noniso_properties_t& pipe,
    const oil_parameters_t& oil,
    const std::vector<double>& temperature,
    const std::vector<double>& mean_amb_temperature,
    std::vector<double>* _Tiso,
    const thermal_model_ident_parameters_t& ident = thermal_model_ident_parameters_t()
)
{
    const std::vector<double>& T1 = temperature;
    const std::vector<double>& Tsr = mean_amb_temperature;
    std::vector<double>& Tiso = *_Tiso;

    equivalent_thermal_model_t heat_dynamic_model = pipe.get_heat_eqivalent_model(oil, ident);
    double A3 = heat_dynamic_model.A3;
    double A4 = heat_dynamic_model.A4;
    double C = heat_dynamic_model.C;
    double F = heat_dynamic_model.F;

    // Расчет распределения qT
    for (size_t i = 0; i < Tiso.size(); i++)
    {
        double DT1 = T1[i] - pipe.heat.ambientTemperature;
        double DTsr = Tsr[i] - pipe.heat.ambientTemperature;

        Tiso[i] = (A3 + A4 * C) * DT1 + A4 * F * DTsr + pipe.heat.ambientTemperature;
    }
}


}