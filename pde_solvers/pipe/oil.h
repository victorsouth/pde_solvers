#pragma once

namespace pde_solvers {
;



/// @brief Универсальная реконструкция вязкости по таблице
struct viscosity_table_model_t {
    static constexpr std::array<double, 3> viscosity_temperatures{
        KELVIN_OFFSET + 0, KELVIN_OFFSET + 20, KELVIN_OFFSET + 50 };

    /// @brief Корректирует базовую таблицу вязкости под рабочую вязкость
    /// @param viscosities Базовая таблица вязкости
    /// @param viscosity_working Рабочая вязкость
    /// @param viscosity_working_temperature Температура при которой измеряется вязкость ()
    /// @return Откорректированная таблица вязкости
    static inline array<double, 3> adapt(array<double, 3> viscosities,
        double viscosity_working, double viscosity_working_temperature)
    {
        auto model = reconstruct(viscosities);

        double visc_calc = calc(model, viscosity_working_temperature);
        double mult = viscosity_working / visc_calc;

        for (auto& v : viscosities) {
            v *= mult;
        }

#ifdef _DEBUG
        auto model_adapt = reconstruct(viscosities);
        double visc_adapt = calc(model_adapt, viscosity_working_temperature);
#endif

        return viscosities;
    }


    /// @brief Восстанавливает аппроксимацию по таблице вязкости
    /// @param viscosities 
    /// @return 
    static inline array<double, 3> reconstruct(const array<double, 3>& viscosities)
    {
        array<double, 3> coeffs{
            std::numeric_limits<double>::quiet_NaN(),
            std::numeric_limits<double>::quiet_NaN(),
            std::numeric_limits<double>::quiet_NaN()
        };

        const auto& v = viscosities;

        constexpr double T1 = viscosity_temperatures[0];
        constexpr double T2 = viscosity_temperatures[1];
        constexpr double T3 = viscosity_temperatures[2];

        constexpr double eps1 = 1e-4; // Проверка равенства вязкостей по их отношению (!!надо исследовать!!)
        constexpr double eps2 = 1e-4; // Проверка 

        constexpr double a1 = (T3 - T2) / (T2 - T1);

        double lognu1div2 = log(v[0] / v[1]);
        double nu2div3 = v[1] / v[2];
        double nu3div1 = v[2] / v[0];

        if (abs(nu3div1 - 1) < eps1 && abs(nu2div3 - 1) < eps1) {
            // Вязкость = const
            coeffs[0] = v[0];
            return coeffs;
        }

        double a = a1 * lognu1div2 / log(nu2div3); // похоже, a > 1

        if (abs(a - 1) < eps2) {
            // Филонов-Рейнольс
            coeffs[0] = v[1]; // при двадцати градусах

            constexpr double DT = T2 - T1;
            coeffs[1] = lognu1div2 / DT; // перепроверить!
            return coeffs;
        }

        // Фогель-Фульчер-Тамман
        coeffs[0] = v[2] * pow(nu3div1, 1.0 / (a - 1)); // nu_inf
        coeffs[1] = (T3 - a * T1) / (1 - a); // theta
        const double& theta = coeffs[1];
        coeffs[2] = (T1 - theta) * (T2 - theta) / (T2 - T1) * lognu1div2; // b

        return coeffs;
    }
    /// @brief Расчет вязкости с после реконструкции
    /// @param coeffs 
    /// @param temperature 
    /// @return 
    static inline double calc(const array<double, 3>& coeffs, double temperature)
    {
        if (!std::isnan(coeffs[2])) {
            // Фогель-Фульчер-Тамман
            const auto& nu_inf = coeffs[0];
            const auto& theta = coeffs[1];
            const auto& b = coeffs[2];
            return nu_inf * exp(b / (temperature - theta));
        }
        else if (!std::isnan(coeffs[1])) {
            // Филонов-Рейнольс
            const auto& nu_nominal = coeffs[0];
            const auto& k = coeffs[1];
            return nu_nominal * exp(-k * (temperature - viscosity_temperatures[1]));
        }
        else {
            // Константа
            return coeffs[0];
        }
    }
};




/// @brief Точка вискограммы
struct viscosity_data_point {
    double temperature;
    double kinematic_viscosity;
};


/// @brief Модель вязкости (с потенциалом на общую модель вязкости)
struct oil_viscosity_parameters_t
{
    /// @brief Номинальная температура для вязкости
    double nominal_temperature{ KELVIN_OFFSET + 20 }; // 20 градусов Цельсия
    /// @brief Кинематическая вязкость при номинальной температуре перекачки
    double nominal_viscosity{ 10e-6 };
    /// @brief Коэффициент в формуле Филонова-Рейнольдса
    double temperature_coefficient{ 0 };

    /// @brief Формула вискограммы Филонова-Рейнольдса
    static double viscosity_Filonov_Reynolds(double default_viscosity,
        double default_temperature, double kinematic_viscosity_temperature_coefficient,
        double temperature)
    {
        const double& k = kinematic_viscosity_temperature_coefficient;
        double viscosity = default_viscosity * exp(-k * (temperature - default_temperature));
        return viscosity;
    }
    /// @brief Рассчитывает температурный коэффициент и температуру при 20 град по двум точкам
    /// (T_1, \nu_1), (T_2, \nu_2)
    static double find_kinematic_viscosity_temperature_coefficient(
        double viscosity1, double temperature1,
        double viscosity2, double temperature2)
    {
        double k = -log(viscosity1 / viscosity2) / (temperature1 - temperature2);
        return k;
    }
    /// @brief Расчет вязкости по текущей температуре Филоновым-Рейнольдсом
    double operator()(double temperature) const
    {
        return viscosity_Filonov_Reynolds(nominal_viscosity, nominal_temperature,
            temperature_coefficient, temperature);
    }
    /// @brief Изотермическая вязкость - просто возвращаем при номинальном режиме
    double operator()() const
    {
        return nominal_viscosity;
    }

    /// @brief Инициализация модели вязкости по вискограмме из двух точек
    /// @param viscogramm 
    oil_viscosity_parameters_t(const std::array<viscosity_data_point, 2>& viscogramm)
    {
        const auto& v = viscogramm;
        temperature_coefficient = find_kinematic_viscosity_temperature_coefficient(
            v[0].kinematic_viscosity, v[0].temperature,
            v[1].kinematic_viscosity, v[1].temperature);

        nominal_temperature = v[0].temperature;
        nominal_viscosity = v[0].kinematic_viscosity;
    }
    /// @brief Инициализирует модель вязкости по Филонову-Рейнольдсу по стандартной таблице
    /// @param viscosity_standard_table Таблица стандартных значений вязкости 0, 20, 50 гр. Цельсия
    oil_viscosity_parameters_t(const std::array<double, 3>& viscosity_standard_table) 
    {
        std::array<double, 3> coeffs = viscosity_table_model_t::reconstruct(viscosity_standard_table);
        bool check_model = std::isfinite(coeffs[0])
                        && std::isfinite(coeffs[1])
                        && !std::isfinite(coeffs[2]);
        if (!check_model) {
            // переданные данные не соответствуют модели Филонова-Рейнольдса
            throw std::runtime_error("Filonov-Reynolds data is required");
        }

        nominal_viscosity = coeffs[0];
        temperature_coefficient = coeffs[1];      /// заданная вискограмма
        nominal_temperature = KELVIN_OFFSET + 20;

    }

    /// @brief Дефолтная инициализация
    oil_viscosity_parameters_t() = default;
};

struct oil_density_parameters_t {
    /// @brief Плотность при номинальных условиях, кг/м3
    double nominal_density{ 760 };
    /// @brief Модуль упругости жидкости, Па
    /// Значение по умолчанию приведено по [Лурье 2017], стр. 77
    double fluid_elasticity_modulus{ 1.5e9 };
    /// @brief Номинальное давление, при котором фиксировалась плотность, Па
    double nominal_pressure{ ATMOSPHERIC_PRESSURE };

    /// @brief Коэффициент сжимаемости для жидкости (1/Па)
    /// Коэффициент, учитывающий изменение плотности жидкости при отклонении давления от номинального
    /// В документах обозначает как \beta_\rho
    double getCompressionRatio() const {
        return 1 / fluid_elasticity_modulus;
    }
    /// @brief Плотность при рабочем давлении, с учетом коэффициента сжимаемости жидкости
    double getDensity(double pressure) const {
        return nominal_density * (1 + getCompressionRatio() * (pressure - nominal_pressure));
    }

    /// @brief Плотность без учета сжимаемости и температурного расширения
    double operator()() const
    {
        return nominal_density;
    }
};

/// @brief Теплофизические свойства нефти (переименовать во флюид)
struct oil_heat_parameters_t {
    /// @brief Коэффициент внутренней теплоотдачи, Вт*м-2*К-1
    double internalHeatTransferCoefficient = 257;
    /// @brief Теплоемкость, Дж*кг-1*К-1
    double HeatCapacity = 2000;
    /// @brief Температура застывания, градусы C (переделать на Кельвины!)
    double pourPointTemperature{ 12 };
};


/// @brief Сущность нефти
struct oil_parameters_t {
    /// @brief Модель плотности
    oil_density_parameters_t density;
    /// @brief Модель вязкости
    oil_viscosity_parameters_t viscosity;
    /// @brief Тепловая модель
    oil_heat_parameters_t heat;
    /// @brief Теплоемкость по формуле Крэга
    /// Формула из документа "Полный вывод НЕизотермических..."
    double get_heat_capacity_kreg(double temperature) const {
        double Cp = 31.56 / sqrt(density.nominal_density) * (762 + 3.39 * temperature);
        return Cp;
    }

};


/// @brief Динамические (пересчитываемые в процессе расчета) параметры нефти
/// @tparam DataBuffer Задается vector<double> для профилей, double для точечного случая
template <typename BufferDensity, typename BufferViscosity>
struct fluid_properties_dynamic {
    /// @brief Плотность при стандартных (нормальных, номинальных) условиях
    BufferDensity nominal_density;
    /// @brief Таблица вязкости
    BufferViscosity viscosity_approximation;

    fluid_properties_dynamic(BufferDensity nominal_density, BufferViscosity viscosity_approximation)
        : nominal_density(nominal_density)
        , viscosity_approximation(viscosity_approximation)
    {

    }
};



/// @brief Статические (неизменные в процессе расчета) параметры нефти
struct fluid_properties_static {
    // === Вязкость
    /// @brief Температуры, при которых задавалась вязкость


    // === Плотность
    /// @brief Модуль упругости жидкости, Па
    /// Значение по умолчанию приведено по [Лурье 2017], стр. 77
    double fluid_elasticity_modulus{ 1.5e9 };
    /// @brief Номинальное давление, при котором фиксировалась плотность, Па
    double nominal_pressure{ ATMOSPHERIC_PRESSURE };
    /// @brief Номинальная температура, при которой фиксировалась плотность, K
    double nominal_temperature{ 20 + KELVIN_OFFSET };


    // === Теплофизика
    /// @brief Теплоемкость, Дж*кг-1*К-1
    double heat_capacity = 2000;
    /// @brief Температура застывания, градусы K
    double pour_temperature{ 12 + KELVIN_OFFSET };

    /// @brief Коэффициент сжимаемости для жидкости (1/Па)
    /// Коэффициент, учитывающий изменение плотности жидкости при отклонении давления от номинального
    /// В документах обозначает как \beta_\rho
    double get_compression_ratio() const {
        return 1 / fluid_elasticity_modulus;
    }

};

/// @brief Профиль свойств флюида
struct fluid_properties_profile_t :
    fluid_properties_dynamic<const vector<double>&, const vector<std::array<double, 3>>&>,
    fluid_properties_static
{
    // здесь все функции зависят от координаты (индекса на профиле)

    fluid_properties_profile_t(const vector<double>& nominal_density,
        const vector<std::array<double, 3>>& viscosity_approximation)
        : fluid_properties_dynamic<const vector<double>&, const vector<std::array<double, 3>>&>(nominal_density, viscosity_approximation)
    {

    }

    /// @brief Возвращает вязкость
    /// @param grid_index 
    /// @param temperature 
    /// @return 
    double get_viscosity(size_t grid_index, double temperature) const {
        double result =
            viscosity_table_model_t::calc(viscosity_approximation[grid_index], temperature);
        return result;
    }

};

struct fluid_properties_t :
    fluid_properties_dynamic<double, double>,
    fluid_properties_static
{
    // здесь все функции не зависят от координаты
};

}