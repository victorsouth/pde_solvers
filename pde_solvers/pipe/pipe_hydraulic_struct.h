#pragma once



/// @brief Сущность профиля трубы
struct PipeProfile {
    /// @brief Координатные отметки, м
    vector<double> coordinates;
    /// @brief Высотные отметки, м
    vector<double> heights;
    /// @brief Несущая способность, Па
    vector<double> capacity;
    /// @brief Длина участка трубы, м
    double getLength() const
    {
        return coordinates.back() - coordinates.front();
    }
    /// @brief Количество границ расчетной сетки
    size_t getPointCount() const
    {
        return coordinates.size();
    }

    static PipeProfile create(size_t segment_count,
        double x_begin, double x_end, double z_begin, double z_end,
        double capacity)
    {
        size_t n = segment_count + 1;
        PipeProfile result;
        result.coordinates = result.heights = vector<double>(n);
        result.capacity = vector<double>(n, capacity);
        double length = x_end - x_begin;
        double dx = length / segment_count;
        for (size_t index = 0; index < n; ++index) {
            double x = x_begin + dx * index;
            result.coordinates[index] = x;

            double alpha = 1 - x / length;
            double z = z_begin * alpha + z_end * (1 - alpha);
            result.heights[index] = z;
        }
        return result;
    }
};

/// @brief Параметры стенки трубы
struct pipe_wall_model_t {
    /// @brief Внутренний диаметр, м
    double diameter{ 0.8 };
    /// @brief Толщина стенок, м
    double wallThickness{ 10e-3 }; // 10 мм
    /// @brief Эквивалентная шероховатость, м (не в мм!)
    double equivalent_roughness{ 0.125e-3 }; // перепроверить по РД
    /// @brief Номинальное давление, при котором фиксировался диаметр, Па
    double nominal_pressure{ 101325 }; // что если номинальное давление для нефти отличается?
    /// @brief Модуль упругости Юнга материала труб 
    /// (значение: Лурье 2017, стр. 92)
    double wall_elasticity_modulus{ 2e11 };
    /// @brief коэффициент упругости Пуассона
    /// Значение по умолчанию взято из [Лурье 2017], стр. 93
    double poissonStrain{ 0.28 };

    /// @brief Коэффициент сжимаемости для трубы (1/Па)
    /// Коэффициент, учитывающий изменение площади сечения при отклонении давления от номинального
    /// В документах обозначает как \beta_S
    double getCompressionRatio() const {
        return diameter * (1 - poissonStrain * poissonStrain) / (wall_elasticity_modulus * wallThickness);
    }
    /// @brief Относительная шероховатость
    double relativeRoughness() const {
        return equivalent_roughness / diameter;
    }
    /// @brief Площадь сечения \S_0, м^2
    double getArea() const
    {
        return (M_PI * diameter * diameter) / 4;
    }

    /// @brief Площадь сечения \S, зависящая от давления, м^2
    /// \param oil
    /// \param area
    /// \param pressure
    /// \return
    double getPressureArea(const OilParameters& oil, double nominalArea, double pressure) const
    {
        double beta_rho = getCompressionRatio();
        double multiplier = (1 + beta_rho * (pressure - nominal_pressure));
        return nominalArea * multiplier;
    }
};

/// @brief Параметры адаптации (в тепловых зонах есть еще!!)
struct adaptation_parameters {
    /// @brief Адаптация трения для коррекции формулы Лямбды
    double friction{ 1 };
    /// @brief Поправка на состояние трубопровода
    double diameter{ 1 };
    /// @brief Поправка на коэффициент теплообмена флюида
    double heat_capacity{ 1 };
};


/// @brief Сущность трубы
template <typename AdaptationParameters>
struct pipe_properties_t
{
    /// @brief Профиль для расчета (сетка). Не исходный профиль
    PipeProfile profile;
    /// @brief Модель для расчета стенок
    pipe_wall_model_t wall;

    /// @brief Параметры адаптации
    AdaptationParameters adaptation;


    double(*resistance_function)(double, double) { hydraulic_resistance_isaev };


    /// @brief Скорость звука в жидкости, м^2/с
    double getSoundVelocity(const OilParameters& oil) const
    {
        double beta_S = wall.getCompressionRatio();
        double beta_rho = oil.density.getCompressionRatio();
        double c = sqrt(1 / (oil.density() * (beta_S + beta_rho)));
        return c;
    }

    /// @brief Скорость звука в жидкости
    double get_sound_velocity(double beta_rho, double nominal_density) const
    {
        double beta_S = wall.getCompressionRatio();
        double c = sqrt(1 / (nominal_density * (beta_S + beta_rho)));
        return c;
    }

    /// @brief Умозрительная площадь сечения \A, которая нужна, чтобы поместить жидкость при номинальных условиях, м^2
    /// \param oil
    /// \param pressureArea
    /// \param density
    /// \return
    double getNominalArea(OilParameters oil, double pressureArea, double density) const
    {
        return (density / oil.density()) * pressureArea;
    }
    /// @brief \teta
    /// \param oil
    /// \param pressure
    /// \return
    double getTeta(OilParameters oil, double pressure) const
    {
        return 1 + (oil.density.getCompressionRatio() + wall.getCompressionRatio())
            * (pressure - wall.nominal_pressure);
    }
};

typedef pipe_properties_t<adaptation_parameters> PipeProperties;
