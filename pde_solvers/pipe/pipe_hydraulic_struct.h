#pragma once
#include <vector>
#include <stdexcept>
namespace pde_solvers {

/// @brief Упрощенные параметры трубы
struct simple_pipe_properties {
    double length = 12000;
    double dx = 1000;
    double elevation = 0;
    double diameter = 0.7;
    /// @brief Количество сегментов при дроблении длины length с шагом dx 
    /// с округлением до ближайшего
    size_t get_segment_count() const {
        return static_cast<size_t>(0.5 + length / dx);
    }
    /// @brief Имитация ТУ для расчетов движения партий
    /// длина 700км, диаметр 514мм, горизонтальная, сетка 100 м
    static simple_pipe_properties sample_district() {
        simple_pipe_properties simple_pipe;
        simple_pipe.length = 700e3; 
        simple_pipe.diameter = 0.514;
        simple_pipe.dx = 100; 
        return simple_pipe;
    }
    /// @brief Имитация участка трубы для синтетических проверок расчетов движения партий
    /// 50км, диаметр 700мм, сетка 100м.
    static simple_pipe_properties sample_section()
    {
        simple_pipe_properties simple_pipe;
        simple_pipe.length = 400e3;
        simple_pipe.diameter = 0.7;
        simple_pipe.dx = 100;
        return simple_pipe;
    }
};

/// @brief Сущность профиля трубы
struct pipe_profile_t {
    /// @brief Координатные отметки, м
    std::vector<double> coordinates;
    /// @brief Высотные отметки, м
    std::vector<double> heights;
    /// @brief Несущая способность, Па
    std::vector<double> capacity;


    double get_height_derivative(ptrdiff_t index, int direction) const {
        ptrdiff_t neighbour_index = index + direction;

        if (neighbour_index < 0 || neighbour_index >= static_cast<ptrdiff_t>(get_point_count()))
            throw std::runtime_error("Wrong neighbour profile index");

        if (index < 0 || index >= static_cast<ptrdiff_t>(get_point_count()))
            throw std::runtime_error("Wrong profile index");

        double dx = coordinates[neighbour_index] - coordinates[index];
        double dz = heights[neighbour_index] - heights[index];
        return dz / dx;

    }

    /// @brief Длина участка трубы, м
    double get_length() const
    {
        return coordinates.back() - coordinates.front();
    }
    /// @brief Количество границ расчетной сетки
    size_t get_point_count() const
    {
        return coordinates.size();
    }
    /// @brief Создает профиль с линейным уклоном по параметрам в начале и в конце
    /// @param segment_count 
    /// @param x_begin 
    /// @param x_end 
    /// @param z_begin 
    /// @param z_end 
    /// @param capacity 
    /// @return 
    static pipe_profile_t create(size_t segment_count,
        double x_begin, double x_end, double z_begin, double z_end,
        double capacity)
    {
        size_t n = segment_count + 1;
        pipe_profile_t result;
        result.coordinates = result.heights = std::vector<double>(n);
        result.capacity = std::vector<double>(n, capacity);
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


    /// @brief Формула расчета гидравлического сопротивления
    double(*resistance_function_pointer)(double, double) { hydraulic_resistance_isaev };
    double resistance_function_adaptation{ 1 };

    /// @brief Формула расчета гидравлического сопротивления
    /// Параметры шероховатости подхватывается из этой структуры
    double resistance_function(double Re) const {
        double eps = relativeRoughness();
        double lambda = resistance_function_pointer(Re, eps);
        lambda *= resistance_function_adaptation;
        return lambda;
    }


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
    double getPressureArea(const oil_parameters_t& oil, double nominalArea, double pressure) const
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
struct pipe_properties
{
    /// @brief Профиль для расчета (сетка). Не исходный профиль
    pipe_profile_t profile;
    /// @brief Модель для расчета стенок
    pipe_wall_model_t wall;
    /// @brief Параметры адаптации
    AdaptationParameters adaptation;
    /// @brief Формула расчета гидравлического сопротивления
    /// Параметры шероховатости подхватывается из wall
    double resistance_function(double Re) const {
        return wall.resistance_function(Re);
    }


    /// @brief Скорость звука в жидкости, м/с^2
    /// TODO: указать источник литературы
    double getSoundVelocity(const oil_parameters_t& oil) const
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
    double getNominalArea(const oil_parameters_t& oil, double pressureArea, double density) const
    {
        return (density / oil.density()) * pressureArea;
    }
    /// @brief Чувствительность к давлению? \theta
    /// \param oil
    /// \param pressure
    /// \return
    double getTeta(const oil_parameters_t& oil, double pressure) const
    {
        return 1 + (oil.density.getCompressionRatio() + wall.getCompressionRatio())
            * (pressure - wall.nominal_pressure);
    }

    /// @brief Трубопровод по умолчанию 
    inline static pipe_properties<AdaptationParameters> build_simple_pipe(
        const simple_pipe_properties& simple)
    {
        pipe_properties<AdaptationParameters> pipe;

        double Pcapacity = 10e6; // несущая способность
        size_t segment_count = simple.get_segment_count();
        pipe.profile = pipe_profile_t::create(segment_count, 0, simple.length, 0, simple.elevation, Pcapacity);
        pipe.wall.equivalent_roughness = 0.0001;

        // это диаметр внутренний
        pipe.wall.diameter = simple.diameter;
        // это толщина одной стенки, к внешнему надо прибавлять удвоенную толщину
        pipe.wall.wallThickness = 0.01;
        return pipe;
    }

};

typedef pipe_properties<adaptation_parameters> pipe_properties_t;

}
