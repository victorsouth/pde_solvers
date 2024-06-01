﻿#pragma once

namespace pde_solvers {
;

using std::string;

/// @brief Проблемно-ориентированный слой для расчета методом конечных объемов 
struct density_viscosity_cell_layer {
    /// @brief Профиль плотности
    std::vector<double> density;
    /// @brief Профиль вязкости
    std::vector<double> viscosity;
    /// @brief Профиль давления
    std::vector<double> pressure;
    /// @brief Дифференциальный профиль давления
    std::vector<double> pressure_delta;
    /// @brief Профиль вспомогательных расчетов для метода конечных объемов (и для вязкости, и для плотности)
    quickest_ultimate_fv_solver_traits<1>::specific_layer specific;
    /// @brief Инициализация профилей
    /// @param point_count Количество точек
    density_viscosity_cell_layer(size_t point_count)
        : density(point_count - 1)
        , viscosity(point_count - 1)
        , specific(point_count)
        , pressure(point_count)
        , pressure_delta(point_count)
    {}

    // @brief Подготовка плотности для расчета методом конечных объемов 
    static quickest_ultimate_fv_wrapper<1> get_density_wrapper(density_viscosity_cell_layer& layer)
    {
        return quickest_ultimate_fv_wrapper<1>(layer.density, layer.specific);
    }
    /// @brief Подготовка вязкости для расчета методом конечных объемов 
    static quickest_ultimate_fv_wrapper<1> get_viscosity_wrapper(density_viscosity_cell_layer& layer)
    {
        return quickest_ultimate_fv_wrapper<1>(layer.viscosity, layer.specific);
    }
};
/// @brief Проблемно-ориентированный слой для расчета методом характеристик
struct density_viscosity_layer_moc {
    /// @brief Профиль плотности
    std::vector<double> density;
    /// @brief Профиль вязкости
    std::vector<double> viscosity;
    /// @brief Профиль давления
    std::vector<double> pressure;
    /// @brief Дифференциальный профиль давления
    std::vector<double> pressure_delta;
    /// @brief Профиль вспомогательных расчетов для МХ (и для вязкости, и для плотности)
    moc_solver<1>::specific_layer specific;
    /// @brief Конструктор на заданное количество точек
    density_viscosity_layer_moc(size_t point_count)
        : density(point_count)
        , viscosity(point_count)
        , specific(point_count)
        , pressure(point_count)
        , pressure_delta(point_count)
    {

    }
    /// @brief Подготовка плотности для расчета методом характеристик
    /// Оборачивает профиль плотности и вспомогательный расчет МХ в обертку для МХ
    static moc_layer_wrapper<1> get_density_wrapper(density_viscosity_layer_moc& layer)
    {
        return moc_layer_wrapper<1>(layer.density, layer.specific);
    }
    /// @brief Подготовка вязкости для расчета методом характеристик
    static moc_layer_wrapper<1> get_viscosity_wrapper(density_viscosity_layer_moc& layer)
    {
        return moc_layer_wrapper<1>(layer.viscosity, layer.specific);
    }
};

/// @brief Структура, созданная для хранения в себе начального профиля давлений и буфера с расчетными данными
template <typename Layer>
struct isothermal_quasistatic_task_buffer_t {
    /// @brief Изначальный профиль давления
    vector<double> pressure_initial;
    /// @brief Буфер профилей давления, плотности, вязкости
    ring_buffer_t<Layer> buffer;
    isothermal_quasistatic_task_buffer_t(size_t point_count)
        : pressure_initial(point_count)
        , buffer(2, point_count)
    {}
};

/// @brief Структура, содержащая в себе начальные условия задачи PQ
struct isothermal_quasistatic_task_boundaries_t {
    /// @brief Изначальный объемный расход
    double volumetric_flow;
    /// @brief Изначальное давление на входе
    double pressure_in;
    /// @brief Изначальная плотность на входе
    double density;
    /// @brief Изначальная вязкость на входе
    double viscosity;

    isothermal_quasistatic_task_boundaries_t() = default;

    isothermal_quasistatic_task_boundaries_t(const vector<double>& values) {
        volumetric_flow = values[0];
        pressure_in = values[1];
        density = values[2];
        viscosity = values[3];
    }

    static isothermal_quasistatic_task_boundaries_t default_values() {
        isothermal_quasistatic_task_boundaries_t result;
        result.volumetric_flow = 0.2;
        result.pressure_in = 6e6;
        result.density = 850;
        result.viscosity = 15e-6;
        return result;
    }
};
/// @brief Расчетная задача (task) для гидравлического изотермического 
/// квазистационарного расчета в условиях движения партий с разной плотностью и вязкостью
/// Расчет партий делается методом характеристик или Quickest-Ultimate
/// @tparam Layer Тип слоя, содержащего профили плотности, вязкости, давления
/// Для партий методом характеристик = density_viscosity_layer_moc, для партий методом Quickest-Ultimate =density_viscosity_cell_layer
/// @tparam Solver Тип солвера партий (moc_solver или quickest_ultimate_fv_solver)
template <typename Layer, typename Solver>
class isothermal_quasistatic_task_t {
    pipe_properties_t pipe;
    isothermal_quasistatic_task_buffer_t<Layer> buffer;

public:
    isothermal_quasistatic_task_t(const pipe_properties_t& pipe)
        : pipe(pipe)
        , buffer(pipe.profile.getPointCount())
    {

    }

    void solve(const isothermal_quasistatic_task_boundaries_t& initial_conditions)
    {
        size_t n = pipe.profile.getPointCount();

        // Инициализация реологии
        auto& current = buffer.buffer.current();

        // Инициализация начального профиля плотности (не важно, ячейки или точки)
        for (double& density : current.density) {
            density = initial_conditions.density;
        }
        // Инициализация начального профиля вязкости (не важно, ячейки или точки)
        for (double& viscosity : current.viscosity) {
            viscosity = initial_conditions.viscosity;
        }

        // Начальный гидравлический расчет
        int euler_direction = +1;
        isothermal_pipe_PQ_parties_t pipeModel(pipe, current.density, current.viscosity, initial_conditions.volumetric_flow, euler_direction);
        solve_euler<1>(pipeModel, euler_direction, initial_conditions.pressure_in, &current.pressure);

        buffer.pressure_initial = current.pressure; // Получаем изначальный профиль давлений
    }
public:
    double get_time_step_assuming_max_speed(double v_max) const {
        const auto& x = pipe.profile.coordinates;
        double dx = x[1] - x[0]; // Шаг сетки
        double dt = abs(dx / v_max); // Постоянный шаг по времени для Куранта = 1
        return dt;
    }
private:
    void make_rheology_step(double dt, const isothermal_quasistatic_task_boundaries_t& boundaries) {
        size_t n = pipe.profile.getPointCount();
        vector<double>Q_profile(n, boundaries.volumetric_flow); // задаем по трубе новый расход из временного ряда

        PipeQAdvection advection_model(pipe, Q_profile);

        auto density_wrapper = buffer.buffer.get_buffer_wrapper(
            &Layer::get_density_wrapper);

        auto viscosity_wrapper = buffer.buffer.get_buffer_wrapper(
            &Layer::get_viscosity_wrapper);

        if constexpr (std::is_same<Solver, moc_solver<1>>::value) {
            // Шаг по плотности
            moc_solver<1> solver_rho(advection_model, density_wrapper);
            solver_rho.step_optional_boundaries(dt, boundaries.density, boundaries.density);
            // Шаг по вязкости
            moc_solver<1> solver_nu(advection_model, viscosity_wrapper);
            solver_nu.step_optional_boundaries(dt, boundaries.viscosity, boundaries.viscosity);

        }
        else {
            // Шаг по плотности
            quickest_ultimate_fv_solver solver_rho(advection_model, density_wrapper);
            solver_rho.step(dt, boundaries.density, boundaries.density);
            // Шаг по вязкости
            quickest_ultimate_fv_solver solver_nu(advection_model, viscosity_wrapper);
            solver_nu.step(dt, boundaries.viscosity, boundaries.viscosity);

        }
    }
    void calc_pressure_layer(const isothermal_quasistatic_task_boundaries_t& boundaries) {
        // Получаем новый профиль давлений

        auto& current = buffer.buffer.current();

        vector<double>& p_profile = current.pressure;
        int euler_direction = +1; // Задаем направление для Эйлера

        isothermal_pipe_PQ_parties_t pipeModel(pipe, current.density, current.viscosity, boundaries.volumetric_flow, euler_direction);
        solve_euler<1>(pipeModel, euler_direction, boundaries.pressure_in, &p_profile);
        // Получаем дифференциальный профиль давлений
        std::transform(buffer.pressure_initial.begin(), buffer.pressure_initial.end(), p_profile.begin(),
            current.pressure_delta.begin(),
            [](double initial, double current) {return initial - current;  });

    }
public:
    void step(double dt, const isothermal_quasistatic_task_boundaries_t& boundaries) {
        make_rheology_step(dt, boundaries);
        calc_pressure_layer(boundaries);
    }
    void advance()
    {
        buffer.buffer.advance(+1);
    }

    auto& get_buffer()
    {
        return buffer.buffer;
    }
protected:
    /// @brief Формирует имя файл для результатов исследования разных численных метов
    /// @tparam Solver Класс солвера
    /// @param path Путь, в котором формируется файл
    /// @return Путь к файлу в заивисимости от указанного класса солвера
    static std::string get_courant_research_filename_for_qsm(const string& path, const string& layer_name)
    {
        std::stringstream filename;
        filename << path << "output " << layer_name << ".csv";
        return filename.str();
    }

    void print(const std::vector<double>& layer, const std::time_t& dt, const std::string& path, const std::string& layer_name)
    {
        std::string filename = get_courant_research_filename_for_qsm(path, layer_name);

        std::ofstream  file(filename, std::ios::app);
        if (file.is_open()) {
            std::tm tm_buf;
            localtime_s(&tm_buf, &dt);
            file << std::put_time(&tm_buf, "%c") << ";";
            for (int j = 0; j < layer.size(); j++)
            {
                file << std::to_string(layer[j]) << ";";
            }
            file << "\n";
            file.close();
        }
    }
public:
    void print_all(const time_t& dt, const string& path) {
        auto& current = buffer.buffer.current();
        print(current.density, dt, path, "density");
        print(current.viscosity, dt, path, "viscosity");
        print(current.pressure, dt, path, "pressure");
        print(current.pressure_delta, dt, path, "pressure_delta");
    }
};


}