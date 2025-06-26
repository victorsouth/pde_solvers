#pragma once
#include <string>
#include <vector>

namespace pde_solvers {
;

/// @brief Варианты расчёта квазистационарной модели
/// Stationary - стационарный расчёт (частный случай квазистаца)
/// DensityQuasi - квазистац по плотности, ост как в стационаре
/// ViscosityQuasi - квазистац по вязкости, ост как в стационаре
/// TempQuasi - квазистац по температуре, ост как в стац
/// FullQuasi - полный квазистац
enum class noniso_qsm_model_type {
    Stationary, DensityQuasi, ViscosityQuasi, FullQuasi, TempQuasi
};

/// @brief Проблемно-ориентированный слой для TQ расчета
/// @tparam CellFlag Флаг расчёта реологии 
/// true - в ячейках для метода конечных объёмов (Quickest-Ultimate)
/// false - в точках для метода характеристик (advection_moc_solver)
template <bool CellFlag>
struct density_viscosity_temp_quasi_layer {
    /// @brief Профиль плотности
    std::vector<double> density;
    /// @brief Профиль вязкости
    std::vector<double> viscosity;
    /// @brief Профиль температуры
    std::vector<double> temperature;
    /// @brief Профиль температуры  Шуховым
    std::vector<double> temperature_shukhov;
    /// @brief Профиль вспомогательных расчетов для метода конечных объемов (и для вязкости, и для плотности)
    quickest_ultimate_fv_solver_traits<1>::specific_layer specific;
    /// @brief Профиль вспомогательных расчетов для МХ (и для температуры)
    moc_solver<1>::specific_layer moc_specific;
    /// @brief Инициализация профилей
    /// @param point_count Количество точек
    density_viscosity_temp_quasi_layer(size_t point_count)
        : density(point_count - static_cast<int>(CellFlag))
        , viscosity(point_count - static_cast<int>(CellFlag))
        , temperature(point_count - static_cast<int>(CellFlag))
        , temperature_shukhov(point_count)
        , specific(point_count)
        , moc_specific(point_count)
    {
    }

    /// @brief Подготовка плотности для расчета методом конечных объемов 
    /// @param layer Слой
    /// @return Обертка над составным слоем
    static quickest_ultimate_fv_wrapper<1> get_density_wrapper(density_viscosity_temp_quasi_layer& layer)
    {
        return quickest_ultimate_fv_wrapper<1>(layer.density, layer.specific);
    }

    /// @brief Подготовка вязкости для расчета методом конечных объемов
    /// @param layer Слой
    /// @return Обертка над составным слоем
    static quickest_ultimate_fv_wrapper<1> get_viscosity_wrapper(density_viscosity_temp_quasi_layer& layer)
    {
        return quickest_ultimate_fv_wrapper<1>(layer.viscosity, layer.specific);
    }

    /// @brief Подготовка температуры для расчета методом конечных объемов
    /// @param layer Слой
    /// @return Обертка над составным слоем
    static quickest_ultimate_fv_wrapper<1> get_temperature_wrapper(density_viscosity_temp_quasi_layer& layer)
    {
        return quickest_ultimate_fv_wrapper<1>(layer.temperature, layer.specific);
    }


    /// @brief Подготовка температуры для расчета по методу характеристик
    /// Оборачивает профиль температуры и вспомогательный расчет МХ в обертку для МХ
    static moc_layer_wrapper<1> get_temperature_moc_wrapper(density_viscosity_temp_quasi_layer& layer)
    {
        return moc_layer_wrapper<1>(layer.temperature, layer.moc_specific);
    }
};

/// @brief Структура, содержащая в себе краевые условия задачи T?
struct nonisothermal_quasistatic_PQ_task_boundaries_t { 
    /// @brief Изначальный объемный расход
    double volumetric_flow;
    /// @brief Изначальное температура на входе
    double temperature;
    /// @brief Изначальная плотность на входе
    double density;
    /// @brief Изначальная вязкость на входе
    double viscosity;
    /// @brief Изначальная темп на входе
    double temperature_shukhov;

    /// @brief Конструктор по умолчанию
    nonisothermal_quasistatic_PQ_task_boundaries_t() = default;

    /// @brief Конструктор краевых условий
    /// @param values Значения краевых условий
    nonisothermal_quasistatic_PQ_task_boundaries_t(const vector<double>& values) {
        volumetric_flow = values[0];
        temperature = values[1];
        density = values[2];
        viscosity = values[3];
        temperature_shukhov = values[1];
    }

    /// @brief Создание структуры со значениями по умолчанию
    static nonisothermal_quasistatic_PQ_task_boundaries_t default_values() {
        nonisothermal_quasistatic_PQ_task_boundaries_t result;
        result.volumetric_flow = 0.2;
        result.temperature = 300;
        result.density = 850;
        result.viscosity = 15e-6;
        result.temperature_shukhov = 300;
        return result;
    }
};


/// @brief Расчетная задача (task) для гидравлического неизотермического 
/// квазистационарного расчета в условиях движения партий с разной плотностью и вязкостью, температурой
/// Расчет партий делается методом характеристик или Quickest-Ultimate
/// @tparam Solver Тип солвера партий (advection_moc_solver или quickest_ultimate_fv_solver)
template <typename Solver = advection_moc_solver>
class nonisothermal_quasistatic_PQ_task_t {  
public:
    /// @brief Тип слоя
    using layer_type = density_viscosity_temp_quasi_layer<std::is_same<Solver, advection_moc_solver>::value ? false : true>;
    /// @brief Тип буфера
    using buffer_type = ring_buffer_t<layer_type>;
    /// @brief Тип граничных условий
    using boundaries_type = nonisothermal_quasistatic_PQ_task_boundaries_t;
private:
    // Модель трубы
    pipe_noniso_properties_t pipe;
    // Нефть
    oil_parameters_t oil;
    // Создаётся буфер, тип слоя которого определяется в зависимости от типа солвера
    buffer_type buffer;
    // Тип расчёта
    noniso_qsm_model_type model_type;

public:
    /// @brief Конструктор
    /// @param pipe Модель трубопровода
    nonisothermal_quasistatic_PQ_task_t(const pipe_noniso_properties_t& pipe, const oil_parameters_t& oil, noniso_qsm_model_type model_type = noniso_qsm_model_type::FullQuasi)
        : pipe(pipe)
        , oil(oil)
        , buffer(2, pipe.profile.get_point_count()) //2
        , model_type(model_type)
    {
    }

    /// @brief Геттер для текущего слоя  
    density_viscosity_temp_quasi_layer<std::is_same<Solver, advection_moc_solver>::value ? false : true> get_current_layer() {
        return buffer.current();
    }

    /// @brief Тип режима расчета на шаге
    enum class step_mode_t {
        Advection,            /// Метод адвеции
        Shukhov,              /// Метод Шухова
        ShukhovWithAdvection  /// Шухов c адвекцией
    };

    /// @brief Начальный стационарный расчёт. 
    /// Ставим по всей трубе реологию из initial_conditions, делаем гидравлический расчет
    /// @param initial_conditions Начальные условия
    void solve(const nonisothermal_quasistatic_PQ_task_boundaries_t& initial_conditions)
    {
        // Количество точек
        size_t n = pipe.profile.get_point_count();

        // Инициализация реологии
        auto& current = buffer.current();

        // Инициализация начального профиля плотности (не важно, ячейки или точки)
        for (double& density : current.density) {
            density = initial_conditions.density;
        }
        // Инициализация начального профиля вязкости (не важно, ячейки или точки)
        for (double& viscosity : current.viscosity) {
            viscosity = initial_conditions.viscosity;
        }
        // Инициализация начального профиля температуры (не важно, ячейки или точки)
        for (double& temperature : current.temperature) {
            temperature = initial_conditions.temperature;
        }
        // Инициализация начального профиля температуры (не важно, ячейки или точки)
        for (double& temperature_shukhov : current.temperature_shukhov) {
            temperature_shukhov = initial_conditions.temperature_shukhov;
        }
    }
public:
    /// @brief Рассчёт шага по времени для Cr = 1
    /// @param v_max Максимальная скорость течение потока в трубопроводе
    double get_time_step_assuming_max_speed(double v_max) const {
        const auto& x = pipe.profile.coordinates;
        double dx = x[1] - x[0]; // Шаг сетки
        double dt = abs(dx / v_max); // Постоянный шаг по времени для Куранта = 1
        return dt;
    }
private:
    /// @brief Проводится рассчёт шага движения партии
    /// @param dt Временной шаг моделирования
    /// @param boundaries Краевые условия
    void make_rheology_step_advection(double dt, const nonisothermal_quasistatic_PQ_task_boundaries_t& boundaries) {
        
        size_t n = pipe.profile.get_point_count();
        vector<double>Q_profile(n, boundaries.volumetric_flow); /// задаем по трубе новый расход из временного ряда
        vector<double> G(pipe.profile.get_point_count(), Q_profile[0] * oil.density.nominal_density);   /// массовый расход
        auto heatModel = std::make_unique<PipeHeatInflowConstArea>(pipe, oil, G);   

        advance(); // Сдвигаем текущий и предыдущий слои

        if constexpr (std::is_same<Solver, advection_moc_solver>::value) {
            // считаем партии методом характеристик
            // температура - квазистац или стац
            if (model_type == noniso_qsm_model_type::FullQuasi || model_type == noniso_qsm_model_type::TempQuasi) {
                // Шаг по temp
                auto temperature_buffer = buffer.get_buffer_wrapper(&density_viscosity_temp_quasi_layer<0>::get_temperature_moc_wrapper);               

                moc_solver<1> solver_tm(*heatModel, temperature_buffer);
                solver_tm.step_optional_boundaries(dt, boundaries.temperature, boundaries.temperature);
            }
            else
                buffer.current().temperature = vector<double>(buffer.current().temperature.size(), boundaries.temperature);
        }
        else {
            // считаем партии с помощью QUICKEST-ULTIMATE       
            // temp - квазистац или стац
            if (model_type == noniso_qsm_model_type::FullQuasi || model_type == noniso_qsm_model_type::TempQuasi) {
                auto temp_wrapper = buffer.get_buffer_wrapper(&density_viscosity_temp_quasi_layer<1>::get_temperature_wrapper);

                // Шаг по вязкости
                quickest_ultimate_fv_solver solver_tm(*heatModel, temp_wrapper);
                solver_tm.step(dt, boundaries.temperature, boundaries.temperature);
            }
            else {
                buffer.current().temperature = vector<double>(buffer.current().temperature.size(), boundaries.temperature);
            }
        }
    }

    /// @brief 
    /// @param dt 
    /// @param boundaries 
    void make_rheology_step_shukhov(double dt, const nonisothermal_quasistatic_PQ_task_boundaries_t& boundaries) {

        vector<double>Q_profile(pipe.profile.get_point_count(), boundaries.volumetric_flow); /// задаем по трубе новый расход из временного ряда
        vector<double> G(pipe.profile.get_point_count(), Q_profile[0] * oil.density.nominal_density);   /// массовый расход
        //pipe.heat.ambient_heat_transfer = 1.4917523388199689;
        PipeHeatInflowConstArea heatModel(pipe, oil, G);         /// один хитмодел для квикеста, второй для Шухова

        advance(); // Сдвигаем текущий и предыдущий слои

        /// Солвер для Шухова
        solve_euler_corrector<1>(heatModel, +1, { boundaries.temperature }, &buffer.current().temperature);
    }

public:
    /// @brief Рассчёт шага моделирования, включающий в себя расчёт шага движения партии и гидравлический расчёт
    /// Функция делат сдвиг буфера (advance) так, что buffer.current после вызова содержит свежерасчитанный слой
    /// @param dt временной шаг моделирования
    /// @param boundaries Краевые условие
    void step(double dt, nonisothermal_quasistatic_PQ_task_boundaries_t& boundaries, step_mode_t mode) {
        
        switch (mode) {
        case step_mode_t::Advection:
            pipe.heat.ambient_heat_transfer = 1.3786917741689342;
            make_rheology_step_advection(dt, boundaries);
            break;

        case step_mode_t::Shukhov:
            pipe.heat.ambient_heat_transfer = 1.4917523388199689;
            make_rheology_step_shukhov(dt, boundaries);
            break;

        case step_mode_t::ShukhovWithAdvection:
            pipe.heat.ambient_heat_transfer = 1.4917523388199689;
            make_rheology_step_shukhov(dt, boundaries);
            boundaries.temperature = buffer.current().temperature.back();
            make_rheology_step_advection(dt, boundaries);
            break;
        }

         //////////
    }

    /// @brief Сдвиг текущего слоя в буфере
    void advance()
    {
        buffer.advance(+1);
    }

    /// @brief Возвращает ссылку на буфер
    auto& get_buffer()
    {
        return buffer;
    }
};


/// @brief Стационарный расчет (с помощью initial boundaries),
/// а затем квазистационарный расчет по краевым условиям (boundary_timeseries)
/// @tparam Solver Численный метод расчета движения партий
/// @tparam Printer Класс для вывода результатов в файл
/// @param path Путь к файлам с результатами
/// @param pipe Модель трубы
/// @param initial_boundaries Начальные условия
/// @param boundary_timeseries Краевые условия
/// !!! Важно, чтобы вектор на заданный момент времени был совместим по порядку параметров с isothermal_quasistatic_task_boundaries_t !!!
/// @param model_type Способ расчёта
/// @param etalon_timeseries Эталонные данные давления в конце трубопровода 
/// @param dt Шаг по времени либо задаётся постоянным, 
/// либо рассчитывается на каждом шаге моделирования для Cr = 1 
template <typename Solver, typename Printer>
inline void perform_noniso_quasistatic_simulation(
    const string& path,
    const pipe_noniso_properties_t& pipe,
    const oil_parameters_t& oil,
    const nonisothermal_quasistatic_PQ_task_boundaries_t& initial_boundaries,
    const vector_timeseries_t& boundary_timeseries,
    const noniso_qsm_model_type& model_type,
    const vector_timeseries_t& etalon_timeseries,
    typename nonisothermal_quasistatic_PQ_task_t<Solver>::step_mode_t step_mode,
    double dt
)
{
    time_t t = boundary_timeseries.get_start_date(); // Момент времени начала моделирования
    Printer layer_printer;

    nonisothermal_quasistatic_PQ_task_t<Solver> task(pipe, oil, model_type);
    task.solve(initial_boundaries);

    // Печатаем профиль трубы и первый слой к нему
    write_profile(pipe.profile, path + "pipe_coord_heights");
    // Вывод начального расчёта
    layer_printer.print_t_all(path, t, pipe, task.get_current_layer());

    do
    {
        // Интерполируем значения параметров в заданный момент времени
        vector<double> values_in_time_model = boundary_timeseries(t);
        nonisothermal_quasistatic_PQ_task_boundaries_t boundaries(values_in_time_model);
        
        double time_step = dt;
        if (std::isnan(time_step)) {
            const auto& vec = boundary_timeseries.data[0].second;
            auto max_it = std::max_element(vec.begin(), vec.end());
            double v = *max_it / pipe.wall.getArea();       // Cr превышается, костыль                    
            //double v = boundaries.volumetric_flow / pipe.wall.getArea();      
            time_step = task.get_time_step_assuming_max_speed(v);
        }
        t += static_cast<time_t>(time_step);

        // Делаем шаг
        task.step(time_step, boundaries, step_mode);

        // Вывод профилей и временного ряда сравнения с эталонными данными
        if (etalon_timeseries.data.empty())
        {

            layer_printer.print_t_all(path, t, pipe, task.get_current_layer());
        }
        else {

            layer_printer.print_t_all(path, t, pipe, task.get_current_layer(), etalon_timeseries(t), boundaries);
        }

    } while (t <= boundary_timeseries.get_end_date());

}

/// @brief Перегрузка функции для возможности 
/// не передавать начальные условия
template <typename Solver, typename Printer>
inline void perform_noniso_quasistatic_simulation(
    const string& path,
    const pipe_noniso_properties_t& pipe,
    const oil_parameters_t& oil,
    const vector_timeseries_t& boundary_timeseries,
    const noniso_qsm_model_type& model_type,
    const vector_timeseries_t& etalon_timeseries,
    typename nonisothermal_quasistatic_PQ_task_t<Solver>::step_mode_t step_mode,
    double dt=std::numeric_limits<double>::quiet_NaN()
)
{
    time_t t = boundary_timeseries.get_start_date(); // Момент времени начала моделирования
    nonisothermal_quasistatic_PQ_task_boundaries_t initial_boundaries(boundary_timeseries(t));

    perform_noniso_quasistatic_simulation<Solver, Printer>(path, pipe, oil, initial_boundaries, boundary_timeseries, model_type, etalon_timeseries, step_mode, dt);
}


/// @brief Накапливает результаты по выходной температуре
class nonisothermal_qsm_batch_Tout_collector_t
    : public batch_processor_precalculated_times<density_viscosity_temp_quasi_layer<true>>
{
public:
    /// @brief Тип данных слоя буфера изотермической квазистационарной модели
    typedef density_viscosity_temp_quasi_layer<true> layer_type;
protected:
    /// @brief Вектор расчётных значений давления на выходе ЛУ
    vector<double> pipe_temperature_out;
public:
    /// @brief Конструктор обработчика
    /// @param times Предпосчитанная временная сетка моделирования работы ЛУ
    nonisothermal_qsm_batch_Tout_collector_t(const vector<double>& times)
        : pipe_temperature_out(times.size(), std::numeric_limits<double>::quiet_NaN())
    {

    }

    /// @brief Сохранение результатов расчёта давления в конце ЛУ в вектор
    /// @param step_index Текущий шаг моделирования
    /// @param layer Текущий слой
    virtual void process_data(size_t step_index,
        const density_viscosity_temp_quasi_layer<true>& layer) override
    {
        // at() - проверяет выход за границы массива
        //pipe_pressure_out.at(step_index) = layer.pressure.back();
        pipe_temperature_out[step_index] = layer.temperature.back();
    }
    /// @brief Геттер для вектора собранных результатов расчёта давления в конце ЛУ
    /// @return Вектор расчётных значений давления на выходе ЛУ
    const vector<double>& get_temp_out_calculated() const {
        return pipe_temperature_out;
    }
};

/// @brief Пакетный изотермический квазистатический расчет с предподсчитанным временем
/// делает статический расчет task.solve, 
/// а затем столько раз task.step, сколько временных меток в times
/// @tparam Task Тип расчетной задачи
/// @tparam DataProcessor Обработчик данных каждого расчетного слоя
/// @param task 
/// @param times 
/// @param boundary_timeseries 
/// @param data_processor 
template <typename Task, typename DataProcessor>
inline void quasistatic_batch(
    Task& task,
    const std::vector<double>& times,
    const std::vector<std::vector<double>>& boundary_timeseries,
    DataProcessor* data_processor,
    typename Task::step_mode_t step_mode
)
{
    // Вычленение начальных условий
    typename Task::boundaries_type initial_boundaries(boundary_timeseries[0]);
    task.solve(initial_boundaries);
    data_processor->process_data(0, task.get_buffer().current());

    for (size_t step_index = 1; step_index < times.size(); step_index++)
    {
        double time_step = times[step_index] - times[step_index - 1];
        typename Task::boundaries_type boundaries(boundary_timeseries[step_index]);

        task.step(time_step, boundaries, step_mode);

        data_processor->process_data(step_index, task.get_buffer().current());
    }
};

}
