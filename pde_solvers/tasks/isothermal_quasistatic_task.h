#pragma once

namespace pde_solvers {
;

/// @brief Варианты расчёта квазистационарной модели
/// Stationary - стационарный расчёт (частный случай квазистаца)
/// DensityQuasi - квазистац по плотности, вязкость рассчитывается как в стационаре
/// ViscosityQuasi - квазистац по вязкости, плотность рассчитывается как в стационаре
/// FullQuasi - полный квазистац
enum class QuasistaticModelType {
    Stationary, DensityQuasi, ViscosityQuasi, FullQuasi, TempQuasi
};

/// @brief Проблемно-ориентированный слой для гидравлического квазистационарного расчета
/// @tparam CellFlag Флаг расчёта реологии 
/// true - в ячейках для метода конечных объёмов (Quickest-Ultimate)
/// false - в точках для метода характеристик (advection_moc_solver)
template <bool CellFlag>
struct density_viscosity_quasi_layer {
    /// @brief Профиль плотности
    std::vector<double> density;
    /// @brief Профиль вязкости
    std::vector<double> viscosity;
    /// @brief Профиль давления
    std::vector<double> pressure;
    /// @brief Дифференциальный профиль давления
    std::vector<double> pressure_delta;
    /// @brief Изначальный профиль давления
    std::vector<double> pressure_initial;
    /// @brief Профиль вспомогательных расчетов для метода конечных объемов (и для вязкости, и для плотности)
    quickest_ultimate_fv_solver_traits<1>::specific_layer specific;
    /// @brief Инициализация профилей
    /// @param point_count Количество точек
    density_viscosity_quasi_layer(size_t point_count)
        : density(point_count - static_cast<int>(CellFlag))
        , viscosity(point_count - static_cast<int>(CellFlag))
        , specific(point_count)
        , pressure(point_count)
        , pressure_delta(point_count)
    {}

    /// @brief Подготовка плотности для расчета методом конечных объемов 
    /// @param layer Слой
    /// @return Обертка над составным слоем
    static quickest_ultimate_fv_wrapper<1> get_density_wrapper(density_viscosity_quasi_layer& layer)
    {
        return quickest_ultimate_fv_wrapper<1>(layer.density, layer.specific);
    }

    /// @brief Подготовка вязкости для расчета методом конечных объемов
    /// @param layer Слой
    /// @return Обертка над составным слоем
    static quickest_ultimate_fv_wrapper<1> get_viscosity_wrapper(density_viscosity_quasi_layer& layer)
    {
        return quickest_ultimate_fv_wrapper<1>(layer.viscosity, layer.specific);
    }
};

/// @brief Структура, содержащая в себе краевые условия задачи PQ
struct isothermal_quasistatic_PQ_task_boundaries_t {
    /// @brief Изначальный объемный расход
    double volumetric_flow;
    /// @brief Изначальное давление на входе
    double pressure_in;
    /// @brief Изначальная плотность на входе
    double density;
    /// @brief Изначальная вязкость на входе
    double viscosity;

    /// @brief Конструктор по умолчанию
    isothermal_quasistatic_PQ_task_boundaries_t() = default;

    /// @brief Конструктор краевых условий
    /// @param values Значения краевых условий
    isothermal_quasistatic_PQ_task_boundaries_t(const std::vector<double>& values) {
        volumetric_flow = values[0];
        pressure_in = values[1];
        density = values[2];
        viscosity = values[3];
    }

    /// @brief Создание структуры со значениями по умолчанию
    static isothermal_quasistatic_PQ_task_boundaries_t default_values() {
        isothermal_quasistatic_PQ_task_boundaries_t result;
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
/// @tparam Solver Тип солвера партий (advection_moc_solver или quickest_ultimate_fv_solver)
template <typename Solver=advection_moc_solver>
class isothermal_quasistatic_PQ_task_t {
public:
    /// @brief Тип слоя
    using layer_type = density_viscosity_quasi_layer<std::is_same<Solver, advection_moc_solver>::value ? false : true>;
    /// @brief Тип буфера
    using buffer_type = ring_buffer_t<layer_type>;
    /// @brief Тип граничных условий
    using boundaries_type = isothermal_quasistatic_PQ_task_boundaries_t;
private:
    // Модель трубы
    pipe_properties_t pipe;
    // Создаётся буфер, тип слоя которого определяется в зависимости от типа солвера
    buffer_type buffer;
    // Тип расчёта
    QuasistaticModelType model_type;

public:
    /// @brief Конструктор
    /// @param pipe Модель трубопровода
    isothermal_quasistatic_PQ_task_t(const pipe_properties_t& pipe, QuasistaticModelType model_type = QuasistaticModelType::FullQuasi)
        : pipe(pipe)
        , buffer(2, pipe.profile.get_point_count())
        , model_type( model_type )
    {
    }

    /// @brief Геттер для текущего слоя  
    density_viscosity_quasi_layer<std::is_same<Solver, advection_moc_solver>::value ? false : true> get_current_layer() {
        return buffer.current();
    }

    /// @brief Начальный стационарный расчёт. 
    /// Ставим по всей трубе реологию из initial_conditions, делаем гидравлический расчет
    /// @param initial_conditions Начальные условия
    void solve(const isothermal_quasistatic_PQ_task_boundaries_t& initial_conditions)
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

        //// Начальный гидравлический расчет
        calc_pressure_layer(initial_conditions);
        buffer.previous().pressure_initial = current.pressure_initial = current.pressure; // Получаем изначальный профиль давлений
    }
public:
    /// @brief Расчёт шага по времени для Cr = 1
    /// @param v_max Максимальная скорость течение потока в трубопроводе
    double get_time_step_assuming_max_speed(double v_max) const {
        const auto& x = pipe.profile.coordinates;
        double dx = x[1] - x[0]; // Шаг сетки
        double dt = abs(dx / v_max); // Постоянный шаг по времени для Куранта = 1
        return dt;
    }
private:
    /// @brief Проводится расчёт шага движения партии
    /// @param dt Временной шаг моделирования
    /// @param boundaries Краевые условия
    void make_rheology_step(double dt, const isothermal_quasistatic_PQ_task_boundaries_t& boundaries) {
        size_t n = pipe.profile.get_point_count();
        std::vector<double>Q_profile(n, boundaries.volumetric_flow); // задаем по трубе новый расход из временного ряда

        advance(); // Сдвигаем текущий и предыдущий слои

        if constexpr (std::is_same<Solver, advection_moc_solver>::value) {
            // считаем партии методом характеристик

            // плотность - квазистац или стац
            if (model_type == QuasistaticModelType::FullQuasi || model_type == QuasistaticModelType::DensityQuasi) {
                // Шаг по плотности
                advection_moc_solver solver_rho(pipe, Q_profile[0], buffer.previous().density, buffer.current().density);
                solver_rho.step(dt, boundaries.density, boundaries.density);
            } else
                buffer.current().density = std::vector<double>(buffer.current().density.size(), boundaries.density);

            // вязкость - квазистац или стац
            if (model_type == QuasistaticModelType::FullQuasi || model_type == QuasistaticModelType::ViscosityQuasi) {
                // Шаг по вязкости
                advection_moc_solver solver_nu(pipe, Q_profile[0], buffer.previous().viscosity, buffer.current().viscosity);
                solver_nu.step(dt, boundaries.viscosity, boundaries.viscosity);
            } else
                buffer.current().viscosity = std::vector<double>(buffer.current().viscosity.size(), boundaries.viscosity);
        }
        else {
            // считаем партии с помощью QUICKEST-ULTIMATE
            PipeQAdvection advection_model(pipe, Q_profile);

            // плотность - квазистац или стац
            if (model_type == QuasistaticModelType::FullQuasi || model_type == QuasistaticModelType::DensityQuasi) {
                // Шаг по плотности
                auto density_wrapper = buffer.get_buffer_wrapper(
                    &density_viscosity_quasi_layer<1>::get_density_wrapper);
                quickest_ultimate_fv_solver solver_rho(advection_model, density_wrapper);
                solver_rho.step(dt, boundaries.density, boundaries.density);
            }
            else {
                buffer.current().density = std::vector<double>(buffer.current().density.size(), boundaries.density);
            }

            // вязкость - квазистац или стац
            if (model_type == QuasistaticModelType::FullQuasi || model_type == QuasistaticModelType::ViscosityQuasi) {
                auto viscosity_wrapper = buffer.get_buffer_wrapper(
                    &density_viscosity_quasi_layer<1>::get_viscosity_wrapper);
                
                // Шаг по вязкости
                quickest_ultimate_fv_solver solver_nu(advection_model, viscosity_wrapper);
                solver_nu.step(dt, boundaries.viscosity, boundaries.viscosity);
            } else
                buffer.current().viscosity = std::vector<double>(buffer.current().viscosity.size(), boundaries.viscosity);
        }
    }

    /// @brief Рассчёт профиля давления методом Эйлера (задача PQ)
    /// @param boundaries Краевые условия
    void calc_pressure_layer(const isothermal_quasistatic_PQ_task_boundaries_t& boundaries) {

        auto& current = buffer.current();

        std::vector<double>& p_profile = current.pressure;
        int euler_direction = +1; // Задаем направление для Эйлера

        isothermal_pipe_PQ_parties_t pipeModel(pipe, current.density, current.viscosity, boundaries.volumetric_flow, euler_direction);
        solve_euler<1>(pipeModel, euler_direction, boundaries.pressure_in, &p_profile);
        // Получаем дифференциальный профиль давлений
        std::transform(current.pressure_initial.begin(), current.pressure_initial.end(), p_profile.begin(),
            current.pressure_delta.begin(),
            [](double initial, double current) {return initial - current;  });

    }
public:
    /// @brief Рассчёт шага моделирования, включающий в себя расчёт шага движения партии и гидравлический расчёт
    /// Функция делат сдвиг буфера (advance) так, что buffer.current после вызова содержит свежерасчитанный слой
    /// @param dt временной шаг моделирования
    /// @param boundaries Краевые условие
    void step(double dt, const isothermal_quasistatic_PQ_task_boundaries_t& boundaries) {
        make_rheology_step(dt, boundaries);
        calc_pressure_layer(boundaries);
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
inline void perform_quasistatic_simulation(
    const std::string& path,
    const pipe_properties_t& pipe,
    const isothermal_quasistatic_PQ_task_boundaries_t& initial_boundaries,
    const vector_timeseries_t& boundary_timeseries,
    const QuasistaticModelType& model_type,
    const vector_timeseries_t& etalon_timeseries,
    double dt
)
{
    time_t t = boundary_timeseries.get_start_date(); // Момент времени начала моделирования
    Printer layer_printer;

    isothermal_quasistatic_PQ_task_t<Solver> task(pipe, model_type);
    task.solve(initial_boundaries);


    // Печатаем профиль трубы и первый слой к нему
    write_profile(pipe.profile, path + "pipe_coord_heights");
    // Вывод начального расчёта
    layer_printer.print_all(path, t, pipe, task.get_current_layer());

    do
    {
        // Интерполируем значения параметров в заданный момент времени
        std::vector<double> values_in_time_model = boundary_timeseries(t);
        isothermal_quasistatic_PQ_task_boundaries_t boundaries(values_in_time_model);

        double time_step = dt;
        if (std::isnan(time_step)) {
            double v = boundaries.volumetric_flow / pipe.wall.getArea();
            time_step = task.get_time_step_assuming_max_speed(v);
        }
        t += static_cast<time_t>(time_step);

        // Делаем шаг
        task.step(time_step, boundaries);

        // Вывод профилей и временного ряда сравнения с эталонными данными
        if (etalon_timeseries.data.empty())
        {

            layer_printer.print_all(path, t, pipe, task.get_current_layer());
        }
        else {

            layer_printer.print_all(path, t, pipe, task.get_current_layer(), etalon_timeseries(t));
        }

    } while (t <= boundary_timeseries.get_end_date());

}

/// @brief Перегрузка функции для возможности 
/// не передавать начальные условия
template <typename Solver, typename Printer>
inline void perform_quasistatic_simulation(
    const std::string& path,
    const pipe_properties_t& pipe,
    const vector_timeseries_t& boundary_timeseries,
    const QuasistaticModelType& model_type,
    const vector_timeseries_t& etalon_timeseries,
    double dt=std::numeric_limits<double>::quiet_NaN())
{

    time_t t = boundary_timeseries.get_start_date(); // Момент времени начала моделирования
    isothermal_quasistatic_PQ_task_boundaries_t initial_boundaries(boundary_timeseries(t));

    perform_quasistatic_simulation<Solver, Printer>(path, pipe, initial_boundaries, boundary_timeseries, model_type, etalon_timeseries, dt);
}

/// @brief Перегрузка функции для возможности задания постоянного
/// шага по времени без эталонных данных
template <typename Solver, typename Printer>
inline void perform_quasistatic_simulation(
    const std::string& path,
    const pipe_properties_t& pipe,
    const isothermal_quasistatic_PQ_task_boundaries_t& initial_boundaries,
    const vector_timeseries_t& boundary_timeseries,
    const QuasistaticModelType& model_type,
    double dt=std::numeric_limits<double>::quiet_NaN())
{
    perform_quasistatic_simulation<Solver, Printer>(path, pipe, initial_boundaries, boundary_timeseries, model_type, vector_timeseries_t({}), dt);
}




/// @brief Обработчик результатов расчета при пакетном расчете с предподсчитанными временными метками
/// @tparam LayerType Тип слоя (под МХ, под QUICKEST)
template <typename LayerType>
class batch_processor_precalculated_times {
public:
    /// @brief Шаблон для задания функции обработки результатов расчёта
    /// @param step_index Текущий шаг моделирования
    /// @param layer Текущий слой
    virtual void process_data(size_t step_index, const LayerType& layer) = 0;
};

/// @brief Накапливает результаты по выходному давлению
class isothermal_qsm_batch_Pout_collector_t
    : public batch_processor_precalculated_times<density_viscosity_quasi_layer<true>>
{
public:
    /// @brief Тип данных слоя буфера изотермической квазистационарной модели
    typedef density_viscosity_quasi_layer<true> layer_type;
protected:
    /// @brief Вектор расчётных значений давления на выходе ЛУ
    std::vector<double> pipe_pressure_out;
public:
    /// @brief Конструктор обработчика
    /// @param times Предпосчитанная временная сетка моделирования работы ЛУ
    isothermal_qsm_batch_Pout_collector_t(const std::vector<double>& times)
        : pipe_pressure_out(times.size(), std::numeric_limits<double>::quiet_NaN())
    {

    }

    /// @brief Сохранение результатов расчёта давления в конце ЛУ в вектор
    /// @param step_index Текущий шаг моделирования
    /// @param layer Текущий слой
    virtual void process_data(size_t step_index,
        const density_viscosity_quasi_layer<true>& layer) override
    {
        // at() - проверяет выход за границы массива
        //pipe_pressure_out.at(step_index) = layer.pressure.back();
        pipe_pressure_out[step_index] = layer.pressure.back();
    }
    /// @brief Геттер для вектора собранных результатов расчёта давления в конце ЛУ
    /// @return Вектор расчётных значений давления на выходе ЛУ
    const std::vector<double>& get_pressure_out_calculated() const {
        return pipe_pressure_out;
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
    DataProcessor* data_processor
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

        task.step(time_step, boundaries);

        data_processor->process_data(step_index, task.get_buffer().current());
    }
};


}
