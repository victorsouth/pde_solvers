#pragma once

namespace pde_solvers {
;

/// @brief Варианты беспартийного неизотермического расчёта 
enum class noniso_qsm_model_type {
    /// @brief Полноценная динамическая модель
    Dynamic,            
    /// @brief Формула Шухова
    Shukhov,              
    /// @brief Шухов c адвекцией
    ShukhovWithAdvection,  
    /// @brief Температура по всей трубе равна входной
    Isothemal 
};

/// @brief Проблемно-ориентированный слой для TQ расчета
/// @tparam CellFlag Флаг расчёта реологии 
/// true - в ячейках для метода конечных объёмов (Quickest-Ultimate)
/// false - в точках для метода характеристик (advection_moc_solver)
struct qsm_noniso_T_layer {
    /// @brief Профиль температуры (ячейки)
    std::vector<double> temperature;
    /// @brief Профиль температуры  Шуховым (точки)
    std::vector<double> temperature_shukhov;
    /// @brief Профиль вспомогательных расчетов для метода конечных объемов (и для вязкости, и для плотности)
    quickest_ultimate_fv_solver_traits<1>::specific_layer specific;
    /// @brief Инициализация профилей
    /// @param point_count Количество точек
    qsm_noniso_T_layer(size_t point_count)
        : temperature(point_count - 1)
        , temperature_shukhov(point_count)
        , specific(point_count)
    {
    }

    /// @brief Подготовка температуры для расчета методом конечных объемов
    /// @param layer Слой
    /// @return Обертка над составным слоем
    static quickest_ultimate_fv_wrapper<1> get_temperature_wrapper(qsm_noniso_T_layer& layer)
    {
        return quickest_ultimate_fv_wrapper<1>(layer.temperature, layer.specific);
    }
};

/// @brief Структура, содержащая в себе краевые условия задачи T?
struct qsm_noniso_T_task_boundaries_t { 
    /// @brief Изначальный объемный расход
    double volumetric_flow;
    /// @brief Изначальное температура на входе
    double temperature;
    /// @brief Изначальная темп на входе
    double temperature_shukhov;
    /// @brief Конструктор по умолчанию
    qsm_noniso_T_task_boundaries_t() = default;

    /// @brief Конструктор краевых условий
    /// @param values Значения краевых условий
    qsm_noniso_T_task_boundaries_t(const std::vector<double>& values) {
        volumetric_flow = values[0];
        temperature = values[1];
        temperature_shukhov = values[1];
    }

    /// @brief Создание структуры со значениями по умолчанию
    static qsm_noniso_T_task_boundaries_t default_values() {
        qsm_noniso_T_task_boundaries_t result;
        result.volumetric_flow = 0.2;
        result.temperature = 300;
        result.temperature_shukhov = 300;
        return result;
    }
};


/// @brief Расчетная задача (task) для гидравлического неизотермического 
/// квазистационарного расчета в условиях движения партий с разной плотностью и вязкостью, температурой
/// Расчет партий делается методом характеристик или Quickest-Ultimate
/// @tparam Solver Тип солвера партий (advection_moc_solver или quickest_ultimate_fv_solver)
class qsm_noniso_T_task_t {  
public:
    /// @brief Тип слоя
    using layer_type = qsm_noniso_T_layer;
    /// @brief Тип буфера
    using buffer_type = ring_buffer_t<layer_type>;
    /// @brief Тип граничных условий
    using boundaries_type = qsm_noniso_T_task_boundaries_t;

    
private:
    // Модель трубы
    pipe_noniso_properties_t pipe;
    // Нефть
    oil_parameters_t oil;
    // Создаётся буфер, тип слоя которого определяется в зависимости от типа солвера
    buffer_type buffer;
    /// @brief Тепловая модель трубы
    noniso_qsm_model_type model_type;

public:
    /// @brief Конструктор
    /// @param pipe Модель трубопровода
    qsm_noniso_T_task_t(const pipe_noniso_properties_t& pipe, 
        const oil_parameters_t& oil, noniso_qsm_model_type model_type)
        : pipe(pipe)
        , oil(oil)
        , buffer(2, pipe.profile.get_point_count()) //2
        , model_type(model_type)
    {
    }

    /// @brief Геттер для текущего слоя  
    qsm_noniso_T_layer& get_current_layer() {
        return buffer.current();
    }

    /// @brief Начальный стационарный расчёт. 
    /// Ставим по всей трубе реологию из initial_conditions, делаем гидравлический расчет
    /// @param initial_conditions Начальные условия
    void solve(const qsm_noniso_T_task_boundaries_t& initial_conditions)
    {
        // Количество точек
        size_t n = pipe.profile.get_point_count();

        // Инициализация реологии
        auto& current = buffer.current();

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
    /// @brief 
    /// @param dt 
    /// @param boundaries 
    void make_rheology_step_shukhov(double dt, const qsm_noniso_T_task_boundaries_t& boundaries) {

        std::vector<double>Q_profile(pipe.profile.get_point_count(), boundaries.volumetric_flow); /// задаем по трубе новый расход из временного ряда
        std::vector<double> G(pipe.profile.get_point_count(), Q_profile[0] * oil.density.nominal_density);   /// массовый расход
        //pipe.heat.ambient_heat_transfer = 1.4917523388199689;
        PipeHeatInflowConstArea heatModel(pipe, oil, G);         /// один хитмодел для квикеста, второй для Шухова

        advance(); // Сдвигаем текущий и предыдущий слои

    }

    /// @brief Полноценный шаг по температуре. Выполняется в T_quasi_layer::temperature
    void step_dynamic_temperature(double dt, const qsm_noniso_T_task_boundaries_t& boundaries) {
        size_t n = pipe.profile.get_point_count();
        //std::vector<double> G(n-1, boundaries.volumetric_flow * oil.density.nominal_density);   /// массовый расход
        std::vector<double> G(n, boundaries.volumetric_flow * oil.density.nominal_density);   /// массовый расход
        PipeHeatInflowConstArea heatModel(pipe, oil, G);
        auto temperature_wrapper = buffer.get_buffer_wrapper(&qsm_noniso_T_layer::get_temperature_wrapper);
        quickest_ultimate_fv_solver solver_tm(heatModel, temperature_wrapper);
        solver_tm.step(dt, boundaries.temperature, boundaries.temperature);
    }
    /// @brief Шаг по температуре в виде адвекции. 
    /// Подразумевается, что temperature_in - это температура, рассчитанная по Шухову
    void step_advection_temperature(double dt, double temperature_in, double vol_flow) {
        // считаем партии с помощью QUICKEST-ULTIMATE
        size_t n = pipe.profile.get_point_count();
        std::vector<double> Q_profile(n-1, vol_flow); /// задаем по трубе новый расход из временного ряда
        PipeQAdvection advection_model(pipe, Q_profile);

        auto temperature_wrapper = buffer.get_buffer_wrapper(&qsm_noniso_T_layer::get_temperature_wrapper);

        quickest_ultimate_fv_solver solver_rho(advection_model, temperature_wrapper);
        solver_rho.step(dt, temperature_in, temperature_in);

    }
    /// @brief Расчет по Шухову. Возвращает температуру в конце трубы
    double solve_shukhov_temperature(const qsm_noniso_T_task_boundaries_t& boundaries) {
        size_t n = pipe.profile.get_point_count();
        std::vector<double> G(n, boundaries.volumetric_flow * oil.density.nominal_density);   /// массовый расход
        PipeHeatInflowConstArea heatModel(pipe, oil, G);
        solve_euler_corrector<1>(heatModel, +1, { boundaries.temperature }, &buffer.current().temperature_shukhov);
        double Tout = buffer.current().temperature_shukhov.back();
        return Tout;

    }

public:
    /// @brief Рассчёт шага моделирования, включающий в себя расчёт шага движения партии и гидравлический расчёт
    /// Функция делат сдвиг буфера (advance) так, что buffer.current после вызова содержит свежерасчитанный слой
    /// @param dt временной шаг моделирования
    /// @param boundaries Краевые условие
    void step(double dt, qsm_noniso_T_task_boundaries_t& boundaries) {
        advance(); // Сдвигаем текущий и предыдущий слои


        switch (model_type) {
        case noniso_qsm_model_type::Dynamic: {
            // TODO: здесь нельзя это прописывать, должно быть в исходных данных!!!
            // см. NonisothermalQuasistaticModelWithRealData, DynamicTemperature
            //pipe.heat.ambient_heat_transfer = 1.3786917741689342;
            step_dynamic_temperature(dt, boundaries);
            break;
        }
        case noniso_qsm_model_type::Shukhov: {
            //pipe.heat.ambient_heat_transfer = 1.4917523388199689;
            double T_out_shukhov = solve_shukhov_temperature(boundaries);
            auto& T_layer = buffer.current().temperature;
            std::fill(T_layer.begin(), T_layer.end(), T_out_shukhov); // тупо копируем температур по всему слою
            break;
        }

        case noniso_qsm_model_type::ShukhovWithAdvection: {
            //pipe.heat.ambient_heat_transfer = 1.4917523388199689;
            double T_out_shukhov = solve_shukhov_temperature(boundaries);
            // запускаем адвекцию, на вход которой даем температуру с выхода, посчитанную по Шухову
            step_advection_temperature(dt, T_out_shukhov, boundaries.volumetric_flow); 
            break;
        }
        }
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
//template <typename Solver, typename Printer>
//inline void perform_noniso_quasistatic_simulation(
//    const std::string& path,
//    const pipe_noniso_properties_t& pipe,
//    const oil_parameters_t& oil,
//    const nonisothermal_quasistatic_PQ_task_boundaries_t& initial_boundaries,
//    const vector_timeseries_t& boundary_timeseries,
//    const noniso_qsm_model_type& model_type,
//    const vector_timeseries_t& etalon_timeseries,
//    noniso_qsm_model_type step_mode,
//    double dt
//)
//{
//    time_t t = boundary_timeseries.get_start_date(); // Момент времени начала моделирования
//    Printer layer_printer;
//
//    noniso_quasistatic_PQ_task_t<Solver> task(pipe, oil, model_type);
//    task.solve(initial_boundaries);
//
//    // Печатаем профиль трубы и первый слой к нему
//    write_profile(pipe.profile, path + "pipe_coord_heights");
//    // Вывод начального расчёта
//    layer_printer.print_t_all(path, t, pipe, task.get_current_layer());
//
//    do
//    {
//        // Интерполируем значения параметров в заданный момент времени
//        std::vector<double> values_in_time_model = boundary_timeseries(t);
//        nonisothermal_quasistatic_PQ_task_boundaries_t boundaries(values_in_time_model);
//        
//        double time_step = dt;
//        if (std::isnan(time_step)) {
//            const auto& vec = boundary_timeseries.data[0].second;
//            auto max_it = std::max_element(vec.begin(), vec.end());
//            double v = *max_it / pipe.wall.getArea();       // Cr превышается, костыль                    
//            //double v = boundaries.volumetric_flow / pipe.wall.getArea();      
//            time_step = task.get_time_step_assuming_max_speed(v);
//        }
//        t += static_cast<time_t>(time_step);
//
//        // Делаем шаг
//        task.step(time_step, boundaries, step_mode);
//
//        // Вывод профилей и временного ряда сравнения с эталонными данными
//        if (etalon_timeseries.data.empty())
//        {
//
//            layer_printer.print_t_all(path, t, pipe, task.get_current_layer());
//        }
//        else {
//
//            layer_printer.print_t_all(path, t, pipe, task.get_current_layer(), etalon_timeseries(t), boundaries);
//        }
//
//    } while (t <= boundary_timeseries.get_end_date());
//
//}
//
///// @brief Перегрузка функции для возможности 
///// не передавать начальные условия
//template <typename Solver, typename Printer>
//inline void perform_noniso_quasistatic_simulation(
//    const std::string& path,
//    const pipe_noniso_properties_t& pipe,
//    const oil_parameters_t& oil,
//    const vector_timeseries_t& boundary_timeseries,
//    const noniso_qsm_model_type& model_type,
//    const vector_timeseries_t& etalon_timeseries,
//    noniso_qsm_model_type step_mode,
//    double dt=std::numeric_limits<double>::quiet_NaN()
//)
//{
//    time_t t = boundary_timeseries.get_start_date(); // Момент времени начала моделирования
//    nonisothermal_quasistatic_PQ_task_boundaries_t initial_boundaries(boundary_timeseries(t));
//
//    perform_noniso_quasistatic_simulation<Solver, Printer>(path, pipe, oil, initial_boundaries, boundary_timeseries, model_type, etalon_timeseries, step_mode, dt);
//}


/// @brief Накапливает результаты по выходной температуре
class nonisothermal_qsm_batch_Tout_collector_t
    : public batch_processor_precalculated_times<qsm_noniso_T_layer>
{
public:
    /// @brief Тип данных слоя буфера изотермической квазистационарной модели
    typedef qsm_noniso_T_layer layer_type;
protected:
    /// @brief Вектор расчётных значений давления на выходе ЛУ
    std::vector<double> pipe_temperature_out;
public:
    /// @brief Конструктор обработчика
    /// @param times Предпосчитанная временная сетка моделирования работы ЛУ
    nonisothermal_qsm_batch_Tout_collector_t(const std::vector<double>& times)
        : pipe_temperature_out(times.size(), std::numeric_limits<double>::quiet_NaN())
    {

    }

    /// @brief Сохранение результатов расчёта давления в конце ЛУ в вектор
    /// @param step_index Текущий шаг моделирования
    /// @param layer Текущий слой
    virtual void process_data(size_t step_index,
        const qsm_noniso_T_layer& layer) override
    {
        // at() - проверяет выход за границы массива
        //pipe_pressure_out.at(step_index) = layer.pressure.back();
        pipe_temperature_out[step_index] = layer.temperature.back();
    }
    /// @brief Геттер для вектора собранных результатов расчёта давления в конце ЛУ
    /// @return Вектор расчётных значений давления на выходе ЛУ
    const std::vector<double>& get_temp_out_calculated() const {
        return pipe_temperature_out;
    }
};


}
