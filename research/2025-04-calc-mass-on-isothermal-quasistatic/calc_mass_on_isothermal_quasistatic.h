#pragma once
// Среднее значение коэффициента сжимаемости нефти 
#define AVG_OIL_COMPRESSIBILITY_COEFF 0.00078e-6 // МПа-1 - надо ли переводить?

class CalcMassQSM : public ::testing::Test {

protected:
    // Модель трубы
    pipe_properties_t pipe;
    // Краевые условия по плотности
    vector<double> density;
    // Краевые условия по давлению в начале ЛУ
    vector<double> pressure_in;
    // Краевые условия по давлению в конце ЛУ
    vector<double> pressure_out;
    // Краевые условия по объёмному расходу
    vector<double> Q;
    // Временная сетка моделирования
    vector<double> times;

    /// @brief Подготовка к расчету
    virtual void SetUp() override {
        double Pcapacity = 10e6;

        double L = 100e3; // 100 км
        double dx = 200;
        size_t n = static_cast<size_t>(L / dx + 0.000001);

        pipe.profile = pipe_profile_t::create(n, 0, L, 0, 0, Pcapacity);
        pipe.wall.diameter = 0.7;
        pipe.wall.wallThickness = 0.01;

        // Значения плотности на входе трубы
        density =
        {
            900,
            880,
            880,
            890,
            890,
            880,
            880,
            870
        };

        //значения давления на входе трубы
        pressure_in =
        {
            6e6,
            5.8e6,
            5.8e6,
            5.9e6,
            5.9e6,
            5.8e6,
            5.8e6,
            5.7e6
        };

        //значения давления на выходе трубы
        pressure_out =
        {
            3e6,
            2.8e6,
            2.8e6,
            2.9e6,
            2.9e6,
            2.8e6,
            2.8e6,
            2.7e6
        };

        //значения расхода
        Q =
        {
            0.192325,
            0.200000,
            0.210000,
            0.200000,
            0.180000,
            0.210000,
            0.210000,
            0.210000
        };

        double dt = 60;

        for (size_t index = 0; index < Q.size(); index++) {
            // Заполняем временную сетку
            times.push_back(dt * index);
        }
    }
};

/// @brief  Проблемно-ориентированный слой для расчёта профиля массы с учётом движения партий на QUICKEST-ULTIMATE
struct quasi_layer_for_mass_calculation : public density_viscosity_quasi_layer<true>
{
    /// @brief Профиль массы в ячейках ЛУ
    vector<double> mass;
    /// @brief Инициализация профилей
    /// @param point_count Кол-во точек
    quasi_layer_for_mass_calculation(size_t point_count)
        : density_viscosity_quasi_layer<true>(point_count)
        , mass(point_count - 1)
    {}
};

/// @brief Задача расчёта движения партий на QUICKEST-ULTIMATE 
/// и расчёта профиля массы с учётом растяжения стенок трубы и сжимаемости нефти 
class mass_calculation_on_qsm_task_t {
public:
    /// @brief Тип буфера
    typedef ring_buffer_t<quasi_layer_for_mass_calculation> buffer_type;
private:
    // Модель трубы
    pipe_properties_t pipe;
    // Буфер
    buffer_type buffer;

public:
    /// @brief Конструктор
    /// @param pipe Модель трубопровода
    mass_calculation_on_qsm_task_t(const pipe_properties_t& pipe)
        : pipe(pipe)
        , buffer(2, pipe.profile.get_point_count())
    {
    }

    /// @brief Геттер для текущего слоя  
    quasi_layer_for_mass_calculation get_current_layer() {
        return buffer.current();
    }

    /// @brief Интерполяция давления в промежуточных точках 
    /// по начальному и конечному давлениям
    /// @param pressure_start Начальное давление
    /// @param pressure_end Конечное давление
    /// @return Профиль давления
    vector<double> interpolate_pressure_profile(double pressure_start, double pressure_end) {
        // Количество точек
        size_t n = pipe.profile.get_point_count();
        double min_coord = pipe.profile.coordinates.front();
        double length = pipe.profile.get_length();

        vector<double> pressure_profile;
        pressure_profile.reserve(n);

        for (size_t index = 0; index < n; index++) {
            double current_coordinate = pipe.profile.coordinates[index];
            double pressure = pressure_start + (pressure_end - pressure_start) * (current_coordinate - min_coord) / length;
            pressure_profile.push_back(pressure);
        }

        return pressure_profile;
    }

    /// @brief Начальный расчёт
    /// Считаем профиль плотности партии при нормальных условиях и профиль массы нефти в ячейках
    /// @param density_init Плотность начальной партии внутри трубопровода
    /// @param p_in Давление на входе
    /// @param p_out Давление на выходе
    void solve(double density_init, double p_in, double p_out)
    {
        // Количество точек
        size_t n = pipe.profile.get_point_count();

        auto& current = buffer.current();
        // Получаем профиль давления
        current.pressure = interpolate_pressure_profile(p_in, p_out);

        // Инициализация начального профиля плотности и массы
        for (size_t index = 0; index < n - 1; index++) {
            double cell_pressure = (current.pressure[index] + current.pressure[index + 1]) / 2;
            // Считаем умозрительную площадь сечения
            double A = calc_nominal_area(cell_pressure);
            current.density[index] = density_init / calc_density_ratio(cell_pressure);
            current.mass[index] = current.density[index] * A * (pipe.profile.coordinates[index + 1] - pipe.profile.coordinates[index]);
        }
    }
private:
    /// @brief Расчёт нового профиля плотности методом QUICKEST_ULTIMATE
    /// @param dt Временной шаг моделирования
    /// @param volumetric_flow Объёмный расход
    /// @param working_density Плотность партии на входе в трубопровод
    void make_rheology_step(double dt, double volumetric_flow, double density) {
        size_t n = pipe.profile.get_point_count();
        vector<double>Q_profile(n, volumetric_flow); // задаем по трубе новый расход из временного ряда

        advance(); // Сдвигаем текущий и предыдущий слои

        // считаем партии с помощью QUICKEST-ULTIMATE
        PipeQAdvection advection_model(pipe, Q_profile);

        // Шаг по плотности
        auto density_wrapper = buffer.get_buffer_wrapper(
            &density_viscosity_quasi_layer<1>::get_density_wrapper);

        quickest_ultimate_fv_solver solver_rho(advection_model, density_wrapper);
        solver_rho.step(dt, density, density);
       
    }

    /// @brief Расчёт отношения фактической площади сечения и пл-ди при нормальных условиях
    /// @param pressure Давление в ячейке
    /// @return Отношение фактической площади сечения и пл-ди при нормальных условиях
    double calc_area_ratio(double pressure) {
        double Bs = pipe.wall.getCompressionRatio();
        return 1 + Bs * (pressure - ATMOSPHERIC_PRESSURE);
    }

    /// @brief Расчёт отношения фактической плотности и плотности при нормальных условиях
    /// @param pressure Давление в ячейке
    /// @return Отношение фактической плотности и плотности при нормальных условиях
    double calc_density_ratio(double pressure) {
        return 1 + AVG_OIL_COMPRESSIBILITY_COEFF * (pressure - ATMOSPHERIC_PRESSURE);
    }
    
    /// @brief Расчёт площади сечения при нормальных условиях
    /// @param pressure Давление в ячейке
    /// @return Умозрительная площадь - пл-дь при нормальных условиях
    double calc_nominal_area(double pressure) {

        double S0 = pipe.wall.getArea();
        double A = calc_density_ratio(pressure) * S0 * calc_area_ratio(pressure);
        return A;
    }

public:
    /// @brief Расчёт шага моделирования, включающий в себя расчёт шага движения партии и массовый расчёт
    /// Функция делает сдвиг буфера (advance) так, что buffer.current после вызова содержит свежерасчитанный слой
    /// @param dt Временной шаг моделирования
    /// @param volumetric_flow Объёмный расход
    /// @param working_density Фактическая плотность на входе в линейный участок
    /// @param p_in Давление на входе
    /// @param p_out Давление на выходе
    void step(double dt, double volumetric_flow, double working_density, double p_in, double p_out) {
        auto& current = buffer.current();
        current.pressure = interpolate_pressure_profile(p_in, p_out);

        double nominal_density = working_density / calc_density_ratio(p_in);
        make_rheology_step(dt, volumetric_flow, nominal_density);

        for (size_t index = 0; index < pipe.profile.get_point_count() - 1; index++) {
            double cell_pressure = (current.pressure[index] + current.pressure[index + 1]) / 2;
            double A = calc_nominal_area(cell_pressure);
            double dx = pipe.profile.coordinates[index + 1] - pipe.profile.coordinates[index];
            current.mass[index] = current.density[index] * A * dx;
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

/// @brief Накапливает результаты по профилям массы
class mass_collector_from_calculation_on_qsm_t
    : public batch_processor_precalculated_times<quasi_layer_for_mass_calculation>
{
protected:
    /// @brief Вектор расчётных профилей массы ЛУ
    vector<vector<double>> mass_profile;
public:

    /// @brief Конструктор обработчика
    /// @param steps_count Количество шагов моделирования
    /// @param point_count Кол-во точек в профиле
    mass_collector_from_calculation_on_qsm_t(size_t steps_count, size_t point_count)
        : mass_profile(steps_count, vector<double>(point_count))
    {

    }
    virtual void process_data(size_t step_index,
        const quasi_layer_for_mass_calculation& layer) override
    {
        // Копирование данных
        std::copy(layer.mass.begin(), layer.mass.end(), mass_profile[step_index].begin());
    }

    /// @brief Геттер для вектора собранных результатов расчёта профилей массы ЛУ
    /// @return Вектор расчётных профилей массы ЛУ
    const vector<vector<double>>& get_mass_profile_calculated() const {
        return mass_profile;
    }
};

/// @brief Массовый расчёт нефти внутри трубопровода
/// с учётом растяжения стенок и сжимаемости самой нефти,
/// а также с учётом движения партий на QUICKEST-ULTIMATE
/// @param task Задача расчёта движения партий и расчёта профиля массы
/// @param times Предпосчитанная временная сетка моделирования работы ЛУ
/// @param working_density Краевые условия по плотности
/// @param p_in Краевые условия по давлению в начале ЛУ
/// @param p_out Краевые условия по давлению в конце ЛУ
/// @param volumetric_flow Краевые условия по объёмному расходу
/// @param data_processor Обработчик результатов расчета
void perform_mass_calculation_on_qsm(
    mass_calculation_on_qsm_task_t& task,
    const vector<double>& times,
    const vector<double>& density,
    const vector<double>& p_in,
    const vector<double>& p_out,
    const vector<double>& volumetric_flow,
    mass_collector_from_calculation_on_qsm_t* data_processor
) {
    // Проводим начальный расчёт
    task.solve(density.front(), p_in.front(), p_out.front());

    data_processor->process_data(0, task.get_buffer().current());

    for (size_t step_index = 1; step_index < times.size(); step_index++)
    {
        double time_step = times[step_index] - times[step_index - 1];

        task.step(time_step, volumetric_flow[step_index], density[step_index], p_in[step_index], p_out[step_index]);

        data_processor->process_data(step_index, task.get_buffer().current());
    }
}


/// @brief Расчёт профилей массы
TEST_F(CalcMassQSM, CalcMass)
{
    mass_calculation_on_qsm_task_t task(pipe);
    mass_collector_from_calculation_on_qsm_t collector(times.size(), pipe.profile.get_point_count());

    perform_mass_calculation_on_qsm(
        task,
        times,
        density,
        pressure_in,
        pressure_out,
        Q,
        &collector
    );

    vector<vector<double>> mass_profiles = collector.get_mass_profile_calculated();

}


