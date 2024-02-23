#include <random>
#include <algorithm>
const double g = 9.81, pi = M_PI;
using namespace std;

struct my_pipe_parameters
{
    /// @brief Длина трубы, (м)
    double length = 80e3;
    /// @brief  Внешний диаметр трубы, (м)
    double external_diameter = 0.72;
    /// @brief Толщина стенки трубы, (м)
    double delta_d = 0.01;
    /// @brief Абсолютная шероховатость, (м)
    double abs_roughness = 15e-6;
    /// @brief Внутренний диаметр трубы, (м)
    double internal_diameter = external_diameter - 2 * delta_d;
    /// @brief Шероховатость
    double roughness = abs_roughness / internal_diameter;
    /// @brief Начальная высотная отметка, (м)
    double z_0 = 100;
    /// @brief Конечная высотная отметка, (м)
    double z_L = 50;
    /// @brief Шаг сетки, (м)
    double h = 1e3;
    /// @brief Количество шагов
    size_t n;
    my_pipe_parameters(double length, double external_diameter, double delta_d, double abs_roughness, double z_0, double z_L, double h) :
        length{ length }, external_diameter{ external_diameter }, delta_d{ delta_d }, abs_roughness{ abs_roughness }, z_0{ z_0 }, z_L{ z_L }, h{ h }
    {
        n = static_cast<int>(length / h + 0.5) + 1;
    }
};

struct my_task_parameters
{
    my_pipe_parameters& pipe;
    /// @brief Плотность жидкости, (кг/м^3)
    double rho = 870;
    /// @brief Кинематическая вязкость, (м^2/с)
    double nu = 15e-6;
    /// @brief Давление в начале участка, (Па)
    double p_0 = 5e6;
    /// @brief Давление в конце участка, (Па)
    double p_L = 0.6e6;
    /// @brief Расход жидкости, (м^3/с)
    double Q = 0.972;

    vector<double> rho_profile = vector<double>(pipe.n, rho);
    vector<double> nu_profile = vector<double>(pipe.n, nu);
    my_task_parameters(my_pipe_parameters& pipe, double rho, double nu, double p_0, double p_L, double Q, vector<double> rho_profile, vector<double> nu_profile) :
        pipe{ pipe }, rho{ rho }, nu{ nu }, p_0{ p_0 }, p_L{ p_L }, Q{ Q }, rho_profile{ rho_profile }, nu_profile{ nu_profile }
    {
    }
};

struct print_data {
    vector<double> time;           // время с
    vector<double> density;        // вытесняющая плотность кг/м3
    vector<double> viscosity;      // вытесняющая вязкость Ст
    vector<double> inputPressure;  // давление на входе Па
    vector<double> outputPressure; // давление на выходе Па
    vector<double> flowRate;       // расход м3/с
};

/// @brief Функция расчета скорости из расхода
/// @param Q Расход, (м^3/с)
/// @param internal_diameter Внутренний диаметр, (м)
/// @return Скорость, (м/с)
double calc_speed(double Q, double internal_diameter)
{
    double speed = (4 * Q) / (pow(internal_diameter, 2) * pi);
    return speed;
}

/// @brief Функция расчета расхода из скорости
/// @param speed Скорость, (м/с)
/// @param internal_diameter Внутренний диаметр, (м) 
/// @return Расход, (м^3/с)
double calc_flow(double speed, double internal_diameter)
{
    double flow = (pow(internal_diameter, 2) * pi * speed) / 4;
    return flow;
}

/// @brief Функция расчета числа Рейнольдса
/// @param speed Скорость, (м/с)
/// @param internal_diameter Внутренний диаметр, (м)
/// @param nu Кинетическая вязкость, (м^2/с)
/// @return Число Рейнольдса
double calc_Re(double speed, double internal_diameter, double nu)
{
    double Re = (speed * internal_diameter) / nu;
    return Re;
}

/// @brief Функция расчета касательного напряжения трения
/// @param hydraulic_resistance гидравлическое сопротивление
/// @param rho Плотность, (кг/см^2)
/// @param speed Скорость, (м/с)
/// @return касательное напряжение трения
double calc_tau(double hydraulic_resistance, double rho, double speed)
{
    double tau = (hydraulic_resistance / 8) * rho * pow(speed, 2);
    return tau;
}

/// @brief Класс, решающий задачи PQ, QP методом Эйлера
class euler_solver_with_MOC
{
    const my_pipe_parameters& pipe;
    const my_task_parameters& task;
public:
    euler_solver_with_MOC(const my_pipe_parameters& pipe, const my_task_parameters& task) :
        pipe(pipe), task(task)
    {
    }
    /// @brief Метод нахождения входного давления
    /// @return Pвх
    vector<double> euler_solver_PQ()
    {
        double delta_z = (pipe.z_L - pipe.z_0) / (pipe.n - 1);
        double speed = calc_speed(task.Q, pipe.internal_diameter);
        vector<double> p_profile = vector<double>(pipe.n);
        p_profile[0] = task.p_0;
        for (int i = 1; i < pipe.n; i++)
        {
            double Re = calc_Re(speed, pipe.internal_diameter, task.nu_profile[i]);
            double hydraulic_resistance = hydraulic_resistance_altshul(Re, pipe.roughness);
            double tau = calc_tau(hydraulic_resistance, task.rho_profile[i], speed);
            p_profile[i] = p_profile[i - 1] + pipe.h * ((-4 / pipe.internal_diameter) * tau - task.rho_profile[i] * g * (delta_z / pipe.h));

        }
        return p_profile;
    }
    /// @brief Метод нахождения выходного давления
    /// @return Pвых
    vector<double> euler_solver_QP()
    {
        double delta_z = (pipe.z_L - pipe.z_0) / (pipe.n - 1);
        double speed = calc_speed(task.Q, pipe.internal_diameter);
        vector<double> p_profile = vector<double>(pipe.n);
        p_profile[pipe.n - 1] = task.p_L;
        for (int i = pipe.n - 2; i >= 0; i--)
        {
            double Re = calc_Re(speed, pipe.internal_diameter, task.nu_profile[i]);
            double hydraulic_resistance = hydraulic_resistance_altshul(Re, pipe.roughness);
            double tau = calc_tau(hydraulic_resistance, task.rho_profile[i], speed);
            p_profile[i] = p_profile[i + 1] - pipe.h * ((-4 / pipe.internal_diameter) * tau - task.rho_profile[i] * g * (delta_z / pipe.h));

        }
        return p_profile;
    }
};

class newton_solver_PP_with_euler_with_MOC : public fixed_system_t<1>
{
    const my_pipe_parameters& pipe;
    const my_task_parameters& task;
    using fixed_system_t<1>::var_type;
public:
    newton_solver_PP_with_euler_with_MOC(const my_pipe_parameters& pipe, const my_task_parameters& task) :
        pipe(pipe), task(task)
    {
    }
    /// @brief Задание функции невязок
    /// @param x Искомый расход
    /// @return Функция невязок
    var_type residuals(const var_type& x)
    {
        my_task_parameters temp_task = task; // Временная структура
        temp_task.Q = x; // во временной структуре используем Q для нашего уравнения невязки, эта Q будет идти в солвер
        euler_solver_with_MOC e_solver(pipe, temp_task); // Объявляем переменную класса солвера Эйлером

        p_profile = e_solver.euler_solver_PQ(); // Считаем профиль давлений Эйлером
        return (p_profile.back() - task.p_L);
    }
    vector<double> get_p_profile() const {
        return p_profile;
    }
private:
    vector<double> p_profile;
};

double linear_interpolator(vector<double> original_time, vector<double> original_value, double new_time_step)
{
    size_t index1 = 0;
    size_t index2 = 1;
    while (original_time[index2] < new_time_step)
    {
        ++index1;
        ++index2;
    }
    double t1 = original_time[index1];
    double t2 = original_time[index2];
    double value1 = original_value[index1];
    double value2 = original_value[index2];
    return value1 + (value2 - value1) * (new_time_step - t1) / (t2 - t1);
}

void print_data_to_csv(const vector<double>& time_rho,
    const vector<double>& rho_and_nu_in_0,
    const vector<double>& time_nu,
    const vector<double>& rho_and_nu_in_1,
    const vector<double>& time_p_in,
    const vector<double>& p_in,
    const vector<double>& time_p_out,
    const vector<double>& p_out,
    const vector<double>& time_Q,
    const vector<double>& Q,
    const wstring& filename)
{

    // Определяем максимальную длину вектора
    size_t maxLength = max({ time_rho.size(), rho_and_nu_in_0.size(), time_nu.size(), rho_and_nu_in_1.size(),
                                      time_p_in.size(), p_in.size(), time_p_out.size(), p_out.size(), Q.size(), time_Q.size() });

    // Открываем файл для записи
    ofstream file(filename);

    // Проверяем, открыт ли файл успешно
    if (file.is_open()) {
        // Записываем заголовки столбцов
        file << "Time Density; Density; Time Viscosity; Viscosity; Time Pressure In; Pressure In; Time Pressure Out; Pressure Out; Time Flow Rate; Flow Rate\n";

        // Записываем данные из векторов
        for (size_t i = 0; i < maxLength; ++i) {
            // Если индекс находится в пределах длины вектора, записываем значение, иначе записываем пустую ячейку
            file << (i < time_rho.size() ? to_string(time_rho[i]) : "") << ";"
                << (i < rho_and_nu_in_0.size() ? to_string(rho_and_nu_in_0[i]) : "") << ";"
                << (i < time_nu.size() ? to_string(time_nu[i]) : "") << ";"
                << (i < rho_and_nu_in_1.size() ? to_string(rho_and_nu_in_1[i]) : "") << ";"
                << (i < time_p_in.size() ? to_string(time_p_in[i]) : "") << ";"
                << (i < p_in.size() ? to_string(p_in[i]) : "") << ";"
                << (i < time_p_out.size() ? to_string(time_p_out[i]) : "") << ";"
                << (i < p_out.size() ? to_string(p_out[i]) : "") << ";"
                << (i < time_Q.size() ? to_string(time_Q[i]) : "") << ";"
                << (i < Q.size() ? to_string(Q[i]) : "") << "\n";
        }
        // Закрываем файл
        file.close();
    }
}


void print_layers(const double dt,
    const vector<double>& layer,
    const wstring& filename)
{
    ofstream  file(filename, ios::app);
    if (dt == 0)
    {
        // Если файл существует, очищаем его содержимое
        if (file.is_open()) {
            file.close();
            file.open(filename, ios::out | ios::trunc);
        }
        // Если файл не существует, создаем новый
        else {
            file.open(filename, ios::out);
        }
        // Файл существует, но file.is_open() возвращает False
        if (!file.is_open()) {
            file.clear();  // Очищаем флаг ошибки
            file.open(filename, ios::out | ios::trunc);
        }
    }
    if (file.is_open()) {
        file << to_string(dt) << ";";
        for (int j = 0; j < layer.size(); j++)
        {
            file << to_string(layer[j]) << ";";
        }
        file << "\n";
        file.close();
    }
}


/// @brief Тесты для солвера quickest_ultimate_fv_solver
class QuickWithQuasiStationaryModel : public ::testing::Test {
protected:
    // Профиль переменных
    typedef quickest_ultimate_fv_solver_traits<1>::var_layer_data target_var_t;
    typedef quickest_ultimate_fv_solver_traits<1>::specific_layer specific_data_t;

    // Слой: переменных Vars + сколько угодно служебных Specific
    typedef composite_layer_t<target_var_t, specific_data_t> layer_t;
protected:
    /// @brief Параметры трубы
    pipe_properties_t pipe;
    /// @brief Профиль расхода
    vector<double> Q;
    std::unique_ptr<PipeQAdvection> advection_model;
    std::unique_ptr<ring_buffer_t<layer_t>> buffer;
protected:

    /// @brief Подготовка к расчету для семейства тестов
    virtual void SetUp() override {
        // Упрощенное задание трубы - 50км, с шагом разбиения для расчтной сетки 1км, диаметром 700мм
        simple_pipe_properties simple_pipe;
        //simple_pipe.length = 50e3;
        simple_pipe.length = 700e3; // тест трубы 700км
        //simple_pipe.diameter = 0.7;
        simple_pipe.diameter = 0.514; // тест трубы 700км
        //simple_pipe.dx = 100;
        simple_pipe.dx = 100; // тест трубы 700км
        pipe = pipe_properties_t::build_simple_pipe(simple_pipe);

        Q = vector<double>(pipe.profile.getPointCount(), 0.5);
        advection_model = std::make_unique<PipeQAdvection>(pipe, Q);
        buffer = std::make_unique<ring_buffer_t<layer_t>>(2, pipe.profile.getPointCount());

        layer_t& prev = buffer->previous();
        prev.vars.cell_double[0] = vector<double>(prev.vars.cell_double[0].size(), 850);
    }
};

/// @brief Пример вывода в файл через
TEST_F(QuickWithQuasiStationaryModel, TimeRowsTask)
{
    string path = prepare_test_folder();
    // Упрощенное задание трубы - 50км, с шагом разбиения для расчтной сетки 1км, диаметром 700мм
    simple_pipe_properties simple_pipe;
    //simple_pipe.length = 50e3;
    simple_pipe.length = 700e3; // тест трубы 700км
    //simple_pipe.diameter = 0.7;
    simple_pipe.diameter = 0.514; // тест трубы 700км
    const auto& x = advection_model->get_grid();
    double dx = x[1] - x[0];
    double v = advection_model->getEquationsCoeffs(0, 0);
    double delta_d = 0.01, abs_roughness = 15e-6, z_0 = 50, z_L = 100,
        rho = 850, nu = 15e-6, p_0 = 6e6, p_L = 0.0;
    my_pipe_parameters my_pipe{ simple_pipe.length, simple_pipe.diameter, delta_d, abs_roughness, z_0, z_L, dx };
    double Q = calc_flow(v, my_pipe.internal_diameter);
    my_task_parameters my_task{ my_pipe, rho, nu, p_0, p_L, Q, {}, {} };

    std::srand(std::time(nullptr));
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dis(7.0, 10.0);

    double rho_in = 860;
    double rho_out = 870;
    double T = 350000; // период моделирования

    vector<double> time_row_for_rho;
    for (double time = 0.0; time <= T; time += dis(gen))
        time_row_for_rho.push_back(time);
    vector<double> time_rho_in_row = vector<double>(time_row_for_rho.size(), rho);
    for (size_t i = time_row_for_rho.size()/2; i <= time_row_for_rho.size() - 1; i++)
    {
        time_rho_in_row[i] = rho_in;//rho - abs(rho * 0.001 * std::rand() / RAND_MAX);
    }

    vector<double> time_row_for_nu;
    for (double time = 0.0; time <= T; time += dis(gen))
        time_row_for_nu.push_back(time);
    vector<double> time_nu_in_row = vector<double>(time_row_for_nu.size(), nu);
    for (size_t i = 1; i <= time_row_for_nu.size() - 1; i++)
    {
        time_nu_in_row[i] = nu;//nu - abs(nu * 0.1 * std::rand() / RAND_MAX);
    }

    vector<double> time_rho_out_row = vector<double>(time_row_for_rho.size(), rho);
    vector<double> time_nu_out_row = vector<double>(time_row_for_nu.size(), nu);

    vector<double> time_row_for_p_in;
    for (double time = 0.0; time <= T; time += dis(gen))
        time_row_for_p_in.push_back(time);
    vector<double> time_p_in_row = vector<double>(time_row_for_p_in.size(), p_0);
    for (size_t i = 1; i <= time_row_for_p_in.size() - 1; i++)
    {
        time_p_in_row[i] = p_0;//p_0 - abs(p_0 * 2e-4 * std::rand() / RAND_MAX);
    }

    vector<double> time_row_for_Q;
    for (double time = 0.0; time <= T; time += dis(gen))
        time_row_for_Q.push_back(time);
    vector<double> time_Q_row = vector<double>(time_row_for_Q.size(), v);
    for (size_t i = 1; i <= time_row_for_Q.size() - 1; i++)
    {
        time_Q_row[i] = v - abs(v * 0.1 * std::rand() / RAND_MAX);
    }
    
    double Q_max = *std::max_element(time_Q_row.begin(), time_Q_row.end());
    double v_max = calc_speed(Q_max, simple_pipe.diameter);

    double dt_ideal = abs(dx / v_max);
    double Cr = v_max * dt_ideal / dx; // равен 1, взяли максимальную скорость


    ///////////////////////////////////////////////////////////////////////////////////
    

    advection_model = std::make_unique<PipeQAdvection>(pipe, Q);
    buffer = std::make_unique<ring_buffer_t<layer_t>>(2, pipe.profile.getPointCount());

    layer_t& prev = buffer->previous();
    prev.vars.cell_double[0] = vector<double>(prev.vars.cell_double[0].size(), 850);

    vector<double> new_time_row, new_time_p_in_row, new_time_p_out_row, new_time_Q_row, p_profile, initial_p_profile;

    double t = 0; // текущее время
    double dt = Cr * dt_ideal; // время в долях от Куранта

    std::stringstream filename;
    filename << path << "output Cr=" << Cr << ".csv";
    std::ofstream output(filename.str());

    do {
        new_time_row.push_back(t);
        Cr = time_v_row[0] * dt_ideal / dx;
        if (t == 0) {
            layer_t& prev = buffer->previous();
            prev.vars.print(t, output);
        }

        t += dt;

        quickest_ultimate_fv_solver solver(*advection_model, *buffer);
        solver.step(dt, rho_in, rho_out);

        layer_t& next = buffer->current();
        next.vars.print(t, output);


        buffer->advance(+1);
            

        

        //все готово для следующего шага, интерполированная скорость есть
    } while (t < T);
    output.flush();
    output.close();
    ///////////////////////////////////////////////////////////////////////////////////
}