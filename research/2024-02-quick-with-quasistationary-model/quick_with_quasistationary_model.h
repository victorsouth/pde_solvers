#include <random>
#include <algorithm>


const double g = 9.81, pi = M_PI;
using namespace std;



struct my_task_parameters
{
    /// @brief Плотность жидкости, (кг/м^3)
    double rho;
    /// @brief Кинематическая вязкость, (сСТ)
    double nu;
    /// @brief Давление в начале участка, (Па)
    double p_0;
    /// @brief Давление в конце участка, (Па)
    double p_L;
    /// @brief Расход жидкости, (м^3/с)
    double Q;
    /// @brief Профиль плотностей по всей трубе, (кг/м^3)
    vector<double>& rho_profile;
    /// @brief Профиль вязкостей по всей трубе, (сСТ)
    vector<double>& nu_profile;
    my_task_parameters(double rho, double nu, double p_0, double p_L, double Q, vector<double> rho_profile, vector<double> nu_profile) :
        rho{ rho }, nu{ nu }, p_0{ p_0 }, p_L{ p_L }, Q{ Q }, rho_profile{ rho_profile }, nu_profile{ nu_profile }
    {
    }
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
    const string& filename)
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
    const string& filename)
{
    ofstream  file(filename, ios::app);
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

void clear_directory(const filesystem::path& dir_path) {
    for (const auto& entry : filesystem::directory_iterator(dir_path)) {
        filesystem::remove_all(entry.path());
    }
}


/// @brief Тесты для солвера quickest_ultimate_fv_solver
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
    vector<double> Q_profile;

    
    std::unique_ptr<PipeQAdvection> advection_model;
    std::unique_ptr<ring_buffer_t<layer_t>> buffer, buffer_nu;
protected:

    /// @brief Подготовка к расчету для семейства тестов
    virtual void SetUp() override {

        simple_pipe_properties simple_pipe;
        simple_pipe.length = 200e3; // тест трубы 700км
        simple_pipe.diameter = 0.514; // тест трубы 700км
        //pipe.wall.diameter = 0.514;
        simple_pipe.dx = 100; // тест трубы 700км
        
        pipe = pipe_properties_t::build_simple_pipe(simple_pipe);

        std::srand(std::time(nullptr));
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<double> dis(57.0, 67.0);


        Q_profile = vector<double>(pipe.profile.getPointCount(), 0.2);


        advection_model = std::make_unique<PipeQAdvection>(pipe, Q_profile);
        buffer = std::make_unique<ring_buffer_t<layer_t>>(2, pipe.profile.getPointCount());
        buffer_nu = std::make_unique<ring_buffer_t<layer_t>>(2, pipe.profile.getPointCount());
    }
};

/// @brief Уравнение трубы для задачи PQ
class Pipe_model_for_PQ_t : public ode_t<1>
{
public:
    using ode_t<1>::equation_coeffs_type;
    using ode_t<1>::right_party_type;
    using ode_t<1>::var_type;
protected:
    vector<double>& rho_profile;
    vector<double>& nu_profile;
    pipe_properties_t& pipe;
    double flow;
    mutable size_t flag = 0;
public:
    /// @brief Констуктор уравнения трубы
    /// @param pipe Ссылка на сущность трубы
    /// @param oil Ссылка на сущность нефти
    /// @param flow Объемный расход
    Pipe_model_for_PQ_t(pipe_properties_t& pipe, vector<double>& rho_profile, vector<double>& nu_profile, double flow)
        : pipe(pipe)
        , rho_profile(rho_profile)
        , nu_profile(nu_profile)
        , flow(flow)
    {

    }

    /// @brief Возвращает известную уравнению сетку
    virtual const vector<double>& get_grid() const override {
        return pipe.profile.coordinates;
    }

    /// @brief Возвращает значение правой части ДУ
    /// @param grid_index Обсчитываемый индекс расчетной сетки
    /// @param point_vector Начальные условия
    /// @return Значение правой части ДУ в точке point_vector
    virtual right_party_type ode_right_party(
        size_t grid_index, const var_type& point_vector) const override
    {
        double rho, S_0, v, Re, lambda, tau_w, height_derivative;
        
        size_t reo_index;
        if (flag % 2)
            reo_index = grid_index - 1;
        else
            reo_index = grid_index;
        flag++;
        rho = rho_profile[reo_index];
        S_0 = pipe.wall.getArea();
        v = flow / (S_0);
        Re = v * pipe.wall.diameter / nu_profile[reo_index];
        lambda = pipe.resistance_function(Re, pipe.wall.relativeRoughness());
        tau_w = lambda / 8 * rho * v * abs(v);
        /// Обработка индекса в случае расчетов на границах трубы
        /// Чтобы не выйти за массив высот, будем считать dz/dx в соседней точке
        grid_index = grid_index == 0 ? grid_index + 1 : grid_index;
        grid_index = grid_index == pipe.profile.heights.size() - 1 ? grid_index - 1 : grid_index;

        height_derivative = (pipe.profile.heights[grid_index] - pipe.profile.heights[grid_index - 1]) /
            (pipe.profile.coordinates[grid_index] - pipe.profile.coordinates[grid_index - 1]);
        return { ((-4) / pipe.wall.diameter) * tau_w - rho * M_G * height_derivative };  
    }
};

class newton_solver_PP_with_euler : public fixed_system_t<1>
{
    const my_task_parameters& task;
    using fixed_system_t<1>::var_type;
protected:
    pipe_properties_t& pipe;
    double Q_approx;
public:
    newton_solver_PP_with_euler(pipe_properties_t& pipe, const my_task_parameters& task, double Q_approx)
        : pipe(pipe)
        , task(task)
        , Q_approx(Q_approx)
    {

    }
    /// @brief Задание функции невязок
    /// @param x Искомый расход
    /// @return Функция невязок
    var_type residuals(const var_type& x)
    {
        Q_approx = x; // во временной структуре используем Q для нашего уравнения невязки, эта Q будет идти в солвер
        // Объявляем переменную класса солвера Эйлером
        Pipe_model_for_PQ_t pipeModel(pipe, task.rho_profile, task.nu_profile, Q_approx);

        solve_euler_corrector<1>(pipeModel, 1, task.p_0, &p_profile);
        p_profile.pop_back();
        return (p_profile.back() - task.p_L);
    }
    vector<double> get_p_profile() const {
        return p_profile;
    }
private:
    vector<double> p_profile;
};

/// @brief Пример вывода в файл через
TEST_F(QuickWithQuasiStationaryModel, TimeRowsTask)
{

    const auto& x = advection_model->get_grid();
    double dx = x[1] - x[0];
    double v = advection_model->getEquationsCoeffs(0, 0);

    double rho = 800, nu = 15e-6, p_0 = 6e6;

    // Задаем объемнй расход нефти, [м3/с]
    double Q = calc_flow(v,pipe.wall.diameter);
    double T = 450000; // период моделирования

    std::srand(std::time(nullptr));
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dis(57.0, 67.0);

    double rho_in = 900;
    double rho_out = 870;

    vector<double> time_row_for_rho;
    for (double time = 0.0; time <= T; time += dis(gen))
        time_row_for_rho.push_back(time);
    vector<double> time_rho_in_row = vector<double>(time_row_for_rho.size(), rho);
    for (size_t i = 0; i <= time_row_for_rho.size() - 1; i++)
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

    vector<double> time_Q_row = vector<double>(time_row_for_Q.size(), p_0);

    for (size_t i = 0; i <= 14; i++)
    {
        time_Q_row[i] = 0.2 - abs(0.2 * 0.001 * std::rand() / RAND_MAX);
    }
    for (size_t i = 15; i <= time_Q_row.size() - 1; i++)
    {
        time_Q_row[i] = 0.1 - abs(0.1 * 0.001 * std::rand() / RAND_MAX);
    }

    double v_max = calc_speed(*std::max_element(time_Q_row.begin(), time_Q_row.end()), pipe.wall.diameter);

    double dt = abs(dx / v_max);// Cr равен 1, взяли максимальную скорость

    // Вектор для хранения профиля давления
    layer_t& prev = buffer->previous();
    prev.vars.cell_double[0] = vector<double>(prev.vars.cell_double[0].size(), rho);

    layer_t& prev_nu = buffer_nu->previous();
    prev_nu.vars.cell_double[0] = vector<double>(prev.vars.cell_double[0].size(), nu);

    vector<double> new_time_row, new_time_p_in_row, new_time_p_out_row, new_time_Q_row, initial_p_profile, rho_profile, nu_profile;
    vector<double> diff_p_profile = vector<double>(pipe.profile.getPointCount(), 0);
    vector<vector<double>> rho_and_nu_in = vector<vector<double>>(2);
    vector<vector<double>> rho_and_nu_out = vector<vector<double>>(2);

    double t = 0; // текущее время


    string path_for_quick = string("../research_out/QuickWithQSM/");
    string p_profile_file = string("../research_out/QuickWithQSM/p_profile.csv");
    string rho_profile_file = string("../research_out/QuickWithQSM/rho_profile.csv");
    string nu_profile_file = string("../research_out/QuickWithQSM/nu_profile.csv");
    string diff_p_profile_file = string("../research_out/QuickWithQSM/diff_p_profile.csv");

    filesystem::create_directories(path_for_quick);
    clear_directory(path_for_quick);
    do {
        vector<double> p_profile(pipe.profile.getPointCount());
        new_time_row.push_back(t);
        Q = linear_interpolator(time_row_for_Q, time_Q_row, t);
        new_time_Q_row.push_back(Q);
        rho_and_nu_in[0].push_back(linear_interpolator(time_row_for_rho, time_rho_in_row, t));
        rho_and_nu_in[1].push_back(linear_interpolator(time_row_for_nu, time_nu_in_row, t));

        rho_and_nu_out[0].push_back(linear_interpolator(time_row_for_rho, time_rho_out_row, t));
        rho_and_nu_out[1].push_back(linear_interpolator(time_row_for_nu, time_nu_out_row, t));

        new_time_p_in_row.push_back(linear_interpolator(time_row_for_p_in, time_p_in_row, t));

        double Cr = calc_speed(Q, pipe.wall.diameter) * dt / dx;


        if (t == 0) {
            Pipe_model_for_PQ_t pipeModel(pipe, prev.vars.cell_double[0], prev_nu.vars.cell_double[0], Q);

            solve_euler_corrector<1>(pipeModel, 1, p_0, &p_profile);
            
            initial_p_profile = p_profile;
            print_layers(t, p_profile, p_profile_file);
            print_layers(t, prev.vars.cell_double[0], rho_profile_file);
            print_layers(t, prev_nu.vars.cell_double[0], nu_profile_file);
            print_layers(t, diff_p_profile, diff_p_profile_file);
            new_time_p_out_row.push_back(p_profile.back());
            p_profile = vector<double>(pipe.profile.getPointCount());
        }
        t += dt;
        Q_profile = vector<double>(time_row_for_Q.size(), new_time_Q_row.back());
        advection_model = std::make_unique<PipeQAdvection>(pipe, Q_profile);

        quickest_ultimate_fv_solver solver_rho(*advection_model, *buffer);
        solver_rho.step(dt, rho_and_nu_in[0].back(), rho_and_nu_out[0].back());
        layer_t& next = buffer->current();
        rho_profile = next.vars.cell_double[0];

        quickest_ultimate_fv_solver solver_nu(*advection_model, *buffer_nu);
        solver_nu.step(dt, rho_and_nu_in[1].back(), rho_and_nu_out[1].back());
        layer_t& next_nu = buffer_nu->current();
        nu_profile = next_nu.vars.cell_double[0];

        Pipe_model_for_PQ_t pipeModel(pipe, rho_profile, nu_profile, Q);

        solve_euler_corrector<1>(pipeModel, 1, new_time_p_in_row.back(), &p_profile);

        new_time_p_out_row.push_back(p_profile.back());
        std::transform(initial_p_profile.begin(), initial_p_profile.end(), p_profile.begin(), diff_p_profile.begin(),
            [](double initial, double current) {return initial - current;  });
        print_layers(t, p_profile, p_profile_file);
        print_layers(t, rho_profile, rho_profile_file);
        print_layers(t, nu_profile, nu_profile_file);
        print_layers(t, diff_p_profile, diff_p_profile_file);

        buffer->advance(+1);
        buffer_nu->advance(+1);
    } while (t < T - dt);
    string filename_initial = string("../research_out/QuickWithQSM/initial_data.csv");
    print_data_to_csv(
        time_row_for_rho,
        time_rho_in_row,
        time_row_for_nu,
        time_nu_in_row,
        time_row_for_p_in,
        time_p_in_row,
        {},
        {},
        time_row_for_Q,
        Q_profile,
        filename_initial
    );
    string filename_final = string("../research_out/QuickWithQSM/final_data.csv");
    print_data_to_csv(
        new_time_row,
        rho_and_nu_in[0],
        new_time_row,
        rho_and_nu_in[1],
        new_time_row,
        new_time_p_in_row,
        new_time_row,
        new_time_p_out_row,
        new_time_row,
        new_time_Q_row,
        filename_final
    );
}