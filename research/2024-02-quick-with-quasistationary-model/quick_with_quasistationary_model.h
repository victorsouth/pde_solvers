#include <random>
#include <algorithm>


const double g = 9.81, pi = M_PI;
using namespace std;

struct my_pipe_parameters
{
    /// @brief ����� �����, (�)
    double length = 80e3;
    /// @brief  ������� ������� �����, (�)
    double external_diameter = 0.72;
    /// @brief ������� ������ �����, (�)
    double delta_d = 0.01;
    /// @brief ���������� �������������, (�)
    double abs_roughness = 15e-6;
    /// @brief ���������� ������� �����, (�)
    double internal_diameter = external_diameter - 2 * delta_d;
    /// @brief �������������
    double roughness = abs_roughness / internal_diameter;
    /// @brief ��������� �������� �������, (�)
    double z_0 = 100;
    /// @brief �������� �������� �������, (�)
    double z_L = 50;
    /// @brief ��� �����, (�)
    double h = 1e3;
    /// @brief ���������� �����
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
    /// @brief ��������� ��������, (��/�^3)
    double rho = 870;
    /// @brief �������������� ��������, (�^2/�)
    double nu = 15e-6;
    /// @brief �������� � ������ �������, (��)
    double p_0 = 5e6;
    /// @brief �������� � ����� �������, (��)
    double p_L = 0.6e6;
    /// @brief ������ ��������, (�^3/�)
    double Q = 0.972;

    vector<double> rho_profile = vector<double>(pipe.n, rho);
    vector<double> nu_profile = vector<double>(pipe.n, nu);
    my_task_parameters(my_pipe_parameters& pipe, double rho, double nu, double p_0, double p_L, double Q, vector<double> rho_profile, vector<double> nu_profile) :
        pipe{ pipe }, rho{ rho }, nu{ nu }, p_0{ p_0 }, p_L{ p_L }, Q{ Q }, rho_profile{ rho_profile }, nu_profile{ nu_profile }
    {
    }
};

struct print_data {
    vector<double> time;           // ����� �
    vector<double> density;        // ����������� ��������� ��/�3
    vector<double> viscosity;      // ����������� �������� ��
    vector<double> inputPressure;  // �������� �� ����� ��
    vector<double> outputPressure; // �������� �� ������ ��
    vector<double> flowRate;       // ������ �3/�
};

/// @brief ������� ������� �������� �� �������
/// @param Q ������, (�^3/�)
/// @param internal_diameter ���������� �������, (�)
/// @return ��������, (�/�)
double calc_speed(double Q, double internal_diameter)
{
    double speed = (4 * Q) / (pow(internal_diameter, 2) * pi);
    return speed;
}

/// @brief ������� ������� ������� �� ��������
/// @param speed ��������, (�/�)
/// @param internal_diameter ���������� �������, (�) 
/// @return ������, (�^3/�)
double calc_flow(double speed, double internal_diameter)
{
    double flow = (pow(internal_diameter, 2) * pi * speed) / 4;
    return flow;
}

/// @brief ������� ������� ����� ����������
/// @param speed ��������, (�/�)
/// @param internal_diameter ���������� �������, (�)
/// @param nu ������������ ��������, (�^2/�)
/// @return ����� ����������
double calc_Re(double speed, double internal_diameter, double nu)
{
    double Re = (speed * internal_diameter) / nu;
    return Re;
}

/// @brief ������� ������� ������������ ���������� ������
/// @param hydraulic_resistance �������������� �������������
/// @param rho ���������, (��/��^2)
/// @param speed ��������, (�/�)
/// @return ����������� ���������� ������
double calc_tau(double hydraulic_resistance, double rho, double speed)
{
    double tau = (hydraulic_resistance / 8) * rho * pow(speed, 2);
    return tau;
}

/// @brief �����, �������� ������ PQ, QP ������� ������
class euler_solver
{
    const my_pipe_parameters& my_pipe;
    const my_task_parameters& my_task;
public:
    euler_solver(const my_pipe_parameters& my_pipe, const my_task_parameters& my_task) :
        my_pipe(my_pipe), my_task(my_task)
    {
    }
    /// @brief ����� ���������� �������� ��������
    /// @return P��
    vector<double> euler_solver_PQ()
    {
        double delta_z = (my_pipe.z_L - my_pipe.z_0) / (my_pipe.n - 1);
        double speed = calc_speed(my_task.Q, my_pipe.internal_diameter);
        vector<double> p_profile = vector<double>(my_pipe.n);
        p_profile[0] = my_task.p_0;
        for (int i = 1; i < my_pipe.n - 1; i++)
        {
            double Re = calc_Re(speed, my_pipe.internal_diameter, my_task.nu_profile[i]);
            double hydraulic_resistance = hydraulic_resistance_isaev(Re, my_pipe.roughness);
            double tau = calc_tau(hydraulic_resistance, my_task.rho_profile[i], speed);
            p_profile[i] = p_profile[i - 1] + my_pipe.h * ((-4 / my_pipe.internal_diameter) * tau - my_task.rho_profile[i] * g * (delta_z / my_pipe.h));

        }
        return p_profile;
    }
    /// @brief ����� ���������� ��������� ��������
    /// @return P���
    vector<double> euler_solver_QP()
    {
        double delta_z = (my_pipe.z_L - my_pipe.z_0) / (my_pipe.n - 1);
        double speed = calc_speed(my_task.Q, my_pipe.internal_diameter);
        vector<double> p_profile = vector<double>(my_pipe.n);
        p_profile[my_pipe.n - 1] = my_task.p_L;
        for (int i = my_pipe.n - 2; i >= 0; i--)
        {
            double Re = calc_Re(speed, my_pipe.internal_diameter, my_task.nu_profile[i]);
            double hydraulic_resistance = hydraulic_resistance_isaev(Re, my_pipe.roughness);
            double tau = calc_tau(hydraulic_resistance, my_task.rho_profile[i], speed);
            p_profile[i] = p_profile[i + 1] - my_pipe.h * ((-4 / my_pipe.internal_diameter) * tau - my_task.rho_profile[i] * g * (delta_z / my_pipe.h));

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
    /// @brief ������� ������� �������
    /// @param x ������� ������
    /// @return ������� �������
    var_type residuals(const var_type& x)
    {
        my_task_parameters temp_task = task; // ��������� ���������
        temp_task.Q = x; // �� ��������� ��������� ���������� Q ��� ������ ��������� �������, ��� Q ����� ���� � ������
        euler_solver e_solver(pipe, temp_task); // ��������� ���������� ������ ������� �������

        p_profile = e_solver.euler_solver_PQ(); // ������� ������� �������� �������
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

    // ���������� ������������ ����� �������
    size_t maxLength = max({ time_rho.size(), rho_and_nu_in_0.size(), time_nu.size(), rho_and_nu_in_1.size(),
                                      time_p_in.size(), p_in.size(), time_p_out.size(), p_out.size(), Q.size(), time_Q.size() });

    // ��������� ���� ��� ������
    ofstream file(filename);

    // ���������, ������ �� ���� �������
    if (file.is_open()) {
        // ���������� ��������� ��������
        file << "Time Density; Density; Time Viscosity; Viscosity; Time Pressure In; Pressure In; Time Pressure Out; Pressure Out; Time Flow Rate; Flow Rate\n";

        // ���������� ������ �� ��������
        for (size_t i = 0; i < maxLength; ++i) {
            // ���� ������ ��������� � �������� ����� �������, ���������� ��������, ����� ���������� ������ ������
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
        // ��������� ����
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
        // ���� ���� ����������, ������� ��� ����������
        if (file.is_open()) {
            file.close();
            file.open(filename, ios::out | ios::trunc);
        }
        // ���� ���� �� ����������, ������� �����
        else {
            file.open(filename, ios::out);
        }
        // ���� ����������, �� file.is_open() ���������� False
        if (!file.is_open()) {
            file.clear();  // ������� ���� ������
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




/// @brief ����� ��� ������� quickest_ultimate_fv_solver
/// @brief ����� ��� ������� quickest_ultimate_fv_solver
class QuickWithQuasiStationaryModel : public ::testing::Test {
protected:
    // ������� ����������
    typedef quickest_ultimate_fv_solver_traits<1>::var_layer_data target_var_t;
    typedef quickest_ultimate_fv_solver_traits<1>::specific_layer specific_data_t;

    // ����: ���������� Vars + ������� ������ ��������� Specific
    typedef composite_layer_t<target_var_t, specific_data_t> layer_t;
protected:
    /// @brief ��������� �����
    pipe_properties_t pipe;
    /// @brief ������� �������
    vector<double> Q_profile;
    std::unique_ptr<PipeQAdvection> advection_model;
    std::unique_ptr<ring_buffer_t<layer_t>> buffer, buffer_nu;
protected:

    /// @brief ���������� � ������� ��� ��������� ������
    virtual void SetUp() override {
        // ���������� ������� ����� - 50��, � ����� ��������� ��� �������� ����� 1��, ��������� 700��
        simple_pipe_properties simple_pipe;
        //simple_pipe.length = 50e3;
        simple_pipe.length = 700e3; // ���� ����� 700��
        //simple_pipe.diameter = 0.7;
        simple_pipe.diameter = 0.514; // ���� ����� 700��
        //simple_pipe.dx = 100;
        simple_pipe.dx = 100; // ���� ����� 700��
        pipe = pipe_properties_t::build_simple_pipe(simple_pipe);

        std::srand(std::time(nullptr));
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<double> dis(7.0, 10.0);

        Q_profile = vector<double>(pipe.profile.getPointCount(), 0.2);
        for (size_t i = 1; i <= Q_profile.size()/2; i++)
        {
            Q_profile[i] = 0.2 - abs(0.2 * 0.001 * std::rand() / RAND_MAX);
        }

        for (size_t i = Q_profile.size() / 2 + 1; i <= Q_profile.size() - 1; i++)
        {
            Q_profile[i] = 0.1 - abs(0.1 * 0.001 * std::rand() / RAND_MAX);
        }

        advection_model = std::make_unique<PipeQAdvection>(pipe, Q_profile);
        buffer = std::make_unique<ring_buffer_t<layer_t>>(2, pipe.profile.getPointCount());
        buffer_nu = std::make_unique<ring_buffer_t<layer_t>>(2, pipe.profile.getPointCount());

        layer_t& prev = buffer->previous();
        prev.vars.cell_double[0] = vector<double>(prev.vars.cell_double[0].size(), 850);

        layer_t& prev_nu = buffer_nu->previous();
        prev_nu.vars.cell_double[0] = vector<double>(prev_nu.vars.cell_double[0].size(), 15e-6);
    }
};

/// @brief ��������� ����� ��� ������ PQ
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

public:
    /// @brief ���������� ��������� �����
    /// @param pipe ������ �� �������� �����
    /// @param oil ������ �� �������� �����
    /// @param flow �������� ������
    Pipe_model_for_PQ_t(pipe_properties_t& pipe, vector<double>& rho_profile, vector<double>& nu_profile, double flow)
        : pipe(pipe)
        , rho_profile(rho_profile)
        , nu_profile(nu_profile)
        , flow(flow)
    {

    }

    /// @brief ���������� ��������� ��������� �����
    virtual const vector<double>& get_grid() const override {
        return pipe.profile.coordinates;
    }

    /// @brief ���������� �������� ������ ����� ��
    /// @param grid_index ������������� ������ ��������� �����
    /// @param point_vector ��������� �������
    /// @return �������� ������ ����� �� � ����� point_vector
    virtual right_party_type ode_right_party(
        size_t grid_index, const var_type& point_vector) const override
    {
        double rho = rho_profile[grid_index];
        double S_0 = pipe.wall.getArea();
        double v = flow / (S_0);
        double Re = v * pipe.wall.diameter / nu_profile[grid_index];
        double lambda = pipe.resistance_function(Re, pipe.wall.relativeRoughness());
        double tau_w = lambda / 8 * rho * v * abs(v);
        /// ��������� ������� � ������ �������� �� �������� �����
        /// ����� �� ����� �� ������ �����, ����� ������� dz/dx � �������� �����
        grid_index = grid_index == 0 ? grid_index + 1 : grid_index;
        grid_index = grid_index == pipe.profile.heights.size() - 1 ? grid_index - 1 : grid_index;

        double height_derivative = (pipe.profile.heights[grid_index] - pipe.profile.heights[grid_index - 1]) /
            (pipe.profile.coordinates[grid_index] - pipe.profile.coordinates[grid_index - 1]);

        return { ((-4) / pipe.wall.diameter) * tau_w - rho * M_G * height_derivative };
    }
};

/// @brief ������ ������ � ���� �����
TEST_F(QuickWithQuasiStationaryModel, TimeRowsTask)
{
    // ���������� ������� ����� - 50��, � ����� ��������� ��� �������� ����� 1��, ��������� 700��
    simple_pipe_properties simple_pipe;
    //simple_pipe.length = 50e3;
    simple_pipe.length = 700e3; // ���� ����� 700��
    // ������� �������� �����
    pipe_properties_t pipe = pipe_properties_t::build_simple_pipe(simple_pipe);
    // ������� �������� �����
    oil_parameters_t oil;

    const double len = pipe.profile.getLength();
    //simple_pipe.diameter = 0.7;
    simple_pipe.diameter = 0.514; // ���� ����� 700��
    const auto& x = advection_model->get_grid();
    double dx = x[1] - x[0];
    double v = advection_model->getEquationsCoeffs(0, 0);

    double delta_d = 0.01, abs_roughness = 15e-6, z_0 = 50, z_L = 100,
        rho = 850, nu = 15e-6, p_0 = 6e6, p_L = 0.0;
    my_pipe_parameters my_pipe{ simple_pipe.length, simple_pipe.diameter, delta_d, abs_roughness, z_0, z_L, dx };

    my_task_parameters my_task{ my_pipe, rho, nu, p_0, p_L, {}, {}, {} };

    // ������ ������� ������ �����, [�3/�]
    double Q = calc_flow(v,simple_pipe.diameter);

    layer_t& prev = buffer->previous();
    prev.vars.cell_double[0] = vector<double>(prev.vars.cell_double[0].size(), rho);

    layer_t& prev_nu = buffer_nu->previous();
    prev_nu.vars.cell_double[0] = vector<double>(prev.vars.cell_double[0].size(), nu);

    Pipe_model_for_PQ_t pipeModel(pipe, prev.vars.cell_double[0], prev_nu.vars.cell_double[0], Q);
    // ������ �������� ��������
    double Pin = 6e6;

    // ������ ��� �������� ������� ��������
    
    int points = pipe.profile.getPointCount();
    vector<double> profile(points);
    /// ���������������� ����� ������ ��� ������ pipeModel,
    /// ������ ������� ������-������ ������������ �����,
    /// ��������� ������� Pout, 
    /// ���������� ������� ��������� � ����, �� ������� ��������� start_layer
    solve_euler_corrector<1>(pipeModel, 1, Pin, &profile);
    double pin_debug = profile[0];
    double pout_debug = profile.back();
}