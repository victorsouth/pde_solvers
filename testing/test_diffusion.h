#pragma once

/// @brief ��������������� ���������
/// ���������� ����� ����� ������� ���������������� �������������
/// ������������� 2018 �-�� #7
/// @param ts ������� ������ ��������� (������������ �����)
/// @param tau ������� ������������ �����
/// @param xs ������������ ����������, ��� ������� � �������� ����� ����� ��������� ������� ��������� �� ������������� �������
/// @param C_tau �������� ��������� ��� �������� ������������� �������
/// @param Pe ����� ����� Pe = v*L/K
/// @return �������� ���������������� ���������
double function_under_integral(double ts, double tau, double xs, double C_tau, double Pe)
{
    //�������� �������
    double R = (xs / (pow(ts - tau, 1.5))) * (exp((-Pe / 4) * (pow(xs - (ts - tau), 2) / (ts - tau)))) * C_tau;
    return R;
}

/// @brief ������ �������� �������� ��������� ��� �������������� ������� ������� ��� �������� ����������
/// ��������� �������������� ������� (7) �� 
/// ���������� ����� ����� ������� ���������������� ������������� ������������� 2018
/// @param ts �������� ����� � ������������ �����
/// @param xs �������� ���������� � ������������ �����
/// @param Pe ����� �����
/// @param h ��� ������� ��������� �� ������������� �������
/// @param C_tau ������� ������� �� �����
/// @return �������� �������� ��������� 
double get_C_x_t(double ts, double xs, double Pe, double h, const VectorXd& C_tau)
{
    size_t N = C_tau.size();
    double C = 0;
    for (size_t i = 0; i < N - 1; i++)
    {
        double tau = (i + 1) * h;
        double tau_p = (i)*h;
        double I2 = function_under_integral(ts, tau, xs, C_tau(i), Pe);
        double I1 = function_under_integral(ts, tau_p, xs, C_tau(i), Pe);
        C = C + (I2 + I1) * h / 2; // ����� ��������
        //C = C + I1 * h; // ����� ���������������, ������ ��� 
    }
    C = C * sqrt(Pe) / (2 * sqrt(M_PI));
    return C;
}


/// ����������� � ���������� ������� ������������ 
VectorXd compute_C_experiment(VectorXd t, const VectorXd& x, double Pe, size_t T, size_t L, double del, size_t len)
{
    double h = del / T; // ��� �� ��� � ��������� ������ 
    VectorXd CL = VectorXd::Zero(len);

    for (size_t i = 0; i < len; i++)
    {
        double t_p = t(i) / T;
        double x_p = x(i) / L;

        // VectorXd C_linear = VectorXd::Zero(len);
        double N = t_p / h;
        double length_C0 = t_p / h;
        VectorXd C0_func = VectorXd::Zero(length_C0);

        for (size_t j = 0; j < N - 1; j++)
        {
            double t_j = j * del;
            if (t_j > 3000)
            {
                C0_func(j) = 0.5;   // ������������ ������� ���������� � 0 �� 0,5 ����� 3000-�� �������
            }
        }

        CL(i) = get_C_x_t(t_p, x_p, Pe, h, C0_func);
    }

    return CL;
}

/// @brief ������ ���� ������� ������������ �� ������ ����� (�� prac4)
TEST(TransportCourant, Test2_Concentration)
{
    size_t dt = 60;   // ��� ���������� �������� ������������ (���)
    size_t L = 400000; // ������� ������ 400 ��
    size_t T = 1108320; //����� �������� �� ������� ����� L, (������ �� ����� ����� � �������� ��������)
    // double v = 0.3609;  // �������� ������

    // size_t len = (T + (3600*2)) / dt;  // ������� ������� (��������������� �������� �������)
    size_t len = (3600 * 2) / dt;  // ������� ������� (��������������� �������� �������)
    double Pe = 1438000; // ����� �����
    // size_t len = (T + (3600*2)) / dt;  // ������� ������� (��������������� �������� �������)
    //size_t len = (3600 * 2) / dt;  // ������� ������� (��������������� �������� �������)
    double N = 1.2 * T / dt;
    // ������� ������� ������� ��� ����������� �������
    VectorXd t = VectorXd::Zero(len);
    for (size_t i = 0; i < len; i++)
    {
        t(i) = (i + 1 + 18480) * dt; //����� �����
    }

    // ������� ������� ��� ���������� ����������
    VectorXd x = VectorXd::Zero(len);
    for (size_t i = 0; i < len; i++)
    {
        // x ������ ���������� ��� ������� ����� �������� ������ ���������� t � ����������� �� �������� ��� �������
        x(i) = L; // ����� �����
    }

    // ������ ������������
    //double Pe = 1438000; // ����� �����
    double del = 1; // �������� ���������� ������ �������������� (���)
    VectorXd CL = compute_C_experiment(t, x, Pe, T, L, del, len);

    //����� � �����
    std::ofstream fout;
    fout.open("time.txt");
    fout << t;
    fout.close();

    fout.open("CL.txt");
    fout << CL;
    fout.close();
}


// ������� �������� ������� ������������� ������ �� ������������ (������ ��������)
inline VectorXd create_boundary(double initial_value, double final_value,
    size_t length, size_t start_change, size_t end_change)
{
    VectorXd C0 = VectorXd::Zero(length);
    for (size_t index = 0; index < start_change; ++index) {
        C0(index) = initial_value;
    }

    end_change = std::max(start_change, end_change);

    for (size_t index = start_change; index < end_change; ++index) {
        C0(index) = initial_value + (final_value - initial_value) * (index - start_change) / (end_change - start_change);
    }

    for (size_t index = end_change; index < length; ++index) {
        C0(index) = final_value;
    }
    return C0;

}


/// @brief ������ �������� ������ ��� �������� �������� ������� � �������� ���������
VectorXd compute_C_experiment2(
    const VectorXd& t, const VectorXd& C0_params,
    double v, double L, double K, double delta_t)
{
    size_t N = t.size();
    VectorXd CL = VectorXd::Zero(N); // ��������� ������� - ��������
    size_t n_change = static_cast<size_t>(C0_params(2) / delta_t + 0.5) + 1; // ������ ������� ������ �� ��������� �� �����

    double T = L / v;
    double Pe = v * L / K;

    /*if (t.size() != C0_func.size())
        throw logic_error("t and C0_func must have equal size");*/

    double x_p = 1; // ������������� ����������
    double h = delta_t / T; // ��� �� ��� � ��������� ������ 

#pragma omp parallel for
    for (int i = 0; i < static_cast<int>(N); i++)
    {
        //       wcout << L"Calculating i = " << i << std::endl;
        double t_p = t(i) / T;
        VectorXd C0_func = create_boundary(0, C0_params(1) - C0_params(0), t_p / h, n_change, n_change);
        CL(i) = get_C_x_t(t_p, x_p, Pe, h, C0_func) + C0_params(0);
    }

    return CL;
}

/// @brief �������� ������� ������������ �� ������ �����, ��������� �� ������ �������
TEST(TransportCourant, Test3_Concentration_Compare)
{
    double L = 700000; // ������� ������ 700 ��
    double v = 2.4096;//0.3609;  // �������� ������
    double diameter = 0.514;
    //size_t T = 1108320; //����� �������� �� ������� ����� L, (������ �� ����� ����� � �������� ��������)
    double T = L / v;
    /// @brief �������������� ��������
    double nu = 6e-7;
    double Re = v * diameter / nu;
    /// �������������
    double delta = 15e-5;
    double lambda = 0.11 * pow(((68 / Re) + (delta / diameter)),0.25);



    double K = 3.211 * sqrt(lambda) * v * diameter;//0.1;
    // �������������� ������������ �������� �� ������ �� ������
    //double Pe = 1438000; // ����� �����
    //double L = 400000; // ������� ������ 400 ��
    //double K = 0.1;// v* L / Pe;

    double Pe = v * L / K;

    double dt = 60;   // ��� ���������� �������� ������������ - ����� �� ��� � ������ �������������

    double t_change = 60;
    //size_t n_change = static_cast<size_t>(t_change / dt + 0.5); // ������ ������� ������ �� ��������� �� �����

    // ���������� ����� ��� 1.5 ������� ������� ����������� ������ �� ����� (����� ������� ��������)
    //size_t N = 120;// 1.5 * T / dt;
    size_t N =  1.2 * T / dt;

    VectorXd CinParams(3);
    CinParams << 0, 0.5, t_change; // ��������� ������ - � ���� �� ��� � �����
    //VectorXd Cin = create_boundary(0, 0.5, N, n_change, n_change);


    // ������� ������� ������� ��� ����������� �������
    VectorXd t = VectorXd::Zero(N);
    for (size_t i = 0; i < N; i++)
    {
        //t(i) = (i + 1 + 18480) * dt; //����� �����
        t(i) = (i + 1 + 18480) * dt; //����� �����
    }

    double delta_t = 1; // �������� ������� ��������� (��� � ��������)
    VectorXd CL = compute_C_experiment2(t, CinParams, v, L, K, delta_t);

    ////����� � �����
    std::ofstream fout;
    fout.open("time2.txt");
    fout << t;
    fout.close();

    fout.open("CL2.txt");
    fout << CL;
    fout.close();
}




template <typename ContainerType, typename StreamType>
inline StreamType& send_stream(StreamType& o, const ContainerType& x)
{
    for (typename ContainerType::const_iterator i = x.begin(); i != x.end(); i++) {
        if (i != x.begin())
            o << "; ";
        o << std::fixed << (*i);
    }
    return o;
}

template <typename StreamType>
inline StreamType& operator<<(StreamType& o, const vector<double>& x) {
    return send_stream<vector<double>>(o, x);
}


class diffusion_transport_solver
{
private:
    const pipe_properties_t& pipe;
    const oil_parameters_t& oil;

public:
    diffusion_transport_solver(
        const pipe_properties_t& pipe,
        const oil_parameters_t& oil
    )
        : pipe(pipe)
        , oil(oil)
    {

    }
private:
    /// @brief ��������������� ���������
    /// ���������� ����� ����� ������� ���������������� �������������
    /// ������������� 2018 �-�� #7
    /// @param ts ������� ������ ��������� (������������ �����)
    /// @param tau ������� ������������ �����
    /// @param xs ������������ ����������, ��� ������� � �������� ����� ����� ��������� ������� ��������� �� ������������� �������
    /// @param C_tau �������� ��������� ��� �������� ������������� �������
    /// @param Pe ����� ����� Pe = v*L/K
    /// @return �������� ���������������� ���������
    static double function_under_integral2(double ts, double tau, double xs, double C_tau, double Pe)
    {
        //�������� �������
        double R = (xs / (pow(ts - tau, 1.5))) * (exp((-Pe / 4) * (pow(xs - (ts - tau), 2) / (ts - tau)))) * C_tau;
        return R;
    }

    /// @brief ������ �������� �������� ��������� ��� �������������� ������� ������� ��� �������� ����������
    /// ��������� �������������� ������� (7) �� 
    /// ���������� ����� ����� ������� ���������������� ������������� ������������� 2018
    /// @param ts �������� ����� � ������������ �����
    /// @param xs �������� ���������� � ������������ �����
    /// @param Pe ����� �����
    /// @param h ��� ������� ��������� �� ������������� �������
    /// @param N ���������� ����� ������� ������� �� �����, ������������ � �������
    /// @param input ������� ������� �� �����
    /// @return �������� �������� ��������� � ����� (ts, xs)
    static double get_C_x_t2(double ts, double xs, double Pe, double h, size_t N, const vector<double>& input)
    {
        double C = 0;
        if (N > input.size())
        {
            throw std::runtime_error("Wrong input length");
        }

        for (size_t i = 0; i < N - 1; i++)
        {
            double tau = (i + 1) * h;
            double tau_p = (i)*h;
            double I2 = function_under_integral2(ts, tau, xs, input[i], Pe);
            double I1 = function_under_integral2(ts, tau_p, xs, input[i], Pe);
            C = C + (I2 + I1) * h / 2; // ����� ��������
            //C = C + I1 * h; // ����� ���������������, ������ ��� 
        }
        C = C * sqrt(Pe) / (2 * sqrt(M_PI));
        return C;
    }


    static double calc_diffusive_transport(double t, double x, double delta_t,
        double v, double L, double K,
        const vector<double>& input)
    {
        double T = L / v;
        double Pe = v * L / K;

        double ts = t / T; // ������������ �����
        double xs = x / L; // ������������ ����������

        double delta_ts = delta_t / T; // ��� �������������� � ������������ �������

        size_t N = static_cast<size_t>(t / delta_t + 0.5);

        return get_C_x_t2(ts, xs, Pe, delta_ts, N, input);
    }


public:
    
    static vector<double> create_boundary(double initial_value, double final_value,
        size_t length, size_t start_change, size_t end_change)
    {
        vector<double> result(length);
        for (size_t index = 0; index < start_change; ++index) {
            result[index] = initial_value;
        }

        end_change = std::max(start_change, end_change);

        for (size_t index = start_change; index < end_change; ++index) {
            result[index] = initial_value + (final_value - initial_value) * (index - start_change) / (end_change - start_change);
        }

        for (size_t index = end_change; index < length; ++index) {
            result[index] = final_value;
        }
        return result;

    }



    static double calc_diffusion_coefficient(
        const pipe_properties_t& pipe, const oil_parameters_t& oil, double v
    )
    {
        double nu = oil.viscosity.nominal_viscosity;
        double Re = v * pipe.wall.diameter / nu;
        double relative_roughness = pipe.wall.relativeRoughness();
        double lambda = hydraulic_resistance_altshul(Re, relative_roughness);
        double K = 3.211 * sqrt(lambda) * v * pipe.wall.diameter;//0.1
        return K;
    }

    /// @brief ������ ���������� ���� ��������� �� ������ ������������ ��� �������� ������
    /// ��������� ������ ���� ����� ��������, ������������� ����������� 
    /// ������� ������������ �������� ������� ��� �������� ���������� (����� ��� ���������)
    /// ����������� ����������� ������������� K ������� �� ��������, ������� ��������� ��� ������
    /// @param t_output ������� �������, ��� ������� ��������� �������� �� ������ ������������
    /// @param delta_t ������ ������������� �������� ���������� ����, �� �� ������������ � ���������
    /// @param input ��������� ��� ��������� �� �����
    /// @param v �������� ������
    /// @param use_offset_trick ������������ ����� �������� � �������� ���������� ��������� 
    /// (��� ���������, �.�. ������ �������)
    /// @return �������� ��������� ��� ��������� ��� �������� ������� t_output
    vector<double> solve(
        const vector<double>& t_output,
        double delta_t, vector<double> input,
        double v, bool use_offset_trick)
    {
        double pipe_length = pipe.profile.getLength();
        double K = calc_diffusion_coefficient(pipe, oil, v);

        double offset = 0;
        if (use_offset_trick) {
            offset = input[0];
            // ��������� input � ���������� ������������ ���������� ��������
            for (double& in : input) {
                in -= offset;
            }
        }


        vector<double> output(t_output.size());

#pragma omp parallel for
        for (int i = 0; i < static_cast<int>(t_output.size()); i++)
        {
            double ti = t_output[i];
            double x = pipe_length;
            double value_at_ti = calc_diffusive_transport(ti, x, delta_t, v, pipe_length, K, input);
            output[i] = value_at_ti;
        }

        if (use_offset_trick) {
            // ���������, ��� 0 �� ������ ������������� �������� offset
            for (double& out : output) {
                out += offset;
            }
        }

        return output;
    }

};

/// @brief �������� ������� ������������ �� ������ �����, ��������� �� ������ �������
TEST(DiffusionSolver, UseCase)
{
    auto simple_pipe = simple_pipe_properties::sample_district();
    auto pipe = pipe_properties_t::build_simple_pipe(simple_pipe);
    pipe.wall.equivalent_roughness = 15e-5;
    oil_parameters_t oil;
    oil.viscosity.nominal_viscosity = 6e-7;

    double v = 2.4096; 
    size_t T = 290505; //�������� ���������

    // ������� ������� �������� ������� ��� ������� ��������� ���������
    double dt = 60;   // ��� ���������� �������� ������������ - ����� �� ��� � ������ �������������
    size_t N = 1.2 * T / dt; //�������� ��������� 
    N = (3600 * 2) / dt; //�������� ��������� 
    vector<double> t(N);
    for (size_t i = 0; i < N; i++)
    {
        //t(i) = (i + 1 + 18480) * dt; //����� �����
        t[i] = (i + 1 + 4800) * dt; //����� �����
    }

    // �������� ������� ��������� (��� � ��������)
    // � ������������ ������ ������������� ��� ������� ������� ��������� �������
    double delta_t = 1; 

    double t_change = 60;
    size_t n_change = static_cast<size_t>(t_change / delta_t + 0.5) + 1;
    size_t input_size = static_cast<size_t>(t.back() / delta_t + 0.5);
    vector<double> input = 
        diffusion_transport_solver::create_boundary(850, 860, input_size, n_change, n_change);

    diffusion_transport_solver solver(pipe, oil);
    vector<double> output = solver.solve(t, delta_t, input, v, true);

    ////����� � �����
    std::ofstream fout;
    fout.open("time2.txt");
    fout << t;
    fout.close();

    fout.open("CL2.txt");
    fout << output;
    fout.close();
}
