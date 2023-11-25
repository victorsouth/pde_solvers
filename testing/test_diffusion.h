#pragma once


/// @brief ������ ��������� �� ������ �������
double compute_integral(double ts, double tau, double xs, double C_tau, double Pe)
{
//�������� �������
    double R = (xs / (pow(ts - tau, 1.5))) * (exp((-Pe / 4) * (pow(xs - (ts - tau), 2) / (ts - tau)))) * C_tau;
    return R;
}

/// @brief ������� ������� ������������ ���������� ����� ����� ������� ���������������� �������������
/// @brief ������������� 2018 �-�� #7
double get_C_x_t(double ts, double xs, double Pe, double h, const VectorXd& C_tau)
{
    double C;
    double tau;
    double tau_p;
    C = 0;
    size_t N = C_tau.size();

    for (size_t i = 0; i < N - 1; i++)
    {
        tau = (i + 1) * h;
        tau_p = (i)*h;
        double I2 = compute_integral(ts, tau, xs, C_tau(i), Pe);
        double I1 = compute_integral(ts, tau_p, xs, C_tau(i), Pe);
        C = C + (I2 + I1) * h / 2;
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
    double L = 700000; // ������� ������ 700 ��
    double v = 2.4096;//0.3609;  // �������� ������
    double diameter = 0.514;
    double T = L / v;
    /// @brief �������������� ��������
    double nu = 6e-7;
    double Re = v * diameter / nu;
    /// �������������
    double delta = 15e-5;
    double lambda = 0.11 * pow(((68 / Re) + (delta / diameter)), 0.25);
    double K = 3.211 * sqrt(lambda) * v * diameter;//0.1;
    double Pe = v * L / K;

    double dt = 60;   // ��� ���������� �������� ������������ (���)

    // size_t len = (T + (3600*2)) / dt;  // ������� ������� (��������������� �������� �������)
    //size_t len = (3600 * 2) / dt;  // ������� ������� (��������������� �������� �������)
    double N = 1.2 * T / dt;
    // ������� ������� ������� ��� ����������� �������
    VectorXd t = VectorXd::Zero(N);
    for (size_t i = 0; i < N; i++)
    {
        t(i) = (i + 1 + 18480) * dt; //����� �����
    }

    // ������� ������� ��� ���������� ����������
    VectorXd x = VectorXd::Zero(N);
    for (size_t i = 0; i < N; i++)
    {
        // x ������ ���������� ��� ������� ����� �������� ������ ���������� t � ����������� �� �������� ��� �������
        x(i) = L; // ����� �����
    }

    // ������ ������������
    //double Pe = 1438000; // ����� �����
    double del = 1; // �������� ���������� ������ �������������� (���)
    VectorXd CL = compute_C_experiment(t, x, Pe, T, L, del, N);

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