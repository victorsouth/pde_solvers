#pragma once


/// @brief Расчет интеграла по методу трапецй
double compute_integral(double ts, double tau, double xs, double C_tau, double Pe)
{
    /*if (C_tau == 0)
        return 0;
    if (ts == 0)
        return 0;*/

        //double a1 = pow(ts - tau, 1.5);
        //double A = xs / a1;
        //double b1 = (-Pe / 4) * (pow(xs - (ts - tau), 2) / (ts - tau));

        //double B = exp(b1);
        //double R1 = (A * B) * C_tau;
        //return R1;

         //Исходная формула
    double R = (xs / (pow(ts - tau, 1.5))) * (exp((-Pe / 4) * (pow(xs - (ts - tau), 2) / (ts - tau)))) * C_tau;
    return R;
}

/// @brief Формула расчета концентрации Дидковская Новый метод расчета многопродуктовых магистральных
/// @brief трубопроводов 2018 ф-ла #7
double get_C_x_t(double ts, double xs, double Pe, double h, const VectorXd& C_tau)
{
    double C;
    double tau;
    double tau_p;
    C = 0;
    size_t N = C_tau.size();
    //size_t N = ts / h;

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

/// Эксперимент с изменением входной концентрации 
VectorXd compute_C_experiment(VectorXd t, const VectorXd& x, double Pe, size_t T, size_t L, double del, size_t len)
{
    double h = del / T; // Шаг по тау в интеграле метода 
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
                C0_func(j) = 0.5;   // Концентрация скачков изменяется с 0 до 0,5 после 3000-ой секунды
            }
        }

        CL(i) = get_C_x_t(t_p, x_p, Pe, h, C0_func);
    }

    return CL;
}

// Функция создания линейно возрастающего скачка по концентраций (вектор значений)
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

/// @brief Расчет движения партии для заданных моментов времени и заданных координат
VectorXd compute_C_experiment2(
    const VectorXd& t, const VectorXd& C0_params,
    double v, double L, double K, double delta_t)
{
    size_t N = t.size();
    VectorXd CL = VectorXd::Zero(N); // результат расчета - значение
    size_t n_change = static_cast<size_t>(C0_params(2) / delta_t + 0.5) + 1; // Индекс времени скачка по параметру на входе

    double T = L / v;
    double Pe = v * L / K;

    /*if (t.size() != C0_func.size())
        throw logic_error("t and C0_func must have equal size");*/

    double x_p = 1; // относительная координата
    double h = delta_t / T; // Шаг по тау в интеграле метода 

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

/// @brief Старый тест расчета концентрации на выходе трубы (из prac4)
TEST(TransportCourant, Test2_Concentration)
{
    size_t dt = 60;   // Шаг временного подсчёта концентрации (сек)
    size_t L = 400000; // Участок длиной 400 км
    size_t T = 1108320; //Время движения по участку трубы L, (Расчет от длины трубы и скорости продукта)
    // double v = 0.3609;  // Скорость потока

    // size_t len = (T + (3600*2)) / dt;  // Глубина расчёта (рассматриваемый горизонт времени)
    size_t len = (3600 * 2) / dt;  // Глубина расчёта (рассматриваемый горизонт времени)

    // Задание массива времени или конкретного времени
    VectorXd t = VectorXd::Zero(len);
    for (size_t i = 0; i < len; i++)
    {
        t(i) = (i + 1 + 18480) * dt; //Конец трубы
    }

    // Задание профиля или конкретная координаты
    VectorXd x = VectorXd::Zero(len);
    for (size_t i = 0; i < len; i++)
    {
        // x задано постоянным при расчете можно наоборот задать постоянным t в зависимости от динамики или профиля
        x(i) = L; // Конец трубы
    }

    // Расчет концентрации
    double Pe = 1438000; // Число Пекле
    double del = 1; // Точность численного метода интегрирования (сек)
    VectorXd CL = compute_C_experiment(t, x, Pe, T, L, del, len);

    //вывод в файлы
    std::ofstream fout;
    fout.open("time.txt");
    fout << t;
    fout.close();

    fout.open("CL.txt");
    fout << CL;
    fout.close();
}

/// @brief Проверка расчета концентрации на выходе трубы, сравнение со старой моделью
TEST(TransportCourant, Test3_Concentration_Compare)
{
    double L = 400000; // Участок длиной 400 км
    double v = 0.3609;  // Скорость потока
    //size_t T = 1108320; //Время движения по участку трубы L, (Расчет от длины трубы и скорости продукта)
    double T = L / v;

    double K = 0.1;
    // Восстановление коэффициента диффузии по данным из статьи
    //double Pe = 1438000; // Число Пекле
    //double L = 400000; // Участок длиной 400 км
    //double K = 0.1;// v* L / Pe;

    double Pe = v * L / K;

    double dt = 60;   // Шаг временного подсчёта концентрации - такой же как в методе характеристик

    double t_change = 3000;
    //size_t n_change = static_cast<size_t>(t_change / dt + 0.5); // Момент времени скачка по параметру на входе

    // Количество шагов для 1.5 периода времени прохождения партии по трубе (чтобы увидеть динамику)
    size_t N = 120;// 1.5 * T / dt;
    //size_t N =  1.5 * T / dt;

    VectorXd CinParams(3);
    CinParams << 0, 0.5, t_change; // Параметры скачка - с чего на что и когда
    //VectorXd Cin = create_boundary(0, 0.5, N, n_change, n_change);


    // Задание массива времени или конкретного времени
    VectorXd t = VectorXd::Zero(N);
    for (size_t i = 0; i < N; i++)
    {
        t(i) = (i + 1 + 18480) * dt; //Конец трубы
    }

    double delta_t = 1; // точность расчета интеграла (шаг в секундах)
    VectorXd CL = compute_C_experiment2(t, CinParams, v, L, K, delta_t);

    ////вывод в файлы
    std::ofstream fout;
    fout.open("time2.txt");
    fout << t;
    fout.close();

    fout.open("CL2.txt");
    fout << CL;
    fout.close();
}