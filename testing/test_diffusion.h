#pragma once



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



/// @brief Проверка расчета концентрации на выходе трубы, сравнение со старой моделью
TEST(DiffusionSolver, UseCase)
{
    auto simple_pipe = simple_pipe_properties::sample_district();
    auto pipe = pipe_properties_t::build_simple_pipe(simple_pipe);
    pipe.wall.equivalent_roughness = 15e-5;
    oil_parameters_t oil;
    oil.viscosity.nominal_viscosity = 6e-7;

    double v = 2.4096; 
    size_t T = 290505; //проверка алгоритма

    // Задание массива моментов времени для расчета выходного параметра
    double dt = 60;   // Шаг временного подсчёта концентрации - такой же как в методе характеристик
    size_t N = static_cast<size_t>(1.2 * T / dt); //проверка алгоритма 
    N = static_cast<size_t>(3600 * 2 / dt); //проверка алгоритма 
    vector<double> t(N);
    for (size_t i = 0; i < N; i++)
    {
        //t(i) = (i + 1 + 18480) * dt; //Конец трубы
        t[i] = (i + 1 + 4800) * dt; //Конец трубы
    }

    // точность расчета интеграла (шаг в секундах)
    // и одновременно период дискретизации для задания входных граничных условий
    double delta_t = 1; 

    double t_change = 60;
    size_t n_change = static_cast<size_t>(t_change / delta_t + 0.5) + 1;
    size_t input_size = static_cast<size_t>(t.back() / delta_t + 0.5);
    vector<double> input = 
        diffusion_transport_solver::create_boundary(850, 860, input_size, n_change, n_change);

    diffusion_transport_solver solver(pipe, oil);
    vector<double> output = solver.solve(t, delta_t, input, v, true);

    ////вывод в файлы
    std::ofstream fout;
    fout.open("time2.txt");
    fout << t;
    fout.close();

    fout.open("CL2.txt");
    fout << output;
    fout.close();
}
