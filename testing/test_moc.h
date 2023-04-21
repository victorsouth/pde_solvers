#pragma once



/// @brief Базовый пример использования метода характеристик для уравнения адвекции
TEST(MOC_Solver, UseCase_Advection)
{
    simple_pipe_properties simple_pipe;
    simple_pipe.length = 50e3;
    simple_pipe.dx = 1000;

    PipeProperties pipe = PipeProperties::build_simple_pipe(simple_pipe);

    // Одна переменная, и структуры метода характеристик для нее
    typedef composite_layer_t<templated_layer<1>,
        moc_solver<1>::specific_layer> single_var_moc_t;

    custom_buffer_t<single_var_moc_t> buffer(2, pipe.profile.getPointCount());

    buffer.advance(+1);
    single_var_moc_t& prev = buffer.previous();
    single_var_moc_t& next = buffer.current();
    auto& l = prev.vars.point_double[0];
    l = vector<double>(l.size(), 1); // инициализация начальной "концентрации", равной 1

    vector<double> Q(pipe.profile.getPointCount(), 0.5); // задаем по трубе расход 0.5 м3/с
    PipeQAdvection advection_model(pipe, Q);

    moc_solver<1> solver(advection_model, prev, next);

    double dt = solver.prepare_step();
    double c_in = 2; // "концентрация" на входе
    solver.step_optional_boundaries(dt, c_in, c_in);

    auto& c_new = next.vars.point_double[0];
}


/// @brief Расчет уравнений стационарного, затем нестационарного течения слабосжимаемой жидкости
/// методом характеристик
TEST(MOC_Solver, UseCase_Waterhammer)
{
    typedef templated_layer<2> layer_variables_type;
    typedef moc_solver<2>::specific_layer layer_moc_type;

    typedef composite_layer_t<layer_variables_type, layer_moc_type> composite_layer_type;

    custom_buffer_t<composite_layer_type> buffer(2, 3);

    PipeProperties pipe;
    pipe.profile.coordinates = { 0, 1000, 2000 };
    pipe.profile.heights = pipe.profile.capacity = vector<double>(pipe.profile.coordinates.size(), 0);

    OilParameters oil;
    PipeModelPGConstArea pipeModel(pipe, oil);

    profile_wrapper<double, 2> start_layer(get_profiles_pointers(buffer.current().vars.point_double));

    double G = 400;
    double Pout = 5e5;

    solve_euler_corrector<2>(pipeModel, -1, { Pout, G }, &start_layer);

    auto left_boundary = pipeModel.const_mass_flow_equation(G);
    auto right_boundary = pipeModel.const_pressure_equation(Pout);

    vector<vector<double>> Phist, Ghist;

    for (size_t index = 0; index < 100; ++index) {
        Phist.emplace_back(buffer.current().vars.point_double[0]);
        Ghist.emplace_back(buffer.current().vars.point_double[1]);

        buffer.advance(+1);

        auto left_boundary = pipeModel.const_mass_flow_equation(G + 50);

        moc_layer_wrapper<2> moc_current(buffer.current().vars, std::get<0>(buffer.current().specific));
        moc_layer_wrapper<2> moc_previous(buffer.previous().vars, std::get<0>(buffer.previous().specific));

        moc_solver<2> solver(pipeModel, moc_previous, moc_current);
        //double dt = 0.2;
        solver.step(left_boundary, right_boundary);

    }

}

