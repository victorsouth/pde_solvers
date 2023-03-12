#pragma once



TEST(MOC_Solver, UseCase)
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

