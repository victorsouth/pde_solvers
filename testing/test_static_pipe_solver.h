#pragma once


TEST(Static_Hydraulic_Solver, UseCase)
{
    typedef profile_collection_t<2> layer_variables_type;
    typedef moc_solver<2>::specific_layer layer_moc_type;

    typedef composite_layer_t<layer_variables_type, layer_moc_type> composite_layer_type;

    ring_buffer_t<composite_layer_type> buffer(2, 3);

    pipe_properties_t pipe;
    pipe.profile.coordinates = { 0, 1000, 2000 };
    pipe.profile.heights = pipe.profile.capacity = std::vector<double>(pipe.profile.coordinates.size(), 0);

    oil_parameters_t oil;
    PipeModelPGConstArea pipeModel(pipe, oil);

    profile_wrapper<double, 2> start_layer(get_profiles_pointers(buffer.current().vars.point_double));

    double Pin = 5.2e5;
    double Pout = 5e5;

    double G = solve_pipe_PP(pipeModel, Pin, Pout, &start_layer);
}

