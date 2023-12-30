#pragma once

TEST(MOC_Solver, TemplatedLayer1)
{
    profile_collection_t<1> tl1(10); // 1 профиль с точками размерности 10
    profile_collection_t<2> tl2(10); // 2 профиля с точками размерности 10

    profile_collection_t<0, 1> tl_cell1(10); // 1 профиль с ячейками размерности 9 (10 точек = 9 ячеек)
    profile_collection_t<0, 2> tl_cell2(10); // 2 профиля с ячейками размерности 9 (10 точек = 9 ячеек)
}

TEST(MOC_Solver, MOC_Layer)
{
    //Слой переменных
    typedef profile_collection_t<1> var_layer; 
    
    //Слой служебных данных - 1 буфер под точки, 1 буфер под вектор 9 (специфики метода хар-к)
    moc_solver<1>::specific_layer moc_layer(10); 

    //Обобщенный - слой переменных Vars + сколько угодно служебных Specific
    composite_layer_t<var_layer,
        moc_solver<1>::specific_layer> composite_layer(10);

    //Текущий и предыдущий слой, каждый из которых представляет собой composite_layer (Var+Specific)
    ring_buffer_t<composite_layer_t<profile_collection_t<1>, moc_solver<1>::specific_layer>> buffer(2, 10);

    //Получение текущего/предыдущего слоя
    const composite_layer_t<var_layer, moc_solver<1>::specific_layer>& prev = buffer.previous();
    composite_layer_t<var_layer, moc_solver<1>::specific_layer>& next = buffer.current();

    //Жвижение на слой вперед 
    buffer.advance(+1);
}

TEST(MOC_Solver, MOC_Layer_Refactor)
{
    // Профиль переменных
    typedef profile_collection_t<1> target_var_t;

    // Слой: переменных Vars + сколько угодно служебных Specific
    typedef composite_layer_t<target_var_t, moc_solver<1>::specific_layer> layer_t;

    //Текущий и предыдущий слой, каждый из которых представляет собой composite_layer (Var+Specific)
    ring_buffer_t<layer_t> buffer(2, 10);

    //Получение текущего/предыдущего слоя
    const layer_t& prev = buffer.previous();
    layer_t& next = buffer.current();

    //Движение на слой вперед 
    buffer.advance(+1);
}


/// @brief Базовый пример использования метода характеристик для уравнения адвекции
TEST(MOC_Solver, UseCase_Advection)
{
    // Упрощенное задание трубы - 50км, с шагом разбиения для расчтной сетки 1км, диаметром 700мм
    simple_pipe_properties simple_pipe;
    simple_pipe.length = 50e3;
    simple_pipe.diameter = 0.7;
    simple_pipe.dx = 1000;

    pipe_properties_t pipe = pipe_properties_t::build_simple_pipe(simple_pipe);

    // Одна переменная, и структуры метода характеристик для нееm
    typedef composite_layer_t<profile_collection_t<1>, moc_solver<1>::specific_layer> single_var_moc_t;

    ring_buffer_t<single_var_moc_t> buffer(2, pipe.profile.getPointCount());

    auto& rho_initial = buffer.previous().vars.point_double[0];
    rho_initial = vector<double>(rho_initial.size(), 850); // инициализация начальной плотности

    buffer.advance(+1);
    single_var_moc_t& prev = buffer.previous();
    single_var_moc_t& next = buffer.current();

    vector<double> Q(pipe.profile.getPointCount(), -0.5); // задаем по трубе расход 0.5 м3/с
    PipeQAdvection advection_model(pipe, Q);

    moc_solver<1> solver(advection_model, 
        prev.vars.point_double[0], next.vars.point_double[0], prev.specific);

    double dt = solver.prepare_step();
    double rho_in = 840; // плотность нефти, закачиваемой на входе трубы при положительном расходе
    double rho_out = 860; // плотность нефти, закачиваемой с выхода трубы при отрицательном расходе
    solver.step_optional_boundaries(dt, rho_in, rho_out);

    auto& c_new = next.vars.point_double[0];
}

/// @brief Проблемно-ориентированный слой
/// Задача про плотность и вязкость
struct density_viscosity_layer
{
    vector<double> density;
    vector<double> viscosity;
    moc_solver<1>::specific_layer moc_specific;

    density_viscosity_layer(size_t point_count)
        : density(point_count)
        , viscosity(point_count)
        , moc_specific(point_count)
    {

    }
    static moc_layer_wrapper<1> get_density_moc_wrapper(density_viscosity_layer& layer)
    {
        return moc_layer_wrapper<1>(layer.density, layer.moc_specific);
    }
    static moc_layer_wrapper<1> get_viscosity_moc_wrapper(density_viscosity_layer& layer)
    {
        return moc_layer_wrapper<1>(layer.viscosity, layer.moc_specific);
    }
};



/// @brief Базовый пример использования метода характеристик для уравнения адвекции
TEST(MOC_Solver, UseCase_Advection2)
{
    // Упрощенное задание трубы - 50км, с шагом разбиения для расчтной сетки 1км, диаметром 700мм
    simple_pipe_properties simple_pipe;
    simple_pipe.length = 50e3;
    simple_pipe.diameter = 0.7;
    simple_pipe.dx = 1000;

    pipe_properties_t pipe = pipe_properties_t::build_simple_pipe(simple_pipe);
    vector<double> Q(pipe.profile.getPointCount(), -0.5); // задаем по трубе расход 0.5 м3/с
    PipeQAdvection advection_model(pipe, Q);

    ring_buffer_t<density_viscosity_layer> buffer(2, pipe.profile.getPointCount());

    auto& rho_initial = buffer[0].density;
    auto& viscosity_initial = buffer[0].viscosity;
    rho_initial = vector<double>(rho_initial.size(), 850); // инициализация начальной плотности
    viscosity_initial = vector<double>(viscosity_initial.size(), 1e-5); // инициализация начальной плотности


    {
        auto density_buffer = buffer.get_custom_buffer(&density_viscosity_layer::get_density_moc_wrapper);
        moc_solver<1> solver(advection_model, density_buffer);

        double dt = solver.prepare_step();
        double rho_in = 840; // плотность нефти, закачиваемой на входе трубы при положительном расходе
        double rho_out = 860; // плотность нефти, закачиваемой с выхода трубы при отрицательном расходе
        solver.step_optional_boundaries(dt, rho_in, rho_out);
    }

    {
        auto viscosity_buffer = buffer.get_custom_buffer(&density_viscosity_layer::get_viscosity_moc_wrapper);
        moc_solver<1> solver(advection_model, viscosity_buffer);

        double dt = solver.prepare_step();
        double visc_in = 2e-5; 
        double visc_out = 0.5e-5;
        solver.step_optional_boundaries(dt, visc_in, visc_out);
    }


    auto& curr = buffer[0];


}


/// @brief Расчет уравнений стационарного, затем нестационарного течения слабосжимаемой жидкости
/// методом характеристик
TEST(MOC_Solver, UseCase_Waterhammer)
{
    typedef profile_collection_t<2> layer_variables_type;
    typedef moc_solver<2>::specific_layer layer_moc_type;

    typedef composite_layer_t<layer_variables_type, layer_moc_type> composite_layer_type;

    ring_buffer_t<composite_layer_type> buffer(2, 3);

    pipe_properties_t pipe;
    pipe.profile.coordinates = { 0, 1000, 2000 };
    pipe.profile.heights = pipe.profile.capacity = vector<double>(pipe.profile.coordinates.size(), 0);

    oil_parameters_t oil;
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

