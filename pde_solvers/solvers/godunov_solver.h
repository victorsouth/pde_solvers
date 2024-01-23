#pragma once

namespace pde_solvers {

/// @brief ¬спомогательные структуры данных, необходимые дл¤ расчетной задачи по методу √одунова
template <size_t Dimension>
struct godunov_task_traits 
{
    typedef profile_collection_t<Dimension/*переменные*/, 0, 0, 0, 0, 0> var_layer_data;
    typedef profile_collection_t<
        /* альфы слева и справа от границ */ Dimension * 2,
        /* значени¤ параметров, собственные числа в ¤чейках, линейна¤ реконструкци¤ в ¤чейках: */ Dimension * 3,
        /* значени¤ параметра ERP (слева, справа), решение задачи –имана, потоки на границах ¤чеек: */ 4, Dimension,
        /*собств. векторы в ¤чейках: */ Dimension,
        /*размерность собств. векторов в ¤чейках */ Dimension> specific_layer;
};

template <size_t Dimension>
struct godunov_layer_wrapper : layer_wrapper<Dimension> {
    typedef array<array<double, Dimension>, Dimension> matrix_type;
    typedef array<double, Dimension> vector_type;
    typedef typename godunov_task_traits<Dimension>::var_layer_data var_layer_data;
    typedef typename godunov_task_traits<Dimension>::specific_layer specific_layer_data;

    profile_wrapper<double, Dimension> point_values;
    //profile_wrapper<double, Dimension> interface_values;
    //profile_wrapper<double, Dimension> interface_flux_values;

    /// @brief —ила волны (альфа) слева от границ ¤чеек
    profile_wrapper<double, Dimension> strength_left;
    /// @brief —ила волны (альфа) справа от границ ¤чеек
    profile_wrapper<double, Dimension> strength_right;

    vector<vector_type>& erp_left;
    vector<vector_type>& erp_right;
    vector<vector_type>& interface_values;
    vector<vector_type>& interface_flux_values;

    profile_wrapper<double, Dimension> beta;

    profile_wrapper<double, Dimension> cell_values;
    profile_wrapper<double, Dimension> cell_eigenvalues;
    profile_wrapper<array<double, Dimension>, Dimension> cell_eigenvec;

    godunov_layer_wrapper(
        var_layer_data& vars,
        specific_layer_data& specific
    )
        : point_values(get_profiles_pointers(vars.point_double))
        //, interface_values(create_array<Dimension>([&](int dimension) { return &specific.point_double[dimension]; }))
        //, interface_flux_values(create_array<Dimension>([&](int dimension) { return &specific.point_vector[dimension]; }))
        , erp_left(specific.point_vector[0])
        , erp_right(specific.point_vector[1])
        , interface_values(specific.point_vector[2])
        , interface_flux_values(specific.point_vector[3])
        //, beta(specific.point_vector[2])
        , beta(create_array<Dimension>([&](int dimension) { return &specific.cell_double[dimension + 2 * Dimension]; }))
        , cell_values(create_array<Dimension>([&](int dimension) { return &specific.cell_double[dimension]; }))
        , cell_eigenvalues(create_array<Dimension>([&](int dimension) { return &specific.cell_double[dimension + Dimension]; }))
        , cell_eigenvec(get_profiles_pointers(specific.cell_vector))
        , strength_left(create_array<Dimension>([&](int dimension) { return &specific.point_double[dimension]; }))
        , strength_right(create_array<Dimension>([&](int dimension) { return &specific.point_double[dimension + Dimension]; }))
    {
        //, eigenval()

    }
};

//TODO error with eqs(eqs)
//template <typename... Eqs>
//struct equation_collector : fixed_system_t<sizeof...(Eqs)>
//{
//    std::tuple<Eqs&...> eqs;
//    equation_collector(const Eqs&... eqs)
//        : eqs(eqs)
//    {}
//};

}