#pragma once

using std::vector;




template <typename T>
inline static T linear_interpolation(T f1, T f2, double p)
{
    return f1 * (1 - p) + f2 * p;
}

template <typename T>
inline static T _interpolate(T* values, double interpolation_offset)
{
    double p = interpolation_offset;
    if (p == 0)
        return values[0];
    else if (p > 0)
        return linear_interpolation<T>(values[0], values[+1], p);
    else
        return linear_interpolation<T>(values[-1], values[0], 1 + p);
}

/// @brief ¬озвращает указатели на нулевые элементы профилей, хран¤щиес¤ в array
/// @tparam DataType 
/// @param profiles 
/// @return 
template <size_t Dimension, typename DataType>
inline array<vector<DataType>*, Dimension> get_profiles_pointers(array<vector<DataType>, Dimension>& profiles)
{
    return create_array<Dimension>([&](int dimension) { return &profiles[dimension]; });
}


/// @brief —обирает из составной профиль
/// »з скал¤рных профилей собрать векторный профиль: double -> array<double, Dim>
/// »з векторных профилей собрать матричный профиль: array<double, Dim> -> array<array<double, Dim>, Dim>
template <typename T, size_t Dimension>
class profile_wrapper {
protected:
    array<vector<T>*, Dimension> profiles;
public:
    typedef typename fixed_system_types<Dimension>::var_type vector_type;

    const size_t n;
public:
    profile_wrapper(vector<double>& profile);

    profile_wrapper(array<vector<T>*, Dimension>* profiles)
        : profiles(*profiles)
        , n(profiles->front()->size())
    {

    }

    profile_wrapper(array<vector<T>*, Dimension> profiles)
        : profiles(profiles)
        , n(profiles.front()->size())
    {

    }

    /// @brief ƒлина профил¤ (дл¤ совместимости с std::vector)
    size_t size() const {
        return n;
    }

public:
    array_ref<T, Dimension> operator()(size_t profile_index)
    {
        array_ref<T, Dimension> result([&](int dimension) {
            return &profiles[dimension]->at(profile_index);
            });

        return result;
    }
    array_ref<T, Dimension> operator[](size_t profile_index)
    {
        return (*this)(profile_index);
    }
protected:

public:
    T interpolate_dimension(size_t dimension, size_t profile_index, double frac_offset) const
    {
        if (profile_index == 0 && frac_offset < 0) {
            throw std::runtime_error("profile_index == 0 && frac_offset < 0");
        }
        else if (profile_index == n - 1 && frac_offset > 0) {
            throw std::runtime_error("profile_index == profile_size - 1 && frac_offset > 0");
        }
        if (dimension >= Dimension) {
            throw std::runtime_error("dimension >= Dimension");
        }

        auto& profile_ref = *profiles[dimension];
        T* value_ptr = &profile_ref[profile_index];
        return _interpolate(value_ptr, frac_offset);
    }

    vector_type interpolate(size_t profile_index, double frac_offset) const
    {
        if (profile_index == 0 && frac_offset < 0) {
            throw std::runtime_error("profile_index == 0 && frac_offset < 0");
        }
        else if (profile_index == n - 1 && frac_offset > 0) {
            throw std::runtime_error("profile_index == profile_size - 1 && frac_offset > 0");
        }

        vector_type result;
        for (size_t index = 0; index < Dimension; ++index) {
            auto& profile_ref = *profiles[index];
            T* value_ptr = &profile_ref[profile_index];
            result[index] = _interpolate(value_ptr, frac_offset);
        }
        return result;
    }


    double& operator()(size_t dimension, size_t profile_index)
    {
        return profiles[dimension]->at(profile_index);
    }
    const vector<T>& profile(size_t profile_index) const {
        return *profiles[profile_index];
    }
    vector<T>& profile(size_t profile_index) {
        return *profiles[profile_index];
    }

};

template <>
inline profile_wrapper<double, 1>::vector_type
profile_wrapper<double, 1>::interpolate(size_t profile_index, double frac_offset) const
{
    if (profile_index == 0 && frac_offset < 0) {
        throw std::runtime_error("profile_index == 0 && frac_offset < 0");
    }
    else if (profile_index == n - 1 && frac_offset > 0) {
        throw std::runtime_error("profile_index == profile_size - 1 && frac_offset > 0");
    }

    auto& profile_ref = *profiles[0];
    double* value_ptr = &profile_ref[profile_index];
    vector_type result = _interpolate(value_ptr, frac_offset);

    return result;
}


template <>
inline profile_wrapper<double, 1>::profile_wrapper(vector<double>& profile)
    : n(profile.size())
{
    profiles[0] = &profile;
}

template <>
inline profile_wrapper<double, 1>::profile_wrapper(array<vector<double>*, 1>* profiles)
    : profiles(*profiles)
    , n(profiles->front()->size())
{

}


/// @brief Ўаблонный слой, определ¤ющий нужное количество профилей по точкам и ¤чейкам
/// »спользуетс¤ дл¤ генерации сло¤ расчетных переменных и слоев вспомогательных структур
/// —кал¤рный профиль на точках
/// —кал¤рный профиль на ¤чейках
/// ¬екторный профиль (заданной размерности)
template <size_t PointScalar, size_t CellScalar = 0,
    size_t PointVector = 0, size_t PointVectorDimension = 0,
    size_t CellVector = 0, size_t CellVectorDimension = 0>
struct templated_layer
{
    typedef typename fixed_system_types<PointVectorDimension>::var_type point_vector_type;
    typedef typename fixed_system_types<CellVectorDimension>::var_type cell_vector_type;


    /// @brief —писок скал¤рных профилей на границах ¤чеек
    array<vector<double>, PointScalar> point_double;
    /// @brief —писок скал¤рных профилей в ¤чейках
    array<vector<double>, CellScalar> cell_double;
    /// @brief —писок векторных профилей на границах ¤чеек
    array<vector<point_vector_type>, PointVector> point_vector;
    /// @brief —писок векторных профилей в ¤чейках
    array<vector<cell_vector_type>, CellVector> cell_vector;

    templated_layer(size_t point_count)
        : point_double{ array_maker<vector<double>, PointScalar>::make_array(vector<double>(point_count)) }
        , point_vector{ array_maker<vector<point_vector_type>, PointVector>::make_array(vector<point_vector_type>(point_count)) }
        , cell_double{ array_maker<vector<double>, CellScalar>::make_array(vector<double>(point_count - 1)) }
        , cell_vector{ array_maker<vector<cell_vector_type>, CellVector>::make_array(vector<cell_vector_type>(point_count)) }
    {

    }
};

/// @brief —оставной слой, включающий в себ¤ слой переменных и слои со специальными структурам
/// @tparam Variables 
/// @tparam ...Ts 
template <typename VarLayer, typename... SpecificLayers>
struct composite_layer_t {
    /// @brief ÷елевые переменные
    VarLayer vars;
    /// @brief ¬се специальные структуры
    std::tuple<SpecificLayers...> specific;

    /// @brief ¬озвращает ссылку на специальные структуры с заданным номером (SpecificDataNumber) в кортеже
    template <unsigned SpecificDataNumber>
    auto& get_specific() {
        return std::get<SpecificDataNumber>(specific);
    }

    composite_layer_t(size_t point_count)
        : vars(point_count)
        , specific((sizeof(SpecificLayers), point_count)...)
    {}
};
