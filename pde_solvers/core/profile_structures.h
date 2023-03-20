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

/// @brief Возвращает указатели на нулевые элементы профилей, хранящиеся в array
/// @tparam DataType 
/// @param profiles 
/// @return 
template <size_t Dimension, typename DataType>
inline array<vector<DataType>*, Dimension> get_profiles_pointers(array<vector<DataType>, Dimension>& profiles)
{
    return create_array<Dimension>([&](int dimension) { return &profiles[dimension]; });
}


/// @brief Собирает из составной профиль
/// Из скалярных профилей собрать векторный профиль: double -> array<double, Dim>
/// Из векторных профилей собрать матричный профиль: array<double, Dim> -> array<array<double, Dim>, Dim>
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

    /// @brief Длина профиля (для совместимости с std::vector)
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


/// @brief Шаблонный слой, определяющий нужное количество профилей по точкам и ячейкам
/// Используется для генерации слоя расчетных переменных и слоев вспомогательных структур
/// Скалярный профиль на точках
/// Скалярный профиль на ячейках
/// Векторный профиль (заданной размерности)
template <size_t PointScalar, size_t CellScalar = 0,
    size_t PointVector = 0, size_t PointVectorDimension = 0,
    size_t CellVector = 0, size_t CellVectorDimension = 0>
struct templated_layer
{
    typedef typename fixed_system_types<PointVectorDimension>::var_type point_vector_type;
    typedef typename fixed_system_types<CellVectorDimension>::var_type cell_vector_type;


    /// @brief Список скалярных профилей на границах ячеек
    array<vector<double>, PointScalar> point_double;
    /// @brief Список скалярных профилей в ячейках
    array<vector<double>, CellScalar> cell_double;
    /// @brief Список векторных профилей на границах ячеек
    array<vector<point_vector_type>, PointVector> point_vector;
    /// @brief Список векторных профилей в ячейках
    array<vector<cell_vector_type>, CellVector> cell_vector;

    templated_layer(size_t point_count)
        : point_double{ array_maker<vector<double>, PointScalar>::make_array(vector<double>(point_count)) }
        , point_vector{ array_maker<vector<point_vector_type>, PointVector>::make_array(vector<point_vector_type>(point_count)) }
        , cell_double{ array_maker<vector<double>, CellScalar>::make_array(vector<double>(point_count - 1)) }
        , cell_vector{ array_maker<vector<cell_vector_type>, CellVector>::make_array(vector<cell_vector_type>(point_count)) }
    {

    }
};

/// @brief Составной слой, включающий в себя слой переменных и слои со специальными структурам
/// @tparam Variables 
/// @tparam ...Ts 
template <typename VarLayer, typename... SpecificLayers>
struct composite_layer_t {
    /// @brief Целевые переменные
    VarLayer vars;
    /// @brief Все специальные структуры
    std::tuple<SpecificLayers...> specific;

    /// @brief Возвращает ссылку на специальные структуры с заданным номером (SpecificDataNumber) в кортеже
    template <unsigned SpecificDataNumber>
    auto& get_specific() {
        return std::get<SpecificDataNumber>(specific);
    }

    composite_layer_t(size_t point_count)
        : vars(point_count)
        , specific((sizeof(SpecificLayers), point_count)...)
    {}
};
