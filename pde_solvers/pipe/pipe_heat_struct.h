#pragma once


namespace pde_solvers {
;

/// @brief Теплофизические свойства вещества (грунта, металла)
struct thermophysical_properties_t {
    /// @brief Теплопроводность \lambda, Вт*м-1*К-1
    double conductivity = 1.6; // значение для грунта
    /// @brief Удельная массовая теплоемкость C, Дж/кг-К
    double heat_capacity = 2000; // значение для грунта
    /// @brief Плотность \rho, кг/м3
    double density = 1500; // грунт
    /// @brief Температуропроводность
    double get_thermal_diffusivity() const {
        return conductivity / (density * heat_capacity);
    }
};

inline thermophysical_properties_t soil_default;
inline thermophysical_properties_t soil_sandy_fresh = { 0.299, 799.6788, 1594.329231 };
inline thermophysical_properties_t soil_sandy_saturated = { 2.198427746,   	1482.1272,     1998.47181 };
inline thermophysical_properties_t soil_clay_dry = { 0.250606936,   	481.1626678,   1603.563459 };
inline thermophysical_properties_t soil_clay_saturated = { 1.579687861,   	1389.113408,   1998.416222 };
inline thermophysical_properties_t soil_peat_dry = { 0.060491329,   	992.0339708,   305.3528885 };
inline thermophysical_properties_t soil_peat_saturated = { 0.499485549,   	1519.215498,   1105.978757 };
inline vector<pair<const char*, const thermophysical_properties_t*>> soil_list{
    {"soil default", &soil_default},
    {"sandy fresh",    &soil_sandy_fresh},
    {"sandy saturated",    &soil_sandy_saturated},
    {"clay dry",    &soil_clay_dry},
    {"clay saturated",    &soil_clay_saturated},
    {"peat dry",    &soil_peat_dry},
    {"peat saturated", &soil_peat_saturated}
};

inline const std::map<string, const thermophysical_properties_t*> soils(soil_list.begin(), soil_list.end());



}

