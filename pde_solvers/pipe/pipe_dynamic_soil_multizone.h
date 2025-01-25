#pragma once

namespace pde_solvers {
;

/*
    ��������� ������ � ������� ��� ���������� ������.
    ����������� ������� ������� � ����������� �������������:
      �������� ���������, ��������� ������������� ������, ������������ ������ 
    ��� ����� ��� �������������. ������������ ������ ������� �� ������������� ���������� � ������������� �������������   

*/

enum class HeatModelVer { Legacy, V2 };


/// @brief ��������� ������������� �������� ���� ������������� ������ 
/// ��� ����������� ������
struct heat_zone_adaptation_t {
    /// @brief ��������� ���������������� ���� ��������
    double conductivity_isolation{ 1.0 };
    /// @brief ��������� ���������������� ���� ������
    double conductivity_soil{ 1.0 };
    /// @brief ��������� �������� �������� ������������ ������
    double volume_heat_capacity_soil{ 1.0 };
    /// @brief ��������� ����������� ���������� �������������� ������
    double temperature_ambient_soil{ 1.0 };
};


/// @brief ������������ ���������� �������� ������, ��������� �� ������������� ����������
struct equivalent_heat_zone_coefficients_t {
    double al1; // alpha1
    double al2; // alpha2
    double b;

    double A1;
    double A2;
    double A3;
    double A4;

    double A;
    double B;
    double C;
    double F;

    /// @brief ����������� ���������� �������������� ������
    /// ��������� ������������ � ���� ��������� ����������� � ������ ���������� ���������, 
    /// ����� ����������� ������ ���� � ������ �������������
    double temperature_ambient{ KELVIN_OFFSET - 5 };

    /// @brief ����� ���������� T ������
    double get_relaxation_time() const {
        return -1.0 / F;
    }
    /// @brief ����������� �������� ������� T �������� �� T ������
    double get_gain() const {
        return -C / F;
    }
    /// @brief �������� ����� �� ���������� ������������ ������ � �����������
    /// @param T_oil ����������� �������� � ������������
    /// @param T_soil ����������� ���������� ������������ ������
    /// @param ambient_temperature ����������� ���������� �������������� ������
    double get_heat_flow(double T_oil, double T_soil, double ambient_temperature) const {
        double DT = T_oil - ambient_temperature;
        double DTsr = T_soil - ambient_temperature;
        double q = (A + B * C) * DT + F * B * DTsr;
        return q;
        double q_st = A * DT; // ��� ���������
    }
    /// @brief �������� ����� �� ���������� �������������� ������ � �����������
    /// @param T_oil 
    /// @param ambient_temperature 
    /// @return 
    double get_heat_flow(double T_oil, double ambient_temperature) const {
        double DT = T_oil - ambient_temperature;
        double q = A * DT;
        return q;
    }

    /// @brief ������ ������� ���������� ������ ����� ������ ������� ����������� �������� �� ������
    /// @param T_oil ����������� ��������, ������������ �� ������
    /// @param ambient_temperature ����������� ��������� ������
    double get_soil_temperature_static(double T_oil, double ambient_temperature) const {
        //������ ������������ ����� ������� ����������� ������ � ����������� ������
        //double K1 = lambda * (pow(al2, 2) - 2 * log(al2) - 1) / (2 * b * log(al1) * log(al2) * (pow(al2, 2) - 1));
        double K = -C / F;

        // ���������� ����������� �������� ��� ����
        double DT1 = T_oil - ambient_temperature;
        double result = K * DT1 + ambient_temperature;

        return result;
    }

    /// @brief ������ ����������� �������� ���� ��������
    /// @param T_oil ����������� ��������
    /// @param T_soil ������� ����������� ������
    /// @param ambient_temperature ����������� ��������� ������
    double get_isolation_outer_temperature(double T_oil, double T_soil, double ambient_temperature) const {
        double DT1 = T_oil - ambient_temperature;
        double DTsr = T_soil - ambient_temperature;
        double result = (A3 + A4 * C) * DT1 + A4 * F * DTsr + ambient_temperature;
        return result;
    }
    /// @brief ������ ������� ����������� ������ �� ����� ����
    /// @param dt ��������� ���
    /// @param T_oil ����������� �������� � �����
    /// @param T_soil_old ������ ����������� ������
    /// @param ambient_temperature ����������� ��������� ������
    /// @return ������ ����������� ������ �� ����� ����
    double get_next_soil_temperature(double dt, double T_oil, double T_soil_old, double ambient_temperature) const
    {
        double DT1 = T_oil - ambient_temperature;
        double DTsr = T_soil_old - ambient_temperature;
        double dT = dt * (C * DT1 + F * DTsr);
        double T_soil_new = T_soil_old + dT;
        return T_soil_new;
    }
};



/// @brief ��������������� � �������������� ��������� 
/// ������������� ����������� �������� ������
struct equivalent_heat_zone_parameters_t {
    /// @brief ���������� ������ ����. 
    size_t coordinate_begin;
    size_t coordinate_end;

    /// @brief ���������� ������ ������������
    double r1{ std::numeric_limits<double>::quiet_NaN() };
    /// @brief ������� ������ �������������� ���� ��������
    double r1eq{ std::numeric_limits<double>::quiet_NaN() };
    /// @brief ������� ������ �������������� ���� ������
    double r2eq{ std::numeric_limits<double>::quiet_NaN() };
    /// @brief ���������������� ���� ��������
    double conductivity_isolation{ std::numeric_limits<double>::quiet_NaN() };
    /// @brief ��������������� ��������� ������
    thermophysical_properties_t soil;
    /// @brief ����������� ���������� �������������� ������
    double _temperature_ambient{ KELVIN_OFFSET - 5 };

    /// @brief ������������ ��� ��������� ������ ������������� �������
    equivalent_heat_zone_coefficients_t get_coefficients(
        const heat_zone_adaptation_t& ident) const
    {
        equivalent_heat_zone_coefficients_t result;

        result.temperature_ambient = _temperature_ambient + 100 * (ident.temperature_ambient_soil - 1);

        double soilHeatCapacityVolume = soil.density * soil.heat_capacity * ident.volume_heat_capacity_soil;
        double soilConductivity = soil.conductivity * ident.conductivity_soil;
        double isolationCoductivity = conductivity_isolation * ident.conductivity_isolation;

        // ����������� ���������������������� ������, ���������� � �������� ����������� ������� ������������
        double aSqr = soilConductivity / (soilHeatCapacityVolume * pow(r1, 2));
        double a = sqrt(aSqr);

        // ��������� ��������
        double al1 = result.al1
            = r1eq / r1; // alpha1
        double al2 = result.al2
            = r2eq / r1eq; // alpha2
        double b = result.b
            = isolationCoductivity / log(al1)
            + soilConductivity / log(al2);


        // ������������ ��������� ����� ���������� T1 � Tsr
        double& A1 = result.A1 = (pow(al2, 2) - pow(al2, 2) * log(al2) - log(al2) - 1) / (8 * aSqr * log(al2));
        double& A2 = result.A2 = (pow(al2, 2) - 2 * log(al2) - 1) / (2 * log(al2) * (pow(al2, 2) - 1));
        double& A3 = result.A3 = isolationCoductivity / (b * log(al1));
        double& A4 = result.A4 = soilConductivity * (2 * log(al2) - pow(al2, 2) + 1);

        // ������������ ��� ��������� ��������� ������

        result.A = -2 * M_PI * soilConductivity * isolationCoductivity / (b * log(al1) * log(al2));
        result.B = M_PI * soilConductivity * (2 * log(al2) - pow(al2, 2) + 1) *
            (b * log(al2) - soilConductivity) / (2 * aSqr * log(al2));
        result.C = -A3 * A2 / (A1 + A4 * A2);
        result.F = 1 / (A1 + A4 * A2);

        return result;
    }
};


/// @brief ��������� ��������������� ��������� (��� ���� ������������ � ����������� ����������)
struct primary_heat_zone_parameters_t {
    /// @brief ���������� ������ ����. 
    /// ����� ���� ������������� ���� ��������� �����, ���� ������ ������������
    size_t coordinate_begin;
    /// @brief ������� ������, ����� �������� � ��������� ����, �
    /// �������� ������ �����, ������ ���� ��������, ������ ���� ��������, �������� ����
    array<double, 4> isolation_thickness{ 0.010, 0.003, 0.050, 0.005 };
    /// @brief ���������������� ���������� ������������, ��*�-1*�-1
    /// �����, ������ ���� ��������, ������ ���� ��������, �������� ����
    array<double, 4> isolation_conductivity{ 40, 0.2, 0.035, 0.3 };
    /// @brief ������� ��������� - ������� ���� ������, �
    double depth = 1.5;
    /// @brief ��������������� ��������� ������
    thermophysical_properties_t soil;
    /// @brief ����������� ���������� �������������� ������
    double temperature_ambient{ KELVIN_OFFSET - 5 };

    /// @brief ����� ������� ����������� ���� ��������
    double get_total_isolation_thickness() const {
        return std::accumulate(isolation_thickness.begin(), isolation_thickness.end(), 0.0);
    }

    /// @brief ������ ������������ ���������������� \lambda ���� ��������
    /// @param pipe_inner_diameter 
    /// @return 
    double get_total_isolation_conductivity(double pipe_inner_diameter) const {
        double Kt_sum = get_total_isolation_HTC(pipe_inner_diameter);

        double r_inner = pipe_inner_diameter / 2;
        double r_outer = r_inner + get_total_isolation_thickness();

        double result = Kt_sum * r_inner * log(r_outer / r_inner);
        return result;
    }

    /// @brief ����������� ����������� ���� ��������, ����������� � ����������� �������� �����
    /// ����� ���������������� �������� �������� ����� ��� �������
    /// @param pipe_inner_diameter ���������� ������� ������������
    double get_total_isolation_HTC(double pipe_inner_diameter) const {
        // ����� ���������������� �������� �������� ����� ��� �������
        // �������� �� �������: invKt1 = 1 / (a1 * 2 * r1)

        // invKt1 - ��������, �������� ������������ ����������� ���� �������� (temp1), 
        // ����������� � ����������� �������� �����
        double invKt1 = 0;
        double r_prev = pipe_inner_diameter / 2;
        for (size_t i = 0; i < isolation_thickness.size(); i++)
        {
            double r_next = r_prev + isolation_thickness[i];
            invKt1 += log(r_next / r_prev) / (2 * isolation_conductivity[i]);
            r_prev = r_next;
        }
        double Kt1 = 1 / (pipe_inner_diameter * invKt1);
        return Kt1;
    }

    /// @brief ����������� ����������� ������, ����������� � � �������� �������� ���� ��������
    /// @param pipe_inner_diameter ���������� ������� ������������
    double get_soil_HTC(double pipe_inner_diameter) const {

        double r1 = pipe_inner_diameter / 2; // ���������� ������ �����
        double rn = r1 + get_total_isolation_thickness(); // ������� ������ ���� ��������
        double dn = 2 * rn; // ������� ������� ���� ��������

        double delta_h = 0.1; // ������� ��������������� (����������)
        // ���� �������� ��� ����� ����� ����� �� ��������� ����������� ������ � ���������, �

        double h0 = 2 * (depth - rn) + delta_h;

        double invKt2 = log(h0 / dn +
            sqrt(pow((h0) / dn, 2) - 1))
            / (2 * soil.conductivity);

        // invKt2 - ��������, �������� ������������ ����������� ������,
        // ����������� � �������� �������� ���� ��������
        double Kt2 = 1 / (dn * invKt2);
        return Kt2;
    }

    /// @brief ������ ������������� ����������� ������������� ����� Kt1, Kt2
    /// ������ � ������� ��������� dsn, ��������� ��� �������������
    pair<double, double> get_equivalent_heat_transfers_legacy(double pipe_inner_diameter) const
    {
        // ����������� ���������� �����������. �������� �������, ���� ������
        double a1 = 257;
        // ���������� ������ ���� ()
        double r1 = pipe_inner_diameter / 2;

        // ����������� ������������� �������������� ���� ��������
        double temp1 = 1 / (a1 * 2 * r1);
        double temp0 = temp1;

        double ds = 2 * r1;
        for (size_t i = 0; i < isolation_thickness.size(); i++)
        {
            double dsn = ds + isolation_thickness[i];
            temp1 += log(dsn / ds) / (2 * isolation_conductivity[i]);
            ds = dsn;
        }
        double Kt1 = 1 / (2 * r1 * temp1);

        double temp2 = log((2 * (depth - ds / 2) + 0.1) / ds +
            sqrt(pow((2 * (depth - ds / 2) + 0.1) / ds, 2) - 1))
            / (2 * soil.conductivity);

        double Kt2 = 1 / (ds * temp2);

        return std::make_pair(Kt1, Kt2);
    }
    /// @brief ��������� ������� ������������� ����� �������� (r1) � ������ (r2)
    /// @param Kt1 
    /// @param Kt2 
    /// ����������� ���������������� �������������� ���� ��������. 
    /// �� ��������� �� ��������, � ����������� ������ ���������������� ��������� ���� ��������
    /// ���� �����, �� ����� ����� ��������� �������������
    /// @return 
    pair<double, double> get_equivalent_radiuses_legacy(double Kt1, double Kt2,
        double pipe_inner_diameter) const
    {
        // ������������ ���. ���� ��������, ����� ����� ������� ����
        double isolatonConductivity = isolation_conductivity[2];

        // ������������� �������
        double r1 = pipe_inner_diameter / 2;
        double r1eq = r1 * exp(isolatonConductivity / (Kt1 * r1));
        double r2eq = r1eq * exp(soil.conductivity / (Kt2 * r1eq));

        return std::make_pair(r1eq, r2eq);
    }

    /// @brief ������ ���������� ������������� ������ �������� �������� � ������� �����������-�����
    /// �� ��������� � legacy, 
    ///  a) ������ ������ � ������� Kt1 ��� ����� ��������
    ///  b) ����� ���� ����������� � ��������
    ///  �) �������� ������ ������� r1eq
    ///  d) ��� ���������, ������ r2eq ���� ����������
    /// @param fluid ��������� �������������� ��������
    /// @return ���������� ������������� ������ 
    equivalent_heat_zone_parameters_t get_heat_eqivalent_model_alt(double pipe_inner_diameter) const
    {
        equivalent_heat_zone_parameters_t result;
        result.r1 = pipe_inner_diameter / 2;

        // ������������������ ���� �������� �������������� ������� �� ���� ������� ���������������� ������
        result.r1eq = result.r1 + get_total_isolation_thickness();
        result.conductivity_isolation = get_total_isolation_conductivity(pipe_inner_diameter);

        // ������������������ ������ �������������� ������� �� ���� ������� ������� ������
        double Kt2 = get_soil_HTC(pipe_inner_diameter);
        result.r2eq = result.r1eq * exp(soil.conductivity / (Kt2 * result.r1eq));
        result.soil = soil;
        result._temperature_ambient = temperature_ambient;
        return result;
    }

    /// @brief ������ ������������� ������ ��� ������������� � ������ ������� ������
    /// @param pipe_inner_diameter 
    /// @return 
    equivalent_heat_zone_parameters_t get_heat_eqivalent_model_lurie(double pipe_inner_diameter) const
    {
        equivalent_heat_zone_parameters_t result;
        result.r1 = pipe_inner_diameter / 2;

        result.conductivity_isolation = isolation_conductivity[2];

        auto [Kt1, Kt2] = get_equivalent_heat_transfers_legacy(pipe_inner_diameter);
        std::tie(result.r1eq, result.r2eq)
            = get_equivalent_radiuses_legacy(Kt1, Kt2, pipe_inner_diameter);

        result.soil = soil;
        result._temperature_ambient = temperature_ambient;

        return result;
    }

    equivalent_heat_zone_parameters_t get_heat_eqivalent_model(
        HeatModelVer model_version, double pipe_inner_diameter) const
    {
        if (model_version == HeatModelVer::Legacy) {
            return get_heat_eqivalent_model_lurie(pipe_inner_diameter);
        }
        else {
            return get_heat_eqivalent_model_alt(pipe_inner_diameter);
        }
    }

};


/// @brief ����������� � ������������� ���������������� ���������� ������
struct zoned_pipe_properties : public pipe_properties_t {
    /// @brief ��������� ��������� ��� (�������� ������������� ��������)
    vector<primary_heat_zone_parameters_t> heat_zones;
    /// @brief ��������� ������������� (�������� ������������� �� ��������� ��� � �������� �������������)
    vector<heat_zone_adaptation_t> adaptation_zones;
    /// @brief ��������� ������������� ��� 
    /// ����������� �� ��������� ����������, ����� ������ �� ��������
    /// �� ����, ������� �� ��������� ���
    vector<equivalent_heat_zone_parameters_t> eqheat_zones;
    /// @brief ������������ ������ ��� ������������ ��� 
    /// (����������� �� ���������� ������������� ��� � ���������� �������������)
    vector<equivalent_heat_zone_coefficients_t> coeff_zones;
    vector<size_t> pressure_sensor_indices;
    HeatModelVer model_version{ HeatModelVer::V2 };

    void init_pressure_sensors(vector<double> coordinates) {
        std::sort(coordinates.begin(), coordinates.end());

        double dx = profile.coordinates[1] - profile.coordinates[0];
        pressure_sensor_indices = { 0 };
        for (double x : coordinates) {
            size_t index = static_cast<int>(x / dx + 0.5);
            pressure_sensor_indices.push_back(index);
        }
        pressure_sensor_indices.push_back(profile.coordinates.size() - 1);
    }

    /// @brief ������������ ������������� �������, ������������� ������ � ������ ���������� ���������
    void init_zone_coeffs() {
        eqheat_zones = get_equivalent_zones();
        coeff_zones = get_coeff_zones(adaptation_zones);
    }



    const equivalent_heat_zone_coefficients_t& get_heat_zone_coeff(size_t index) const {
        auto iter = std::find_if(eqheat_zones.rbegin(), eqheat_zones.rend(),
            [index](const equivalent_heat_zone_parameters_t& zone)
            {
                return index >= zone.coordinate_begin;
            });

        if (iter == eqheat_zones.rend())
            throw std::logic_error("wrong zone index");

        size_t zone_index = iter - eqheat_zones.rbegin(); // ������ ���� � �����
        return *(coeff_zones.rbegin() + zone_index);
    }

    const equivalent_heat_zone_parameters_t& get_heat_zone(size_t index) const {

        auto iter = std::find_if(eqheat_zones.rbegin(), eqheat_zones.rend(),
            [index](const equivalent_heat_zone_parameters_t& zone)
            {
                return index >= zone.coordinate_begin;
            });

        if (iter == eqheat_zones.rend())
            throw std::logic_error("wrong zone index");

        return *iter;
    }

    /// @brief ���������� ��������� ������������� ��������������� �������
    vector<equivalent_heat_zone_coefficients_t> get_coeff_zones(const vector<heat_zone_adaptation_t>& adaptation_zones) const
    {
        vector<equivalent_heat_zone_coefficients_t> result;
        for (size_t index = 0; index < eqheat_zones.size(); ++index) {
            if (adaptation_zones.size() == eqheat_zones.size()) {
                result.emplace_back(eqheat_zones[index].get_coefficients(adaptation_zones[index]));
            }
            else {
                heat_zone_adaptation_t zone_ident;
                result.emplace_back(eqheat_zones[index].get_coefficients(zone_ident));
            }
        }
        return result;
    }

    /// @brief ��������� �� ��������� - ���� �������
    vector<heat_zone_adaptation_t> get_initial_ident_parameters() const {
        vector<heat_zone_adaptation_t> result(eqheat_zones.size());
        return result;
    }


    /// @brief ���������� ��������� ������������� ��������������� �������
    vector<equivalent_heat_zone_parameters_t> get_equivalent_zones() const {
        vector<equivalent_heat_zone_parameters_t> result;

        //for (const auto& zone : heat_zones) 
        for (size_t index = 0; index < heat_zones.size(); ++index)
        {
            const auto& zone = heat_zones[index];
            auto eqzone = zone.get_heat_eqivalent_model(model_version, wall.diameter);

            eqzone.coordinate_begin = zone.coordinate_begin;
            if (index + 1 < heat_zones.size()) {
                eqzone.coordinate_end = heat_zones[index + 1].coordinate_begin;
            }
            else {
                eqzone.coordinate_end = profile.getPointCount() - 1;
            }

            result.emplace_back(eqzone);
        }

        return result;
    }

};


}