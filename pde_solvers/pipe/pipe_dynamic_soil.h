#pragma once

namespace pde_solvers {
;

/*
  ��������� ������ � ������� ��� ���������� ������.
   - ��������� ������������� ���������� ������
   - �������� ��������� ������� "�����, ��������, �����"
   - ������������� ������ �������: ������������� ������� �����, ������������ �������
   - ��������������� ��������� ������������, � �������� ��������� �������� ������
   - ���������� �� �������� ������ �����-����������
*/

/// @brief ��������� ������������� ���������� ������������� ������
struct thermal_model_ident_parameters_t {
    /// @brief ���������������� ���� ��������
    double conductivity_isolation{ std::numeric_limits<double>::quiet_NaN() };
    /// @brief ���������������� ���� ������
    double conductivity_soil{ std::numeric_limits<double>::quiet_NaN() };
    /// @brief �������� �������� ������������ ������
    double volume_heat_capacity_soil{ std::numeric_limits<double>::quiet_NaN() };
};



/// @brief ��������� �������� ������ ������� "�����, ��������, �����". 
/// ���� ���� ������
struct pipe_heat_model_t {
    /// @brief ������� ������, ����� �������� � ��������� ����, �
    /// ������ ���� ��������, ������ ���� ��������, �������� ����
    array<double, 4> layerThickness{ 0.010, 0.003, 0.050, 0.005 };

    /// @brief ���������������� ���������� ������������, ��*�-1*�-1
    /// �����, ������ ���� ��������, ������ ���� ��������, �������� ����
    array<double, 4> thermalConductivity{ 40, 0.2, 0.035, 0.3 };

    /// @brief ����������� ���������� �����, �
    double ambientTemperature{ KELVIN_OFFSET - 5 };

    // @brief ������� ��������� - ������� ���� ������, �
    double depthOfOccurrence = 1.5;
    /// @brief ��������������� �������� ������
    thermophysical_properties_t soil;
    /// @brief ����������� ����������� � ���������� ������
    double ambient_heat_transfer{ 0 };

    /// @brief ������� ����������� ��������� ������
    ///  q = -Kt(T - T_ground)
    /// �������� ����� 2017, ���. 57 ���������� �������� �������
    /// @param temperature ����������� ������ � �����
    /// @return �������� �����
    double get_heat_flow_newton(double temperature) const {
        double q = -ambient_heat_transfer * (temperature - ambientTemperature);
        return q;
    }
    /// @brief ���������� ��������� ������ ���������� �������������, ������������ ��������� ��������
    thermal_model_ident_parameters_t get_ident_parameters() const {
        thermal_model_ident_parameters_t result;
        result.conductivity_isolation = thermalConductivity[2];
        result.conductivity_soil = soil.conductivity;
        result.volume_heat_capacity_soil = soil.density * soil.heat_capacity;
        return result;
    }
    /// @brief �������� ��������� �������������
    /// @param ident_parameters 
    void set_ident_parameters(const thermal_model_ident_parameters_t& ident_parameters) {
        thermalConductivity[2] = ident_parameters.conductivity_isolation;
        soil.conductivity = ident_parameters.conductivity_soil;
        soil.heat_capacity = ident_parameters.volume_heat_capacity_soil / soil.density;
    }

};

/// @brief ��������� ������������� ������
struct equivalent_thermal_model_t {
    /// @brief ����������� ����������� � ������������� ���� ��������
    double Kt1;
    /// @brief ����������� ����������� � ������������� ���� ������
    double Kt2;
    /// @brief ������� ������ �������������� ���� ��������
    double r1eq;
    /// @brief ������� ������ �������������� ���� ������
    double r2eq;

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
    double get_heat_flow(double T_oil, double T_soil, double ambient_temperature) {
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
    double get_heat_flow(double T_oil, double ambient_temperature) {
        double DT = T_oil - ambient_temperature;
        double q = A * DT;
        return q;
    }
};


/// @brief ����������� � �������� ���������� �������
struct pipe_noniso_properties_t : public pipe_properties_t {
    /// @brief �������� ������
    pipe_heat_model_t heat;



    /// @brief ������ ������������� ����������� ������������� ����� Kt1, Kt2
    /// @param fluid ��������� ������
    /// @return 
    pair<double, double> get_equivalent_heat_transfers(const oil_parameters_t& fluid) const
    {
        // ����������� ���������� �����������
        double a1 = fluid.heat.internalHeatTransferCoefficient;

        // ���������� ������ ���� ()
        //double r1 = (PipeP.wall.diameter - PipeP.heat.layerThickness[0]) / 2;
        double r1 = wall.diameter / 2;

        // ����������� ������������� �������������� ���� ��������
        double temp1 = 1 / (a1 * 2 * r1);
        double temp0 = temp1;

        double ds = 2 * r1;
        for (size_t i = 0; i < heat.layerThickness.size(); i++)
        {
            double dsn = ds + heat.layerThickness[i];
            temp1 += log(dsn / ds) / (2 * heat.thermalConductivity[i]);
            ds = dsn;
        }
        double Kt1 = 1 / (2 * r1 * temp1);

        double temp2 = log((2 * (heat.depthOfOccurrence - ds / 2) + 0.1) / ds +
            sqrt(pow((2 * (heat.depthOfOccurrence - ds / 2) + 0.1) / ds, 2) - 1))
            / (2 * heat.soil.conductivity);

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
    pair<double, double> get_equivalent_radiuses(double Kt1, double Kt2) const
    {
        // ������������ ���. ���� ��������, ����� ����� ������� ����
        double isolatonConductivity = heat.thermalConductivity[2];

        // ������������� �������
        double r1 = wall.diameter / 2;
        double r1eq, r2eq;
        r1eq = r1 * exp(isolatonConductivity / (Kt1 * r1));
        r2eq = r1eq * exp(heat.soil.conductivity / (Kt2 * r1eq));

        return std::make_pair(r1eq, r2eq);
    }

    /// @brief ������ ���������� ������������� ������ �������� �������� � ������� �����������-�����
    /// @param fluid ��������� �������������� ��������
    /// @return ���������� ������������� ������ 
    equivalent_thermal_model_t get_heat_eqivalent_model(const oil_parameters_t& fluid,
        const thermal_model_ident_parameters_t& ident = thermal_model_ident_parameters_t()
    ) const
    {
        equivalent_thermal_model_t result;
        double r1 = wall.diameter / 2;

        std::tie(result.Kt1, result.Kt2) = get_equivalent_heat_transfers(fluid);
        std::tie(result.r1eq, result.r2eq) = get_equivalent_radiuses(result.Kt1, result.Kt2);

        double isolatonConductivity = heat.thermalConductivity[2];
        double soilConductivity = heat.soil.conductivity;
        double soilHeatCapacityVolume = heat.soil.density * heat.soil.heat_capacity;
        if (isfinite(ident.conductivity_isolation))
            isolatonConductivity = ident.conductivity_isolation;
        if (isfinite(ident.conductivity_soil))
            soilConductivity = ident.conductivity_soil;
        if (isfinite(ident.volume_heat_capacity_soil))
            soilHeatCapacityVolume = ident.volume_heat_capacity_soil;


        // ����������� ���������������������� ������, ���������� � �������� ����������� ������� ������������
        double aSqr = soilConductivity / (soilHeatCapacityVolume * pow(r1, 2));
        double a = sqrt(aSqr);

        // ��������� ��������
        double al1 = result.al1
            = result.r1eq / r1; // alpha1
        double al2 = result.al2
            = result.r2eq / result.r1eq; // alpha2
        double b = result.b
            = isolatonConductivity / log(al1)
            + soilConductivity / log(al2);


        // ������������ ��������� ����� ���������� T1 � Tsr
        double& A1 = result.A1;
        double& A2 = result.A2;
        double& A3 = result.A3;
        double& A4 = result.A4;
        A1 = (pow(al2, 2) - pow(al2, 2) * log(al2) - log(al2) - 1) / (8 * aSqr * log(al2));
        A2 = (pow(al2, 2) - 2 * log(al2) - 1) / (2 * log(al2) * (pow(al2, 2) - 1));
        A3 = isolatonConductivity / (b * log(al1));
        A4 = soilConductivity * (2 * log(al2) - pow(al2, 2) + 1);

        // ������������ ��� ��������� ��������� ������

        result.A = -2 * M_PI * soilConductivity * isolatonConductivity / (b * log(al1) * log(al2));
        result.B = M_PI * soilConductivity * (2 * log(al2) - pow(al2, 2) + 1) *
            (b * log(al2) - soilConductivity) / (2 * aSqr * log(al2));
        result.C = -A3 * A2 / (A1 + A4 * A2);
        result.F = 1 / (A1 + A4 * A2);

        return result;
    }

};



/// @brief ������ ������� ������� ���������� ������ �� �����, ��������� 2021, ������� 14
inline void compute_soil_mean_temperature_distribution(
    const pipe_noniso_properties_t& pipe,
    const oil_parameters_t& oil,
    const vector<double>& temperature_oil,
    vector<double>* _result,
    const thermal_model_ident_parameters_t& ident = thermal_model_ident_parameters_t()
)
{
    vector<double>& temperature = *_result;

    double lambda = pipe.heat.thermalConductivity[2]; //����� �����������, ������ �������� ������ ���� ����

    equivalent_thermal_model_t heat_dynamic_model = pipe.get_heat_eqivalent_model(oil, ident);

    double al1 = heat_dynamic_model.al1;
    double al2 = heat_dynamic_model.al2;
    double b = heat_dynamic_model.b;

    //������ ������������ ����� ������� ����������� ������ � ����������� ������
    double K = lambda * (pow(al2, 2) - 2 * log(al2) - 1) / (2 * b * log(al1) * log(al2) * (pow(al2, 2) - 1));
    double K1 = -heat_dynamic_model.C / heat_dynamic_model.F;

    // ������ ������������� ���
    for (size_t i = 0; i < temperature.size(); i++)
    {
        // ������� ���������� ������������ ������� T��(r) �������� ����������� T��� = T��(r_2)
        double DT1 = temperature_oil[i] - pipe.heat.ambientTemperature;
        temperature[i] = K1 * DT1 + pipe.heat.ambientTemperature;
    }
}

/// @brief �������� ������� ���������� ������ �� �����, ��������� 2021, ���. 6
inline void compute_new_fluid_temperature_distribution(
    //const PipeHeatInflowConstArea& heat_pde,
    const pipe_noniso_properties_t& pipe,
    const oil_parameters_t& oil,
    double mass_flow,
    double Tin,
    double dt,
    const vector<double>& _heat_flow,
    const vector<double>& temperature_oil_old,
    vector<double>* temperature_oil_new)
{
    const vector<double>& qT = _heat_flow;
    vector<double>& T1 = *temperature_oil_new;

    double d = pipe.wall.diameter;
    double S_0 = pipe.wall.getArea();
    double density = oil.density.nominal_density;
    double C = oil.heat.HeatCapacity;


    T1[0] = Tin;
    for (size_t i = 0; i < T1.size() - 1; i++)
    {
        //double s1 = heat_pde.getSourceTerm(i, temperature_oil_old[i]);

        double v = mass_flow / (S_0 * density);
        double Re = v * pipe.wall.diameter / oil.viscosity(temperature_oil_old[i]);
        double lambda = hydraulic_resistance_shifrinson(Re, pipe.wall.relativeRoughness());

        double q = 0.5 * (qT[i] + qT[i + 1]);
        double s = ((4 * q / (density * C * d)) + lambda * pow(v, 3) / (2 * C * d));
        double dT = dt * s;
        T1[i + 1] = temperature_oil_old[i] + dT;
    }
}

/// @brief �������� ������� ���������� ������. 
/// ������� ���������� �� ��������� ������
inline void compute_new_fluid_temperature_distribution2(
    //const PipeHeatInflowConstArea& heat_pde,
    const pipe_noniso_properties_t& pipe,
    const oil_parameters_t& oil,
    double mass_flow,
    double Tin,
    double dt,
    const vector<double>& _heat_flow,
    const vector<double>& temperature_oil_old,
    const vector<double>& Tsr,
    vector<double>* temperature_oil_new,
    const thermal_model_ident_parameters_t& ident
)
{
    const vector<double>& qT = _heat_flow;
    vector<double>& T1 = *temperature_oil_new;

    double d = pipe.wall.diameter;
    double S_0 = pipe.wall.getArea();
    double density = oil.density.nominal_density;
    double C_oil = oil.heat.HeatCapacity;

    equivalent_thermal_model_t heat_dynamic_model = pipe.get_heat_eqivalent_model(oil, ident);
    double A = heat_dynamic_model.A;
    double B = heat_dynamic_model.B;
    double C = heat_dynamic_model.C;
    double F = heat_dynamic_model.F;

    T1[0] = Tin;
    for (size_t i = 0; i < T1.size() - 1; i++)
    {
        auto right_party = [&](size_t i, double T) {
            double v = mass_flow / (S_0 * density);
            double Re = v * pipe.wall.diameter / oil.viscosity(T);
            double lambda = hydraulic_resistance_shifrinson(Re, pipe.wall.relativeRoughness());

            double DT1 = T - pipe.heat.ambientTemperature;
            double DTsr = Tsr[i] - pipe.heat.ambientTemperature;

            double q = (A + B * C) * DT1 + F * B * DTsr;

            double s = ((4 * q / (density * C_oil * d)) + lambda * pow(v, 3) / (2 * C_oil * d));
            return s;
        };

        //double s1 = heat_pde.getSourceTerm(i, temperature_oil_old[i]);
        //double q = 0.5 * (qT[i] + qT[i + 1]);
        //double s = ((4 * q / (density * C * d)) + lambda * pow(v, 3) / (2 * C * d));

        double rp1 = right_party(i, temperature_oil_old[i]);
        double Tnext_estimate = temperature_oil_old[i] + dt * rp1;
        double rp2 = right_party(i + 1, Tnext_estimate);

        T1[i + 1] = temperature_oil_old[i] + dt * 0.5 * (rp1 + rp2);
    }
}


/// @brief �������� ������� ������� ���������� ������ �� �����, ��������� 2021, ���. 6
inline void compute_new_soil_mean_temperature_distribution(
    const pipe_noniso_properties_t& pipe,
    const oil_parameters_t& oil,
    double dt,
    const vector<double>& _temperature,
    const vector<double>& _mean_dirt_temperature,
    vector<double>* _result,
    const thermal_model_ident_parameters_t& ident = thermal_model_ident_parameters_t()

)
{
    const vector<double>& T1 = _temperature;
    const vector<double>& Tsr_old = _mean_dirt_temperature;
    vector<double>& Tsr = *_result;

    equivalent_thermal_model_t heat_dynamic_model = pipe.get_heat_eqivalent_model(oil, ident);
    double C = heat_dynamic_model.C;
    double F = heat_dynamic_model.F;

    for (size_t i = 0; i < Tsr_old.size(); i++)
    {
        double DT1 = T1[i] - pipe.heat.ambientTemperature;
        double DTsr = Tsr_old[i] - pipe.heat.ambientTemperature;
        double dT = dt * (C * DT1 + F * DTsr);
        Tsr[i] = Tsr_old[i] + dT;
    }
}

/// @brief ������ ��������� ������ �� �����, ��������� 2021, ������� 13
inline void compute_heat_flow_distribution(
    const pipe_noniso_properties_t& pipe,
    const oil_parameters_t& oil,
    const vector<double>& temperature,
    const vector<double>& mean_amb_temperature,
    vector<double>* _qT,
    const thermal_model_ident_parameters_t& ident = thermal_model_ident_parameters_t()
)
{
    const vector<double>& T1 = temperature;
    const vector<double>& Tsr = mean_amb_temperature;
    vector<double>& qT = *_qT;

    equivalent_thermal_model_t heat_dynamic_model = pipe.get_heat_eqivalent_model(oil, ident);
    double A = heat_dynamic_model.A;
    double B = heat_dynamic_model.B;
    double C = heat_dynamic_model.C;
    double F = heat_dynamic_model.F;

    // ������ ������������� qT
    for (size_t i = 0; i < qT.size(); i++)
    {
        //double qt_furie = pipe.heat.get_heat_flow_newton(T1[i]); // ��� ���������

        double DT1 = T1[i] - pipe.heat.ambientTemperature;
        double DTsr = Tsr[i] - pipe.heat.ambientTemperature;

        double qt_lurie_st = A * DT1;

        double qt_lurie = (A + B * C) * DT1 + F * B * DTsr;
        qT[i] = qt_lurie;
    }
}

/// @brief ������ ��������� ����������� ����� ��������� � �������
inline void compute_isolation_boundary_temperature(
    const pipe_noniso_properties_t& pipe,
    const oil_parameters_t& oil,
    const vector<double>& temperature,
    const vector<double>& mean_amb_temperature,
    vector<double>* _Tiso,
    const thermal_model_ident_parameters_t& ident = thermal_model_ident_parameters_t()
)
{
    const vector<double>& T1 = temperature;
    const vector<double>& Tsr = mean_amb_temperature;
    vector<double>& Tiso = *_Tiso;

    equivalent_thermal_model_t heat_dynamic_model = pipe.get_heat_eqivalent_model(oil, ident);
    double A3 = heat_dynamic_model.A3;
    double A4 = heat_dynamic_model.A4;
    double C = heat_dynamic_model.C;
    double F = heat_dynamic_model.F;

    // ������ ������������� qT
    for (size_t i = 0; i < Tiso.size(); i++)
    {
        double DT1 = T1[i] - pipe.heat.ambientTemperature;
        double DTsr = Tsr[i] - pipe.heat.ambientTemperature;

        Tiso[i] = (A3 + A4 * C) * DT1 + A4 * F * DTsr + pipe.heat.ambientTemperature;
    }
}


}