#pragma once

namespace pde_solvers {
;

/*
   ������ ����������� ��� ����� ��� ����� �������� �������� ������
   ���������� ���������� ����� ���������������� (������� ����) 
     � �������� ������� (��������� ��������)
   ��������� ����� (� ���������) ��� ������ � PipeParameters
*/

/// @brief ������������ ����� ���������������� ������� ���. 
/// �������� ����������� (�� ������ ���������), ���� ��� �������������� �������������������
struct thermal_conductivity_parameters_t {
    /// @brief ������� �������������� ������ �� ����������� �������� 
    double thickness;
    /// @brief ���������������� �������������� ������
    double conductivity;
    /// @brief ����������� ���������������� UA
    /// @param diameter_inner ���������� ������� �����
    double get_heat_transfer_coefficient(double diameter_inner) const {
        const double& lambda = conductivity;
        double d_in = diameter_inner;
        double d_out = d_in + 2 * thickness;
        double Kt = 2 * lambda / (d_in  * log(d_out / d_in));
        return Kt;
    }
    /// @brief �������� ����� ����������������
    double get_heat_flow(double temperature_inner, 
        double temperature_outer, double diameter_inner) const 
    {
        const double Kt = get_heat_transfer_coefficient(diameter_inner);
        const double& Tin = temperature_inner;
        const double& Tout = temperature_outer;
        return -Kt * (Tin - Tout);
    }
};

/// @brief ��������� ��������� �������� ����� ����� ���������������� �������
struct radiative_transport_parameters_t {
    /// @brief ������� ������� �����
    double thickness;
    /// @brief ������� ������� ���������� �����
    double emissivity_factor_inner{ 0.075 };
    /// @brief ������� ������� ������� �����
    double emissivity_factor_outer{ 0.075 };
    /// @brief ����������� ������� �������
    double get_reduced_emissivity_factor(double diameter_inner) const {
        double diameter_outer = diameter_inner + 2 * thickness;
        double denum =
            1 / emissivity_factor_inner +
            diameter_inner / diameter_outer * (1 / emissivity_factor_outer - 1);
        double result = 1 / denum;
        return result;
    }
    /// @brief �������� �������� �����
    /// @param temperature_inner 
    /// @return 
    double get_heat_flow(double temperature_inner,
        double temperature_outer, double diameter_inner) const 
    {
        const double& Tinner = temperature_inner;
        const double& Touter = temperature_outer;
        constexpr double boltzmann_constant = 5.76e-8;
        double eps_red = get_reduced_emissivity_factor(diameter_inner);
        double q = -eps_red * boltzmann_constant * (pow(Tinner, 4) - pow(Touter, 4));
        return q;
    }
};

/// @brief ������� ������ ����������� ��� �����. 
/// �������� ���������� �� ����� ����� ���������, ����������� ���!
struct heat_model_general_t {
    /// @brief ����������� ���������� ������
    double ambient_temperature{ 273 };
    virtual double get_heat_flow(double T_fluid, double diameter_inner) const = 0;
};

/// @brief ��������������� ������ �����������: 
/// ���������������� � ��������
struct heat_model_conductivity_radiative_t : public heat_model_general_t {
    /// @brief ��������� ����������������
    thermal_conductivity_parameters_t conductivity;
    /// @brief ��������� ��������� �����������
    radiative_transport_parameters_t radiative;
    /// @brief ���� ���������������� � ����� ����������� (�� ���������������� � ���������)
    double conductivity_frac{ 0.04 };
    virtual double get_heat_flow(double T_fluid, double diameter_inner) const override {
        double q_cond = conductivity.get_heat_flow(T_fluid, ambient_temperature, diameter_inner);
        double q_rad = radiative.get_heat_flow(T_fluid, ambient_temperature, diameter_inner);
        double q = q_cond * conductivity_frac + q_rad * (1 - conductivity_frac);
        return q;
    }
};

}