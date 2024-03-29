#pragma once

namespace pde_solvers {

/// @brief �������� ������������� ��������� ������� �������������, 
/// ��� ���� ���������, ��� �������� �� ����� ������������ ���������,
/// � ����� ������� ������ ����� �������
class transport_moc_solver
{
public:
    /// @brief ����������� ������������� �������
    /// @param pipe ��������� ������������ 
    /// @param vol_flow �������� ������
    /// @param prev ���������� ����
    /// @param next ����� ����
    transport_moc_solver(const pipe_properties_t& pipe, double vol_flow,
        vector<double>& prev, vector<double>& next)
        : pipe{ pipe }
        , volumetric_flow{ vol_flow }
        , prev{ prev }
        , next{ next }
    {}

    /// @brief ������ ������ ����
    /// @param par_in �������� ��������� �����, ��������� � ������ ������������
    /// @param par_out �������� ��������� �����, ��������� � ����� ������������ ��� �������� ������� 
    void step(double par_in, double par_out)
    {
        int direction = get_eigen_value() > 0 ? 1 : -1;
        size_t start_index = direction > 0 ? 1 : (next.size()) - 2;
        size_t end_index = direction < 0 ? 0 : next.size();
        next[start_index - direction] = direction > 0 ? par_in : par_out;
        for (size_t index = start_index; index != end_index; index += direction)
        {
            next[index] = prev[index - direction];
        }
    }

    /// @brief ������ ���� �� �������, ��� ������� ������ ����� ������� (Cr = 1)
    double prepare_step() const
    {
        const std::vector<double>& grid = pipe.profile.coordinates;

        double max_eigen = get_eigen_value();
        double dx = grid[1] - grid[0];
        double courant_step = dx / max_eigen;

        return courant_step;
    }

protected:

    /// @brief ������ ������������
    const pipe_properties_t& pipe;
    /// @brief �������� ������
    const double volumetric_flow;
    /// @brief ���������� ����
    vector<double>& prev;
    /// @brief ����� ����
    vector<double>& next;


    /// @brief ������ ������������ ��������
    double get_eigen_value() const
    {
        double S = pipe.wall.getArea();
        return volumetric_flow / S;
    }
};


}