#pragma once


/// @brief ���������� ������ ������ ��������� �� [Leonard 1979]
TEST(UpstreamDifferencing, Develop)
{
    // ������� ����������
    typedef templated_layer<0, 1> target_var_t;
    typedef templated_layer<1, 0> specific_data_t;

    // ����: ���������� Vars + ������� ������ ��������� Specific
    typedef composite_layer_t<target_var_t, specific_data_t> layer_t;

    // ���������� ������� ����� - 50��, � ����� ��������� ��� �������� ����� 1��, ��������� 700��
    simple_pipe_properties simple_pipe;
    simple_pipe.length = 50e3;
    simple_pipe.diameter = 0.7;
    simple_pipe.dx = 1000;

    PipeProperties pipe = PipeProperties::build_simple_pipe(simple_pipe);

    vector<double> Q(pipe.profile.getPointCount(), 0.5); // ������ �� ����� ������ 0.5 �3/�
    PipeQAdvection advection_model(pipe, Q);

    //������� � ���������� ����, ������ �� ������� ������������ ����� composite_layer (Var+Specific)
    custom_buffer_t<layer_t> buffer(2, 10);

    //��������� ��������/����������� ����
    layer_t& prev = buffer.previous();
    layer_t& next = buffer.current();

    double rho_in = 860;
    double rho_out = 870;

    prev.vars.cell_double[0] = vector<double>(prev.vars.cell_double[0].size(), 850);

    auto& next_spec = std::get<0>(next.specific);
    auto& F = next_spec.point_double[0]; // ������ �� �������� �����
    auto& U = prev.vars.cell_double[0];
    auto& U_new = next.vars.cell_double[0]; //

    U[0] = rho_in;
    U[U.size() - 1] = rho_out;
    double v_in = advection_model.getEquationsCoeffs(0, U[0]);
    F[0] = v_in * U[0];
    double v_out = advection_model.getEquationsCoeffs(F.size()-1, U[U.size()-1]);
    F[F.size()-1] = v_out * U[U.size()-1];

    for (size_t cell = 0; cell < U.size(); ++cell) {
        double u = U[cell];
        double v_cell = advection_model.getEquationsCoeffs(cell, u);//�� ������ ���������, �������� � ������ ������� �� �������� �� �� ����� �������
        if (v_cell > 0) {
            size_t right_point = cell + 1;
            double v_right = advection_model.getEquationsCoeffs(right_point, u);
            F[right_point] = u * v_right;
        }
        else {
            size_t left_point = cell;
            double v_left = advection_model.getEquationsCoeffs(left_point, u);
            F[left_point] = u * v_left;
        }
        U_new[cell] = U[cell] + ((F[cell] - F[cell + 1]) / v_cell);
    }

    

    //next_spec.point_double[0].

    //�������� �� ���� ������ 
    buffer.advance(+1);
}