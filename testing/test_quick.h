#pragma once


TEST(QUICK_Solver, QUICK_Layer)
{
    // ������� ����������
    typedef templated_layer<0, 1> target_var_t;

    typedef templated_layer<1, 0> specific_data_t;

    // ����: ���������� Vars + ������� ������ ��������� Specific
    typedef composite_layer_t<target_var_t, specific_data_t> layer_t;

    //������� � ���������� ����, ������ �� ������� ������������ ����� composite_layer (Var+Specific)
    custom_buffer_t<layer_t> buffer(4, 10);

    //��������� ��������/����������� ����
    const layer_t& prev = buffer.previous();
    layer_t& next = buffer.current();

    //�������� �� ���� ������ 
    buffer.advance(+1);
}