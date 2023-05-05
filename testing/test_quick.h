#pragma once


TEST(QUICK_Solver, QUICK_Layer)
{
    // Профиль переменных
    typedef templated_layer<0, 1> target_var_t;

    typedef templated_layer<1, 0> specific_data_t;

    // Слой: переменных Vars + сколько угодно служебных Specific
    typedef composite_layer_t<target_var_t, specific_data_t> layer_t;

    //Текущий и предыдущий слой, каждый из которых представляет собой composite_layer (Var+Specific)
    custom_buffer_t<layer_t> buffer(4, 10);

    //Получение текущего/предыдущего слоя
    const layer_t& prev = buffer.previous();
    layer_t& next = buffer.current();

    //Движение на слой вперед 
    buffer.advance(+1);
}