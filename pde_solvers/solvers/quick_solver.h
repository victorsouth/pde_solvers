#pragma once

namespace pde_solvers {

/// @brief Описание типов данных для метода конечных объемов на основе upstream differencing 
template <size_t Dimension>
struct upstream_fv_solver_traits
{
    typedef profile_collection_t<0, Dimension/*переменные - ячейки*/, 0, 0, 0, 0> var_layer_data;
    typedef profile_collection_t<Dimension /*потоки F*/, 0,
        0, 0,
        0, 0> specific_layer;
};

/// @brief Описание типов данных для метода конечных объемов на основе QUICK 
template <size_t Dimension>
struct quick_fv_solver_traits
{
    typedef profile_collection_t<0, Dimension/*переменные - ячейки*/, 0, 0, 0, 0> var_layer_data;
    typedef profile_collection_t<Dimension /*потоки F*/, 0,
        0, 0,
        0, 0> specific_layer;
};

/// @brief Описание типов данных для метода конечных объемов на основе QUICKEST 
template <size_t Dimension>
struct quickest_fv_solver_traits
{
    typedef profile_collection_t<0, Dimension/*переменные - ячейки*/, 0, 0, 0, 0> var_layer_data;
    typedef profile_collection_t<Dimension /*потоки F*/, 0,
        0, 0,
        0, 0> specific_layer;
};

/// @brief Описание типов данных для метода конечных объемов на основе QUICKEST-ULTIMATE 
template <size_t Dimension>
struct quickest_ultimate_fv_solver_traits
{
    typedef profile_collection_t<0, Dimension/*переменные - ячейки*/, 0, 0, 0, 0> var_layer_data;
    typedef profile_collection_t<Dimension /*потоки F*/, 0,
        0, 0,
        0, 0> specific_layer;
};

/// @brief Солвер на основе upstream differencing, только для размерности 1!
/// [Leonard 1979]
class upstream_fv_solver {
public:
    typedef typename upstream_fv_solver_traits<1>::var_layer_data var_layer_data;
    typedef typename upstream_fv_solver_traits<1>::specific_layer specific_layer;
    typedef typename fixed_system_types<1>::matrix_type matrix_type;
    typedef typename fixed_system_types<1>::var_type vector_type;
protected:
    /// @brief ДУЧП
    pde_t<1>& pde;
    /// @brief Сетка, полученная от ДУЧП
    const vector<double>& grid;
    /// @brief Количество точек сетки
    const size_t n;
    /// @brief Предыдущий слой переменных
    const var_layer_data& prev_vars;
    /// @brief Новый (рассчитываемый) слой переменных
    var_layer_data& curr_vars;
    /// @brief Предыдущий специфический слой (сейчас не нужен! нужен ли в будущем?)
    const specific_layer& prev_spec;
    /// @brief Текущий специфический слой
    specific_layer& curr_spec;
public:
    /// @brief Конструктор для буфера в котором простой слой:
    /// когда в слое только один блок целевых переменных и один блок служебных данных
    /// Из буфера берется current() и previous()
    /// @param pde ДУЧП
    /// @param buffer Буфер слоев
    upstream_fv_solver(pde_t<1>& pde,
        ring_buffer_t<composite_layer_t<var_layer_data, specific_layer>>& buffer)
        : upstream_fv_solver(pde, buffer.previous(), buffer.current())
    {}

    /// @brief Конструктор для простых слоев - 
    /// когда в слое только один блок целевых переменных и один блок служебных данных
    /// @param pde ДУЧП
    /// @param prev Предыдущий слой (уже рассчитанный)
    /// @param curr Следующий (новый), для которого требуется сделать расчет
    upstream_fv_solver(pde_t<1>& pde,
        const composite_layer_t<var_layer_data, specific_layer>& prev,
        composite_layer_t<var_layer_data, specific_layer>& curr)
        : pde(pde)
        , grid(pde.get_grid())
        , n(pde.get_grid().size())
        , prev_vars(prev.vars)
        , curr_vars(curr.vars)
        , prev_spec(std::get<0>(prev.specific))
        , curr_spec(std::get<0>(curr.specific))
    {

    }
    /// @brief Расчет шага
    /// @param dt Заданный период времени
    /// @param u_in Левое граничное условие
    /// @param u_out Правое граничное условие
    void step(double dt, double u_in, double u_out) {
        auto& F = curr_spec.point_double[0]; // потоки на границах ячеек
        const auto& U = prev_vars.cell_double[0];
        auto& U_new = curr_vars.cell_double[0];

        // Расчет потоков на границе на основе граничных условий
        double v_in = pde.getEquationsCoeffs(0, U[0]);
        double v_out = pde.getEquationsCoeffs(F.size() - 1, U[U.size() - 1]);
        if (v_in >= 0) {
            F.front() = v_in * u_in;
        }
        if (v_out <= 0) {
            F.back() = v_out * u_out;
        }


        // Расчет потоков на границе по правилу донорской ячейки
        for (size_t cell = 0; cell < U.size(); ++cell) {
            double u = U[cell];
            double v_cell = pde.getEquationsCoeffs(cell, u); // не совсем корректно, скорость в ячейке берется из скорости на ее левой границе
            if (v_cell > 0) {
                size_t right_point = cell + 1;
                double v_right = pde.getEquationsCoeffs(right_point, u);
                F[right_point] = u * v_right;
            }
            else {
                size_t left_point = cell;
                double v_left = pde.getEquationsCoeffs(left_point, u);
                F[left_point] = u * v_left;
            }
        }

        for (size_t cell = 0; cell < U.size(); ++cell) {
            double dx = grid[cell + 1] - grid[cell]; // ячейки обычно одинаковой длины, но мало ли..
            U_new[cell] = U[cell] + dt / dx * ((F[cell] - F[cell + 1]));
        }
    }
};

inline double quick_border_approximation(double U_L, double U_C, double U_R)
{
    double Ub_linear = (U_C + U_R) / 2;
    double Ub_cubic_correction = -(U_L + U_R - 2 * U_C) / 8;
    return Ub_linear + Ub_cubic_correction;
}

inline double quickest_border_approximation(double U_L, double U_C, double U_R, double hi, double dx, double dt, double v)
{
    double Cour = abs((v * dt) / dx);
    double Ub_linear = (U_C + U_R) / 2;
    double Grad = (U_R - U_C) / dx;
    double Curv = (U_L + U_R - 2 * U_C) / (dx * dx);
    double Ub_correction_first_order = -(dx * Cour * Grad) / 2;
    double Ub_correction_second_order = (dx * dx * (hi - (1 - (Cour * Cour)) / 3) * Curv) / 2;
    return Ub_linear + Ub_correction_first_order + Ub_correction_second_order;
}

inline double quickest_ultimate_border_approximation(double U_L, double U_C, double U_R, double hi, double dx, double dt, double v)
{
    double DEL = U_R - U_L;
    double ADEL = abs(DEL);
    double ACURV = abs(U_L + U_R - 2 * U_C);
    if (ACURV >= ADEL) {
        return U_C;
    }
    double Cour = abs((v * dt) / dx);
    double REF = U_L + ((U_C - U_L) / Cour);
    double Ub_linear = (U_C + U_R) / 2;
    double Grad = (U_R - U_C) / dx;
    double Curv = (U_L + U_R - 2 * U_C) / (dx * dx);
    double Ub_correction_first_order = -(dx * Cour * Grad) / 2;
    double Ub_correction_second_order = (dx * dx * (hi - (1 - (Cour * Cour)) / 3) * Curv) / 2;
    double Uf = Ub_linear + Ub_correction_first_order + Ub_correction_second_order;
    if (DEL > 0) {
        if (Uf < U_C) {
            Uf = U_C;
        }
        if (Uf > std::min(REF, U_R)) {
            Uf = std::min(REF, U_R);
        }
    }
    else {
        if (Uf > U_C) {
            Uf = U_C;
        }
        if (Uf < std::max(REF, U_R)) {
            Uf = std::max(REF, U_R);
        }
    }
    return Uf;
}

/// @brief Солвер на основе QUICK, только для размерности 1!
/// [Leonard 1979]
class quick_fv_solver {
public:
    typedef typename quick_fv_solver_traits<1>::var_layer_data var_layer_data;
    typedef typename quick_fv_solver_traits<1>::specific_layer specific_layer;
    typedef typename fixed_system_types<1>::matrix_type matrix_type;
    typedef typename fixed_system_types<1>::var_type vector_type;
protected:
    /// @brief ДУЧП
    pde_t<1>& pde;
    /// @brief Сетка, полученная от ДУЧП
    const vector<double>& grid;
    /// @brief Количество точек сетки
    const size_t n;
    /// @brief Предыдущий слой переменных
    const var_layer_data& prev_vars;
    /// @brief Новый (рассчитываемый) слой переменных
    var_layer_data& curr_vars;
    /// @brief Предыдущий специфический слой (сейчас не нужен! нужен ли в будущем?)
    const specific_layer& prev_spec;
    /// @brief Текущий специфический слой
    specific_layer& curr_spec;
public:
    /// @brief Конструктор для буфера в котором простой слой:
    /// когда в слое только один блок целевых переменных и один блок служебных данных
    /// Из буфера берется current() и previous()
    /// @param pde ДУЧП
    /// @param buffer Буфер слоев
    quick_fv_solver(pde_t<1>& pde,
        ring_buffer_t<composite_layer_t<var_layer_data, specific_layer>>& buffer)
        : quick_fv_solver(pde, buffer.previous(), buffer.current())
    {}

    /// @brief Конструктор для простых слоев - 
    /// когда в слое только один блок целевых переменных и один блок служебных данных
    /// @param pde ДУЧП
    /// @param prev Предыдущий слой (уже рассчитанный)
    /// @param curr Следующий (новый), для которого требуется сделать расчет
    quick_fv_solver(pde_t<1>& pde,
        const composite_layer_t<var_layer_data, specific_layer>& prev,
        composite_layer_t<var_layer_data, specific_layer>& curr)
        : pde(pde)
        , grid(pde.get_grid())
        , n(pde.get_grid().size())
        , prev_vars(prev.vars)
        , curr_vars(curr.vars)
        , prev_spec(std::get<0>(prev.specific))
        , curr_spec(std::get<0>(curr.specific))
    {

    }
    /// @brief Расчет шага
    /// @param dt Заданный период времени
    /// @param u_in Левое граничное условие
    /// @param u_out Правое граничное условие
    void step(double dt, double u_in, double u_out) {
        auto& F = curr_spec.point_double[0]; // потоки на границах ячеек
        const auto& U = prev_vars.cell_double[0];
        auto& U_new = curr_vars.cell_double[0];

        // Расчет потоков на границе на основе граничных условий
        double v_in = pde.getEquationsCoeffs(0, U[0]);
        double v_out = pde.getEquationsCoeffs(F.size() - 1, U[U.size() - 1]);
        if (v_in >= 0) {
            F.front() = v_in * u_in;
        }
        if (v_out <= 0) {
            F.back() = v_out * u_out;
        }

        double v_pipe = pde.getEquationsCoeffs(0, U[0]);//не совсем корректно, скорость в ячейке берется из скорости на ее левой границе
        // Расчет потоков на границе по правилу QUICK
        if (v_pipe >= 0) {
            for (size_t cell = 0; cell < U.size(); ++cell) {
                size_t right_border = cell + 1;
                double Vb = v_pipe; // предположили, что скорость на границе во всех точках трубы одна и та же
                double Ub;
                if (cell == 0) {
                    Ub = quick_border_approximation(U[cell], U[cell], U[cell + 1]); // костыль U_L = U_C
                }
                else if (cell == U.size() - 1) {
                    Ub = quick_border_approximation(U[cell - 1], U[cell], U[cell]); // костыль U_R = U_C
                }
                else {
                    Ub = quick_border_approximation(U[cell - 1], U[cell], U[cell + 1]); // честный расчет
                }
                F[right_border] = Ub * Vb;
            }
        }
        else {
            for (size_t cell = 0; cell < U.size(); ++cell) {
                size_t left_border = cell;
                double Vb = v_pipe; // предположили, что скорость на границе во всех точках трубы одна и та же
                double Ub;
                if (cell == 0) {
                    Ub = quick_border_approximation(U[cell + 1], U[cell], U[cell]); // костыль U_L = U_C
                }
                else if (cell == U.size() - 1) {
                    Ub = quick_border_approximation(U[cell], U[cell], U[cell - 1]); // костыль U_R = U_C
                }
                else {
                    Ub = quick_border_approximation(U[cell + 1], U[cell], U[cell - 1]); // честный расчет
                }
                F[left_border] = Ub * Vb;
            }

        }

        for (size_t cell = 0; cell < U.size(); ++cell) {
            double dx = grid[cell + 1] - grid[cell]; // ячейки обычно одинаковой длины, но мало ли..
            U_new[cell] = U[cell] + dt / dx * ((F[cell] - F[cell + 1]));
        }

    }
};

/// @brief Солвер на основе QUICKEST, только для размерности 1!
/// [Neumann 2011]
class quickest_fv_solver {
public:
    typedef typename quickest_fv_solver_traits<1>::var_layer_data var_layer_data;
    typedef typename quickest_fv_solver_traits<1>::specific_layer specific_layer;
    typedef typename fixed_system_types<1>::matrix_type matrix_type;
    typedef typename fixed_system_types<1>::var_type vector_type;
protected:
    /// @brief ДУЧП
    pde_t<1>& pde;
    /// @brief Сетка, полученная от ДУЧП
    const vector<double>& grid;
    /// @brief Количество точек сетки
    const size_t n;
    /// @brief Предыдущий слой переменных
    const var_layer_data& prev_vars;
    /// @brief Новый (рассчитываемый) слой переменных
    var_layer_data& curr_vars;
    /// @brief Предыдущий специфический слой (сейчас не нужен! нужен ли в будущем?)
    const specific_layer& prev_spec;
    /// @brief Текущий специфический слой
    specific_layer& curr_spec;
public:
    /// @brief Конструктор для буфера в котором простой слой:
    /// когда в слое только один блок целевых переменных и один блок служебных данных
    /// Из буфера берется current() и previous()
    /// @param pde ДУЧП
    /// @param buffer Буфер слоев
    quickest_fv_solver(pde_t<1>& pde,
        ring_buffer_t<composite_layer_t<var_layer_data, specific_layer>>& buffer)
        : quickest_fv_solver(pde, buffer.previous(), buffer.current())
    {}

    /// @brief Конструктор для простых слоев - 
    /// когда в слое только один блок целевых переменных и один блок служебных данных
    /// @param pde ДУЧП
    /// @param prev Предыдущий слой (уже рассчитанный)
    /// @param curr Следующий (новый), для которого требуется сделать расчет
    quickest_fv_solver(pde_t<1>& pde,
        const composite_layer_t<var_layer_data, specific_layer>& prev,
        composite_layer_t<var_layer_data, specific_layer>& curr)
        : pde(pde)
        , grid(pde.get_grid())
        , n(pde.get_grid().size())
        , prev_vars(prev.vars)
        , curr_vars(curr.vars)
        , prev_spec(std::get<0>(prev.specific))
        , curr_spec(std::get<0>(curr.specific))
    {

    }
    /// @brief Расчет шага
    /// @param dt Заданный период времени
    /// @param u_in Левое граничное условие
    /// @param u_out Правое граничное условие
    void step(double dt, double u_in, double u_out) {
        auto& F = curr_spec.point_double[0]; // потоки на границах ячеек
        const auto& U = prev_vars.cell_double[0];
        auto& U_new = curr_vars.cell_double[0];

        // Расчет потоков на границе на основе граничных условий
        double v_in = pde.getEquationsCoeffs(0, U[0]);
        double v_out = pde.getEquationsCoeffs(F.size() - 1, U[U.size() - 1]);
        if (v_in >= 0) {
            F.front() = v_in * u_in;
        }
        if (v_out <= 0) {
            F.back() = v_out * u_out;
        }

        double v_pipe = pde.getEquationsCoeffs(0, U[0]);//не совсем корректно, скорость в ячейке берется из скорости на ее левой границе
        // Расчет потоков на границе по правилу QUICK
        if (v_pipe >= 0) {
            for (size_t cell = 0; cell < U.size(); ++cell) {
                size_t right_border = cell + 1;
                double Vb = v_pipe; // предположили, что скорость на границе во всех точках трубы одна и та же
                double Ub;
                if (cell == 0) {
                    Ub = quickest_border_approximation(U[cell], U[cell], U[cell + 1], 0, grid[cell + 1] - grid[cell], dt, v_pipe); // костыль U_L = U_C
                }
                else if (cell == U.size() - 1) {
                    Ub = quickest_border_approximation(U[cell - 1], U[cell], U[cell], 0, grid[cell + 1] - grid[cell], dt, v_pipe); // костыль U_R = U_C
                }
                else {
                    Ub = quickest_border_approximation(U[cell - 1], U[cell], U[cell + 1], 0, grid[cell + 1] - grid[cell], dt, v_pipe); // честный расчет
                }
                F[right_border] = Ub * Vb;
            }
        }
        else {
            for (size_t cell = 0; cell < U.size(); ++cell) {
                size_t left_border = cell;
                double Vb = v_pipe; // предположили, что скорость на границе во всех точках трубы одна и та же
                double Ub;
                if (cell == 0) {
                    Ub = quickest_border_approximation(U[cell + 1], U[cell], U[cell], 0, grid[cell + 1] - grid[cell], dt, v_pipe); // костыль U_L = U_C
                }
                else if (cell == U.size() - 1) {
                    Ub = quickest_border_approximation(U[cell], U[cell], U[cell - 1], 0, grid[cell + 1] - grid[cell], dt, v_pipe); // костыль U_R = U_C
                }
                else {
                    Ub = quickest_border_approximation(U[cell + 1], U[cell], U[cell - 1], 0, grid[cell + 1] - grid[cell], dt, v_pipe); // честный расчет
                }
                F[left_border] = Ub * Vb;
            }

        }

        for (size_t cell = 0; cell < U.size(); ++cell) {
            double dx = grid[cell + 1] - grid[cell]; // ячейки обычно одинаковой длины, но мало ли..
            U_new[cell] = U[cell] + dt / dx * ((F[cell] - F[cell + 1]));
        }

    }
};

template <size_t Dimension>
struct quickest_ultimate_fv_wrapper;

/// @brief Обертка над составным слоем
template <>
struct quickest_ultimate_fv_wrapper<1> {
    typedef typename quickest_ultimate_fv_solver_traits<1>::specific_layer specific_layer;
    
    /// @brief Значения рассчитываемых параметров
    std::vector<double>& vars;
    /// @brief Специфический слой
    specific_layer& specific;

    quickest_ultimate_fv_wrapper(
        vector<double>& U,
        specific_layer& specific
    )
        : vars(U)
        , specific(specific)
    {}
};

/// @brief Солвер на основе QUICKEST-ULTIMATE, только для размерности 1!
/// [Leonard 1991]
class quickest_ultimate_fv_solver {
public:
    typedef typename quickest_ultimate_fv_solver_traits<1>::var_layer_data var_layer_data;
    typedef typename quickest_ultimate_fv_solver_traits<1>::specific_layer specific_layer;
    typedef typename fixed_system_types<1>::matrix_type matrix_type;
    typedef typename fixed_system_types<1>::var_type vector_type;
protected:
    /// @brief ДУЧП
    pde_t<1>& pde;
    /// @brief Сетка, полученная от ДУЧП
    const vector<double>& grid;
    /// @brief Количество точек сетки
    const size_t n;
    /// @brief Предыдущий слой переменных
    const vector<double>& prev_vars;
    /// @brief Новый (рассчитываемый) слой переменных
    vector<double>& curr_vars;
    /// @brief Предыдущий специфический слой (сейчас не нужен! нужен ли в будущем?)
    const specific_layer& prev_spec;
    /// @brief Текущий специфический слой
    specific_layer& curr_spec;
public:
    /// @brief Конструктор для буфера в котором простой слой:
    /// когда в слое только один блок целевых переменных и один блок служебных данных
    /// Из буфера берется current() и previous()
    /// @param pde ДУЧП
    /// @param buffer Буфер слоев
    quickest_ultimate_fv_solver(pde_t<1>& pde,
        ring_buffer_t<composite_layer_t<var_layer_data, specific_layer>>& buffer)
        : quickest_ultimate_fv_solver(pde, buffer.previous(), buffer.current())
    {}

    /// @brief Конструктор для простых слоев - 
    /// когда в слое только один блок целевых переменных и один блок служебных данных
    /// @param pde ДУЧП
    /// @param prev Предыдущий слой (уже рассчитанный)
    /// @param curr Следующий (новый), для которого требуется сделать расчет
    quickest_ultimate_fv_solver(pde_t<1>& pde,
        const composite_layer_t<var_layer_data, specific_layer>& prev,
        composite_layer_t<var_layer_data, specific_layer>& curr)
        : pde(pde)
        , grid(pde.get_grid())
        , n(pde.get_grid().size())
        , prev_vars(prev.vars.cell_double[0])
        , curr_vars(curr.vars.cell_double[0])
        , prev_spec(std::get<0>(prev.specific))
        , curr_spec(std::get<0>(curr.specific))
    {

    }
    /// @brief Конструктор на основе буфера оберток
    /// (созданного с помощью ring_buffer_t::get_custom_buffer)
    /// @param pde ДУЧП
    /// @param wrapper Буфер оберток
    quickest_ultimate_fv_solver(pde_t<1>& pde,
        ring_buffer_t<quickest_ultimate_fv_wrapper<1>>& wrapper)
        : pde(pde)
        , grid(pde.get_grid())
        , n(pde.get_grid().size())
        , prev_vars(wrapper.previous().vars)
        , curr_vars(wrapper.current().vars)
        , prev_spec(wrapper.previous().specific)
        , curr_spec(wrapper.current().specific)
    {}
    /// @brief Конструктор, заточенный для удобства выдергивания специфического слоя, если он один в буфере
    /// Очень специфический
    quickest_ultimate_fv_solver(pde_t<1>& pde,
        const vector<double>& prev_vars, vector<double>& curr_vars,
        const specific_layer& prev_spec, specific_layer& curr_spec)
        : pde(pde)
        , grid(pde.get_grid())
        , n(pde.get_grid().size())
        , prev_vars(prev_vars)
        , curr_vars(curr_vars)
        , prev_spec(prev_spec)
        , curr_spec(curr_spec)
    {

    }

    /// @brief Расчет шага
    /// @param dt Заданный период времени
    /// @param u_in Левое граничное условие
    /// @param u_out Правое граничное условие
    void step(double dt, double u_in, double u_out) {
        auto& F = curr_spec.point_double[0]; // потоки на границах ячеек
        const auto& U = prev_vars;
        auto& U_new = curr_vars;

        // Расчет потоков на границе на основе граничных условий
        double v_in = pde.getEquationsCoeffs(0, U[0]);
        double v_out = pde.getEquationsCoeffs(F.size() - 1, U[U.size() - 1]);
        if (v_in >= 0) {
            F.front() = v_in * u_in;
        }
        if (v_out <= 0) {
            F.back() = v_out * u_out;
        }

        double v_pipe = pde.getEquationsCoeffs(0, U[0]);//не совсем корректно, скорость в ячейке берется из скорости на ее левой границе
        // Расчет потоков на границе по правилу QUICK
        if (v_pipe >= 0) {
            for (size_t cell = 0; cell < U.size(); ++cell) {
                size_t right_border = cell + 1;
                double Vb = v_pipe; // предположили, что скорость на границе во всех точках трубы одна и та же
                double Ub;
                if (cell == 0) {
                    Ub = quickest_ultimate_border_approximation(U[cell], U[cell], U[cell + 1], 0, grid[cell + 1] - grid[cell], dt, v_pipe); // костыль U_L = U_C
                }
                else if (cell == U.size() - 1) {
                    Ub = quickest_ultimate_border_approximation(U[cell - 1], U[cell], U[cell], 0, grid[cell + 1] - grid[cell], dt, v_pipe); // костыль U_R = U_C
                }
                else {
                    Ub = quickest_ultimate_border_approximation(U[cell - 1], U[cell], U[cell + 1], 0, grid[cell + 1] - grid[cell], dt, v_pipe); // честный расчет
                }
                F[right_border] = Ub * Vb;
            }
        }
        else {
            for (size_t cell = 0; cell < U.size(); ++cell) {
                size_t left_border = cell;
                double Vb = v_pipe; // предположили, что скорость на границе во всех точках трубы одна и та же
                double Ub;
                if (cell == 0) {
                    Ub = quickest_ultimate_border_approximation(U[cell + 1], U[cell], U[cell], 0, grid[cell + 1] - grid[cell], dt, v_pipe); // костыль U_L = U_C
                }
                else if (cell == U.size() - 1) {
                    Ub = quickest_ultimate_border_approximation(U[cell], U[cell], U[cell - 1], 0, grid[cell + 1] - grid[cell], dt, v_pipe); // костыль U_R = U_C
                }
                else {
                    Ub = quickest_ultimate_border_approximation(U[cell + 1], U[cell], U[cell - 1], 0, grid[cell + 1] - grid[cell], dt, v_pipe); // честный расчет
                }
                F[left_border] = Ub * Vb;
            }

        }

        for (size_t cell = 0; cell < U.size(); ++cell) {
            double dx = grid[cell + 1] - grid[cell]; // ячейки обычно одинаковой длины, но мало ли..
            double Cr = v_in * dt / dx;
            if (Cr > 1) {
                throw std::runtime_error("Quickest-ultimate is called with Cr > 1");
            }
            U_new[cell] = U[cell] + dt / dx * ((F[cell] - F[cell + 1]));
        }

    }
};

}