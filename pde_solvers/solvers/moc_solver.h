#pragma once

namespace pde_solvers {

//using std::array;
using std::pair;

///// @brief Структура-родитель всех task_traits
///// \tparam Dimension
//template <size_t Dimension>
//struct task_traits
//{
//    // TODO добавить var_layer_data, specific_layer
//};


template <size_t Dimension>
struct moc_task_traits;

/// @brief Описание типов данных для метода характеристик заданной размерности
template <>
struct moc_task_traits<1>
{
    typedef profile_collection_t<1/*переменные*/, 0, 0, 0, 0, 0> var_layer_data;
    typedef std::vector<double> specific_layer;
};


/// @brief Описание типов данных для метода характеристик заданной размерности
template <size_t Dimension>
struct moc_task_traits 
{
    typedef profile_collection_t<Dimension/*переменные*/, 0, 0, 0, 0, 0> var_layer_data;
    typedef profile_collection_t<Dimension/*собственные числа*/, 0,
        Dimension /*собств. векторы*/, Dimension /*размерность собств. векторов*/,
        0, 0> specific_layer;
};


template <size_t Dimension>
struct layer_wrapper {
    // TODO matrix_type, vector_type
};

template <size_t Dimension>
struct moc_layer_wrapper;

/// @brief Обертка над составным слоем для метода характеристик
template <>
struct moc_layer_wrapper<1> : layer_wrapper<1> {
    typedef typename moc_task_traits<1>::var_layer_data var_layer_data;
    typedef typename moc_task_traits<1>::specific_layer specific_layer_data;

    /// @brief Значения рассчитываемых параметров
    vector<double>& values;
    /// @brief Собственные числа. 
    vector<double>& eigenval;

    moc_layer_wrapper(
        vector<double>& U,
        vector<double>& eigenval
    )
        : values(U)
        , eigenval(eigenval)
    {}
};

/// @brief Обертка над составным слоем для метода характеристик
template <size_t Dimension>
struct moc_layer_wrapper : layer_wrapper<Dimension> {
    typedef typename fixed_system_types<Dimension>::var_type vector_type;
    typedef typename moc_task_traits<Dimension>::var_layer_data var_layer_data;
    typedef typename moc_task_traits<Dimension>::specific_layer specific_layer_data;

    /// @brief Значения рассчитываемых параметров
    profile_wrapper<double, Dimension> values;
    /// @brief Собственные числа. 
    profile_wrapper<double, Dimension> eigenval;
    /// @brief Собственные векторы (левые). 
    profile_wrapper<vector_type, Dimension> eigenvec;

    moc_layer_wrapper(
        profile_wrapper<double, Dimension> U,
        profile_wrapper<double, Dimension> eigenval,
        profile_wrapper<array<double, Dimension>, Dimension> eigenvec
    )
        : values(U)
        , eigenval(eigenval)
        , eigenvec(eigenvec)
    {}
    moc_layer_wrapper(
        array<vector<double>*, Dimension> U,
        array<vector<double>*, Dimension> eigenval,
        array<vector<array<double, Dimension>>*, Dimension> eigenvec
    )
        : values(U)
        , eigenval(eigenval)
        , eigenvec(eigenvec)
    {}
    moc_layer_wrapper(
        var_layer_data& vars,
        specific_layer_data& specific
    )
        : values(get_profiles_pointers(vars.point_double))
        , eigenval(get_profiles_pointers(specific.point_double))
        , eigenvec(get_profiles_pointers(specific.point_vector))
    {

    }

    moc_layer_wrapper(
        profile_wrapper<double, Dimension> U,
        specific_layer_data& specific
    )
        : values(U)
        , eigenval(get_profiles_pointers(specific.point_double))
        , eigenvec(get_profiles_pointers(specific.point_vector))
    {

    }
    moc_layer_wrapper(
        vector<double>& U,
        specific_layer_data& specific
    )
        : values(U)
        , eigenval(get_profiles_pointers(specific.point_double))
        , eigenvec(get_profiles_pointers(specific.point_vector))
    {

    }

};

/// @brief Расчетчик метода характеристик
/// @tparam Dimension Размерность задачи
template <size_t Dimension>
class moc_solver;


/// @brief Расчетчик метода характеристик
/// @tparam Dimension Размерность задачи
template <>
class moc_solver<1>
{
public:
    typedef typename moc_task_traits<1>::specific_layer specific_layer;
    //typedef typename fixed_system_types<1>::matrix_type matrix_type;
    //typedef typename fixed_system_types<1>::var_type vector_type;

protected:
    /// @brief ДУЧП
    pde_t<1>& pde;
    /// @brief Сетка, полученная от ДУЧП
    const vector<double>& grid;
    /// @brief Количество точек сетки
    const size_t n;
    /// @brief Предыдущий слой (начальные условия)
    vector<double>& prev;
    /// @brief Рассчитываемый слой (он же новый, следующий, текущий)
    vector<double>& curr;
    /// @brief Вспомогательный буфер для расчета собственных чисел 
    vector<double>& eigenvals;

public:
    /// @brief Базовый конструктор, наиболее детальный
    /// @param pde ДУЧП
    /// @param prev Предыдуший слой
    /// @param curr Новый слой
    /// @param eigenvals Буфер для расчета собственных чисел (рекомендуется относить к прошлому слою)
    moc_solver(pde_t<1>& pde,
        vector<double>& prev,
        vector<double>& curr,
        vector<double>& eigenvals
        )
        : pde(pde)
        , grid(pde.get_grid())
        , n(pde.get_grid().size())
        , prev(prev)
        , curr(curr)
        , eigenvals(eigenvals)
    { }
    /// @brief Конструктор на основе слове представленных через MOC-обертку
    moc_solver(pde_t<1>& pde,
        moc_layer_wrapper<1>& prev,
        moc_layer_wrapper<1>& curr)
        : moc_solver(pde, prev.values, curr.values, prev.eigenval)
    { }
    /// @brief Конструктор на основе буфера оберток 
    /// (созданного с помощью ring_buffer_t::get_custom_buffer)
    moc_solver(pde_t<1>& pde,
        ring_buffer_t<moc_layer_wrapper<1>>& buffer)
        : moc_solver(pde, buffer[-1], buffer[0])
    { }
    /// @brief Конструктор, заточенный для удобства выдергивания специфического слоя, если он один в буфере
    /// Очень специфический
    moc_solver(pde_t<1>& pde, vector<double>& prev, vector<double>& curr,
        std::tuple<vector<double>>& eigenvals)
        : moc_solver(pde, prev, curr, std::get<0>(eigenvals))
    { }
    /// @brief Еще один специфический конструктор, когда composite_layer_t содержит только одну задачу
    moc_solver(pde_t<1>& pde,
        composite_layer_t<profile_collection_t<1>, moc_solver<1>::specific_layer>& prev,
        composite_layer_t<profile_collection_t<1>, moc_solver<1>::specific_layer>& curr)
        : moc_solver(pde, prev.vars.point_double[0], curr.vars.point_double[0], prev.specific)
    {

    }


    /// @brief Расчет собственных чисел и на их основе расчет шага
    /// по Куранту dtCr (т.е. dt, при котором Cr=1).
    /// Если желаемый шаг превышает шаг по Куранту dtCr, либо не задан (time_step = nan),
    /// то возвращается шаг dtCr
    /// Иначе - возвращается time_step
    double prepare_step(double time_step = std::numeric_limits<double>::quiet_NaN()) {
        auto& values = prev;

        double max_egenval = 0;
        for (size_t grid_index = 0; grid_index < grid.size(); ++grid_index) {
            double eigen_value = eigenvals[grid_index] = pde.getEquationsCoeffs(grid_index, values[grid_index]);

            max_egenval = std::max(max_egenval, std::abs(eigen_value));
        }

        double dx = grid[1] - grid[0];
        double courant_step = dx / max_egenval;
        if (std::isnan(time_step) || time_step > courant_step) {
            time_step = courant_step;
        }
        return time_step;
    }
    /// @brief Расчет смещения при неточном попадании линий характеристик между точками сетки
    static double characteristic_interpolation_offset(double dt, const double* lambda, const double* grid) 
    {
        //  linear interpolation

        if (lambda[0] < 0) {
            double dl = grid[1] - grid[0];
            double l0 = lambda[0];
            double l1 = lambda[+1];
            double p = dt * l0 / (dt * (l0 - l1) - dl);
            return p;
        }
        else if (lambda[0] > 0) {
            double dl = grid[0] - grid[-1];
            double l0 = lambda[0];
            double l1 = lambda[-1];
            double p = dt * l0 / (dt * (l1 - l0) - dl);
            return p;

        }
        else {
            double p = 0;
            return p;
        }
    }
    /// @brief Расчет внутренних точек методом первого порядка
    /// (с учетом наклона характеристик)
    double step_inner(double time_step = std::numeric_limits<double>::quiet_NaN())
    {
        time_step = prepare_step(time_step);

        int index_from = eigenvals[0] > 0
            ? 1
            : 0;
            
        int index_to = eigenvals[n - 1] < 0
            ? static_cast<int>(grid.size() - 2)
            : static_cast<int>(grid.size() - 1);
            

        auto& curr_values = curr;

        profile_wrapper<double, 1> prev_values(this->prev); // оборачиваем только для интерполяции 

        for (int index = index_from; index <= index_to; ++index)
        {
            double p = characteristic_interpolation_offset(
                time_step, &eigenvals[index], &grid[index]);
            double u_old = prev_values.interpolate(index, p);
            double b = pde.getSourceTerm(index, u_old);

            double& u_new = curr_values[index];
            u_new = u_old - time_step * b;
        }

        return time_step;
    }
    /// @brief Опциональный расчет нового слоя, учитываюшего граничные условия, 
    /// в зависимости от наклона характеристик
    void step_optional_boundaries(
        double time_step, double left_value, double right_value)
    {
        step_inner(time_step);
        if (eigenvals[0] > 0) {
            curr[0] = left_value;
        }
        if (eigenvals[n - 1] < 0) {
            curr[n - 1] = right_value;
        }
    }
    /// @brief Расчет внутренних точек нового слоя методом второго порядка
    /// Учитывает, что конфигураций характеристики могут позволить рассчитать граничные точки
    /// Также считает собственные числа, векторы для ВСЕХ точек
    double step2_inner(double time_step = std::numeric_limits<double>::quiet_NaN())
    {
        time_step = prepare_step(time_step);

        auto prev_values = profile_wrapper<double, 1>(prev);
        //auto& curr_values = curr;

        constexpr double eps = 1e-8;

        double dl = grid[1] - grid[0];

        for (int grid_index = 0; grid_index < grid.size(); ++grid_index)
        {
            const double& eigenval = eigenvals[grid_index];
            if (grid_index == 0 && eigenval > 0) {
                // надо брать точку с координатой i = -1
                continue;
            }
            if (grid_index == grid.size() - 1 && eigenval < 0) {
                continue;
            }

            double p = characteristic_interpolation_offset(time_step, &eigenval, &grid[grid_index]);

            // предиктор
            double u_old;
            double rp1;
            double absp = abs(p);
            if (absp < eps || abs(1.0 - absp) < eps) {
                // характеристика точно между двумя точками: либо косая, либо вертикальная
                size_t index = static_cast<size_t>(grid_index + p + 0.5);
                u_old = prev[index];
                rp1 = pde.getSourceTerm(index, u_old);
            }
            else {
                // интерполяция правой части
                size_t grida = static_cast<size_t>(grid_index + sgn(p));
                size_t gridb = grid_index;
                double rp1a = pde.getSourceTerm(grida, prev_values[grida]);
                double rp1b = pde.getSourceTerm(gridb, prev_values[gridb]);
                rp1 = rp1a * absp + rp1b * (1 - absp); // проверка: если p = 0, то берем b(grid_index)
                u_old = prev_values.interpolate(grid_index, p);
            }

            double u_estimate = u_old + time_step * rp1;

            // корректор
            double rp2 = pde.getSourceTerm(grid_index, u_estimate);
            curr[grid_index] = u_old + time_step * 0.5 * (rp1 + rp2);

            if (!isfinite(rp1) || !isfinite(rp2) || !isfinite(u_estimate) || !isfinite(u_old)) {
                throw std::logic_error("infinite value");
            }

        }

        return time_step;
    }
    /// @brief Опциональный расчет граничных условий, в зависимости от наклона характеристик
    /// Метод второго порядка
    void step2_optional_boundaries(double time_step,
        double left_value, double right_value)
    {
        step2_inner(time_step);

        if (eigenvals[0] > 0) {
            curr[0] = left_value;
        }
        if (eigenvals[n - 1] < 0) {
            curr[n - 1] = right_value;
        }
    }
};

/// @brief Расчетчик метода характеристик
/// @tparam Dimension Размерность задачи
template <size_t Dimension>
class moc_solver {
public:
    typedef typename moc_task_traits<Dimension>::specific_layer specific_layer;
    typedef typename fixed_system_types<Dimension>::matrix_type matrix_type;
    typedef typename fixed_system_types<Dimension>::var_type vector_type;

//protected:
    /// @brief ДУЧП
    pde_t<Dimension>& pde;
    /// @brief Сетка, полученная от ДУЧП
    const vector<double>& grid;
    /// @brief Количество точек сетки
    const size_t n;

    moc_layer_wrapper<Dimension> curr;
    moc_layer_wrapper<Dimension> prev;

protected:
    /// @brief Надо очень подробно задокументировать нотацию и ограничения использования
    /// @param lambda 
    /// @param grid 
    /// @return 
    double characteristic_interpolation_offset(double dt, const double* lambda, const double* grid) const
    {
        //  linear interpolation

        if (lambda[0] < 0) {
            double dl = grid[1] - grid[0];
            double l0 = lambda[0];
            double l1 = lambda[+1];
            double p = dt * l0 / (dt * (l0 - l1) - dl);
            return p;
        }
        else if (lambda[0] > 0) {
            double dl = grid[0] - grid[-1];
            double l0 = lambda[0];
            double l1 = lambda[-1];
            double p = dt * l0 / (dt * (l1 - l0) - dl);
            return p;

        }
        else {
            double p = 0;
            return p;
        }
    }
    double characteristic_interpolation_offset(double dt, const double* lambda, double dl) const
    {
        //  linear interpolation

        if (lambda[0] < 0) {
            double l0 = lambda[0];
            double l1 = lambda[+1];
            double p = dt * l0 / (dt * (l0 - l1) - dl);
            return p;
        }
        else if (lambda[0] > 0) {
            double l0 = lambda[0];
            double l1 = lambda[-1];
            double p = dt * l0 / (dt * (l1 - l0) - dl);
            return p;

        }
        else {
            double p = 0;
            return p;
        }
    }

public:
    /// @brief Конструктор для простых слоев - 
    /// когда в слое только один блок целевых переменных и один блок служебных данных
    /// @param pde Экземпляр уравнения
    /// @param prev Прошлый слой (начальные условия)
    /// @param curr Новый, рассчитываемый слой
    moc_solver(pde_t<Dimension>& pde,
        composite_layer_t<profile_collection_t<Dimension>, moc_solver<Dimension>::specific_layer>& prev,
        composite_layer_t<profile_collection_t<Dimension>, moc_solver<Dimension>::specific_layer>& curr);

    moc_solver(pde_t<Dimension>& pde,
        moc_layer_wrapper<Dimension>& prev,
        moc_layer_wrapper<Dimension>& curr)
        : pde(pde)
        , grid(pde.get_grid())
        , n(pde.get_grid().size())
        , prev(prev)
        , curr(curr)
    {

    }

public:
    /// @brief Формирует уравнение на характеристической линии, заданной индексом eigenval_index
    /// li * u_new = li * (u_old + dt*b) [обозначим si = li * (u_old + dt*b)]
    /// li * u_new = si
    /// Выполняет линейную интерполяцию значений с учетом шага dt и величины собственного числа
    /// @param time_step Величина временного шага
    /// @param eigenval_index 
    /// @param grid_index 
    /// @return 
    std::pair<vector_type, double> get_characteristic_equation(
        double time_step, size_t eigenval_index, size_t grid_index) const
    {
        const profile_wrapper<double, Dimension>& eigenvals = prev.eigenval;
        const profile_wrapper<vector_type, Dimension>& eigenvecs = prev.eigenvec;
        const profile_wrapper<double, Dimension>& values = prev.values;

        const double* eigenval = &eigenvals.profile(eigenval_index)[grid_index];
        double p = characteristic_interpolation_offset(time_step, eigenval, &grid[grid_index]);

        vector_type li = eigenvecs.interpolate_dimension(eigenval_index, grid_index, p);
        vector_type u_old = values.interpolate(grid_index, p);

        // тут не совсем логично, grid_index не учитывает интерполяцию
        // может быть добавить туда смещение, вроде: getSourceTerm(grid_index, p, u_old); 
        vector_type b = pde.getSourceTerm(grid_index, u_old);

        //L[eigenval_index] = li;
        double s = inner_prod(li, u_old + time_step * b);

        return std::make_pair(li, s);

    }

    std::pair<matrix_type, vector_type> get_characteristic_equations(double time_step, size_t grid_index) const
    {
        matrix_type L;
        vector_type S;

        for (size_t eigenval_index = 0; eigenval_index < Dimension; ++eigenval_index)
        {
            std::tie(L[eigenval_index], S[eigenval_index])
                = get_characteristic_equation(time_step, eigenval_index, grid_index);
        }

        return make_pair(L, S);
    }

    /// @brief Расчет всех точек нового слоя (работает только с гидроударом)
    /// \param time_step
    /// \param left_boundary
    /// \param right_boundary
    double step(const pair<vector_type, double>& left_boundary,
        const pair<vector_type, double>& right_boundary,
        double time_step = std::numeric_limits<double>::quiet_NaN())
    {
        time_step = step_inner(time_step); // если отдать в step_inner dt = nan, то он его пересчитает в шаг по Куранту!

        pair<vector_type, double> eq_left =
            get_characteristic_equation(time_step, 0, 0);
        pair<vector_type, double> eq_right =
            get_characteristic_equation(time_step, 1, grid.size() - 1);

        curr.values(0) =
            solve_linear_system({ eq_left.first, left_boundary.first }, { eq_left.second, left_boundary.second });
        curr.values(n - 1) =
            solve_linear_system({ eq_right.first, right_boundary.first }, { eq_right.second, right_boundary.second });
        return time_step;
    }

    /// @brief Опциональный расчет граничных условий, в зависимости от наклона характеристик
    /// Реализация только для размерности 1
    /// @param time_step 
    /// @param left_value 
    /// @param right_value 
    void step_optional_boundaries(double time_step,
        const double& left_value,
        const double& right_value);




    static double get_max_abs(double v)
    {
        return abs(v);
    }

    static double get_max_abs(const std::array<double, Dimension>& v)
    {
        double max_egenval = -std::numeric_limits<double>::infinity();
        for (double eval : v) {
            max_egenval = std::max(max_egenval, abs(eval));
        }
        return max_egenval;
    }

    static double get_max(double v)
    {
        return v;
    }
    static double get_max(const std::array<double, Dimension>& v)
    {
        double max_egenval = -std::numeric_limits<double>::infinity();
        for (double eval : v) {
            max_egenval = std::max(max_egenval, eval);
        }
        return max_egenval;
    }

    double prepare_step(double time_step = std::numeric_limits<double>::quiet_NaN()) {
        auto& eigenval = prev.eigenval;
        auto& eigenvec = prev.eigenvec;
        auto& values = prev.values;

        double max_egenval = 0;
        for (size_t grid_index = 0; grid_index < grid.size(); ++grid_index) {
            auto [val, vec] = pde.GetLeftEigens(grid_index, values(grid_index));
            
            max_egenval = std::max(max_egenval, get_max_abs(val));
            eigenval(grid_index) = val;
            eigenvec(grid_index) = vec;
        }

        double dx = grid[1] - grid[0];
        double courant_step = dx / max_egenval;
        if (std::isnan(time_step) || time_step > courant_step) {
            time_step = courant_step;
        }
        return time_step;
    }


    /// @brief Расчет внутренних точек нового слоя
    /// Также считает собственные числа, векторы для ВСЕХ точек
    /// \param time_step
    double step_inner(double time_step = std::numeric_limits<double>::quiet_NaN())
    {
        time_step = prepare_step(time_step);

        auto& eigenval = prev.eigenval;
        auto& eigenvec = prev.eigenvec;
        auto& values = prev.values;


        int index_from = 1;
        int index_to = static_cast<int>(grid.size() - 2);

        auto& curr_values = curr.values;

        for (int index = index_from; index <= index_to; ++index)
        {
            // li * u_new = li * (u_old - dt*b) [обозначим si = li * (u_old - dt*b)]
            // L * u_new = S
            auto [L, S] = get_characteristic_equations(time_step, index);
            
            curr_values(index) = solve_linear_system(L, S);
            
        }

        return time_step;
    }


};


/// @brief Решатель транспортного уравнения методом характеристик, 
/// при этом считается, что скорость по длине трубопровода постоянна,
/// а число Куранта всегда равно единице
class transport_moc_solver
{
public:
    /// @brief Конструктор транспортного солвера
    /// @param pipe Параметры трубопровода 
    /// @param vol_flow Объёмный расход
    /// @param prev Предыдущий слой
    /// @param next Новый слой
    transport_moc_solver(const pipe_properties_t& pipe, double vol_flow,
        vector<double>& prev, vector<double>& next)
        : pipe{ pipe }, volumetric_flow{ vol_flow },
        prev{ prev }, next{ next }
    {}

    /// @brief Расчёт нового слоя
    /// @param par_in Значение параметра среды, втекающей в начало трубопровода
    /// @param par_out Значение параметра среды, втекающей в конец трубопровода при обратном течении 
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

    /// @brief Расчёт шага по времени, при котором Курант равен единице (Cr = 1)
    double prepare_step() const
    {
        const std::vector<double>& grid = pipe.profile.coordinates;

        double max_eigen = get_eigen_value();
        double dx = grid[1] - grid[0];
        double courant_step = dx / max_eigen;

        return courant_step;
    }

protected:

    /// @brief модель трубопровода
    const pipe_properties_t& pipe;
    /// @brief Объемный расход
    const double volumetric_flow;
    /// @brief Предыдущий слой
    vector<double>& prev;
    /// @brief Новый слой
    vector<double>& next;


    /// @brief Расчёт собственного значения
    double get_eigen_value() const
    {
        double S = pipe.wall.getArea();
        return volumetric_flow / S;
    }
};

}

