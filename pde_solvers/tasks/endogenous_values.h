namespace pde_solvers
{
;


/// @brief Эндогенный параметр и его код достоверности
struct endogenous_confident_value_t {
    /// @brief Значение эндогенного свойства
    double value{ std::numeric_limits<double>::quiet_NaN() };
    /// @brief Код достоверности - по умолчанию значение НЕдостоверное
    bool confidence{ false };
};

/// @brief Возвращает значения по умолчанию для эндогенных свойств или их метаданных
/// T = bool - возвращает false (параметр по умолчанию не рассчитывается)
/// T = double - возвращает NaN (значение еще не было расчитано)
/// T = endogenious_confident_value_t - возвращает значения по умолчанию для расчетного знчания + код достоверности
template <typename T>
inline constexpr T default_endogenious_parameter() {
    if constexpr (std::is_same<T, bool>::value) {
        return false;
    }
    else if constexpr (std::is_same<T, double>::value) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    else if constexpr (std::is_same<T, endogenous_confident_value_t>::value) {
        return endogenous_confident_value_t();
    }
    else {
        throw std::runtime_error("Unknown endogenious parameter type");
    }
}

/// @brief Шаблон структуры эндогенных параметров и их метаданных.
/// Например, в зависимости от типа T возможны следующие сценарии использования:
/// T = double - расчетное значение эндогенного свойства
/// T = bool - настройка - выполнять ли расчет свойства
/// T = endogenious_confident_value_t - расчетное значение свойства и его достоверность
template <typename T>
struct endogenous_parameters_template_t {
    /// @brief Номинальная плотность
    T density_std{ default_endogenious_parameter<T>() };
    /// @brief Номинальная вязкость для изотермического случая
    // TODO: не ли здесь дублирования с visc20?
    T viscosity_working{ default_endogenious_parameter<T>() };
    /// @brief Серосодержанние
    T sulfur{ default_endogenious_parameter<T>() };
    /// @brief Концентрация ПТП
    T improver{ default_endogenious_parameter<T>() };
    /// @brief Температура
    T temperature{ default_endogenious_parameter<T>() };
    /// @brief Вязкость при 0С
    T viscosity0{ default_endogenious_parameter<T>() };
    /// @brief Вязкость при 20С
    T viscosity20{ default_endogenious_parameter<T>() };
    /// @brief Вязкость при 50С
    T viscosity50{ default_endogenious_parameter<T>() };
};

/// @brief Селектор эндогенных параметров - используется в моделях и в смесителе как настройка расчета
using endogenous_selector_t = endogenous_parameters_template_t<bool>;
/// @brief Структура расчетных значений эндогенных свойств для распространения по объектам графа. 
/// Значаения эндогенных свойств определяются либо из измерения, либо в результате расчета.
/// Значения полей могут быть в следуюищх состояниях:
/// 1. Не задано из измерения и еще не расчитано
/// 2. Из расчета получено неадевтаное NaN значение
/// 3. Из расчета получено численное значение
using endogenous_double_values_t = endogenous_parameters_template_t<double>;
/// @brief Эндогенные параметры и их достоверность
using endogenous_values_t = endogenous_parameters_template_t<endogenous_confident_value_t>;


/// @brief Сбрасывает код достоверности в "ложь"
inline void reset_confidence(endogenous_values_t* values)
{
    auto invalidate = [](endogenous_confident_value_t& confidence_value) {
        confidence_value.confidence = false;
        };

    invalidate(values->density_std);
    invalidate(values->viscosity_working);
    invalidate(values->sulfur);
    invalidate(values->improver);
    invalidate(values->temperature);
    invalidate(values->viscosity0);
    invalidate(values->viscosity20);
    invalidate(values->viscosity50);
}


/// @brief Интерпретирует степень достоверность как булев флаг достоверности
inline bool discriminate_confidence_level(double confidence_level)
{
    return confidence_level > 0.95; // с константой экспериментируем
}

/// @brief Расчетный слой и его код достоверности
struct confident_layer_t {
    /// @brief Сам расчетный параметр
    std::vector<double> value;
    /// @brief Код достоверности - считается численным методом, 
    /// поэтому не булевый а вещественный
    std::vector<double> confidence;
    /// @brief Принимает количество точек, инициализирует количество ячеек
    confident_layer_t(size_t point_count, double initial_value)
        : value(point_count - 1, initial_value)
        , confidence(point_count - 1, 0.0)
    { }
    bool is_confident_layer() const {
        for (double confidence_value : confidence) {
            if (discriminate_confidence_level(confidence_value) == false)
                return false;
        }
        return true;
    }
};

/// @brief Расчетные целевые профили по трубе 
struct pipe_endogenous_variable_layer_t
{
    /// @brief Номинальный объемный расход
    double std_volumetric_flow{ std::numeric_limits<double>::quiet_NaN() };
    /// @brief Массовый расход
    std::vector<double> mass_flow;
    /// @brief Скорость потока
    std::vector<double> velocity;
    /// @brief Профиль плотности
    confident_layer_t density_std;
    /// @brief Профиль вязкости (рабочая при изотермическом расчете)
    confident_layer_t viscosity_working;
    /// @brief Серосодержанние
    confident_layer_t sulfur;
    /// @brief Концентрация ПТП
    confident_layer_t improver;
    /// @brief Температура
    confident_layer_t temperature;
    /// @brief Вязкость при 0С (сортовая при неизотермическом расчете)
    confident_layer_t viscosity0;
    /// @brief Вязкость при 20С (сортовая при неизотермическом расчете)
    confident_layer_t viscosity20;
    /// @brief Вязкость при 50С (сортовая при неизотермическом расчете)
    confident_layer_t viscosity50;
    /// @brief Принимает количество точек, инициализирует количество ячеек
    /// ВАЖНО!!! 
    /// векторы mass_flow, velocity ничего не знают про точки ячейки
    /// а confident_layer_t знает про это и ожидает в конструкторе ТОЧКИ!!!
    pipe_endogenous_variable_layer_t(size_t point_count)
        : mass_flow(point_count - 1, 0.0)
        , velocity(point_count - 1, 0.0)
        , density_std(point_count, 860)
        , viscosity_working(point_count, 1e-6)
        , sulfur(point_count, 1e-3)
        , improver(point_count, 0.0)
        , temperature(point_count, 300)
        , viscosity0(point_count, 0.0)
        , viscosity20(point_count, 0.20e-6)
        , viscosity50(point_count, 0.50e-6)
    {

    }
};

/// @brief Сбрасывает код достоверности в "ложь"
inline void reset_confidence(pipe_endogenous_variable_layer_t* layer) {
    auto invalidate = [](confident_layer_t& parameter_layer) {
        std::fill(parameter_layer.confidence.begin(), parameter_layer.confidence.end(), false);
        };

    invalidate(layer->density_std);
    invalidate(layer->viscosity_working);
    invalidate(layer->sulfur);
    invalidate(layer->improver);
    invalidate(layer->temperature);
    invalidate(layer->viscosity0);
    invalidate(layer->viscosity20);
    invalidate(layer->viscosity50);
}

/// @brief Расчетный профиль на квикест
struct pipe_endogenous_calc_layer_t : public pipe_endogenous_variable_layer_t
{
    /// @brief Профиль вспомогательных расчетов для метода конечных объемов (и для вязкости, и для плотности)
    quickest_ultimate_fv_solver_traits<1>::specific_layer specific;
    /// @brief Инициализация профилей
    /// @param point_count Количество точек
    pipe_endogenous_calc_layer_t(size_t point_count)
        : pipe_endogenous_variable_layer_t(point_count)
        , specific(point_count)
    {
    }

    /// @brief Подготовка плотности для расчета методом конечных объемов   
    static quickest_ultimate_fv_wrapper<1> get_density_std_wrapper(
        pipe_endogenous_calc_layer_t& layer)
    {
        return quickest_ultimate_fv_wrapper<1>(layer.density_std.value, layer.specific);
    }
    static quickest_ultimate_fv_wrapper<1> get_density_std_confidence_wrapper(
        pipe_endogenous_calc_layer_t& layer)
    {
        return quickest_ultimate_fv_wrapper<1>(layer.density_std.confidence, layer.specific);
    }
    /// @brief Подготовка вязкости для расчета методом конечных объемов 
    static quickest_ultimate_fv_wrapper<1> get_viscosity_working_wrapper(
        pipe_endogenous_calc_layer_t& layer)
    {
        return quickest_ultimate_fv_wrapper<1>(layer.viscosity_working.value, layer.specific);
    }
    static quickest_ultimate_fv_wrapper<1> get_viscosity_working_confidence_wrapper(
        pipe_endogenous_calc_layer_t& layer)
    {
        return quickest_ultimate_fv_wrapper<1>(layer.viscosity_working.confidence, layer.specific);
    }
};


}

