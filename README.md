# Partial Differential Equations Solvers

Библиотека предназначена для реализации составных моделей расчета гидравлических задач с учетом движения партий жидкости по трубопроводу

В readme приведены примеры работы с некоторыми сущностями библиотеки

## Оглавление

* [Реализация метода Эйлера для решения дифференциальных уравнений](#Реализация-метода-Эйлера-для-решения-дифференциальных-уравнений)
  * [Пример использования solve_euler ](#Пример-использования-solve_euler )

------------

## Реализация метода Эйлера для решения дифференциальных уравнений

`solve_euler` это функция, которая решает математическую задачу – находит решение заданной системы 
обыкновенных дифференциальных уравнений (ОДУ) методом Эйлера.

В `pde_solvers` система ОДУ представлена базовым классом `ode_t`. Система уравнений 
в частных производных представлена базовым классом `pde_t`. 

На уровне идеи эти сущности представляют собой матричную запись системы дифференциальных уравнений. 
Эти сущности предоставляют методы, позволяющие численно найти решение ДУ.

Численное решение ОДУ методом Эйлера предполагает:
* расчет производной f'(x<sub>i</sub> ) в i-ом узле расчетной сетки
* оценку приращения функции по значению производной ∆y = f'(x<sub>i</sub>)∙(x<sub>i+1</sub>-x<sub>i</sub>)
* расчет значения функции в следующем узле сетки на основании оценки приращения функции: y<sub>i+1</sub> = y<sub>i</sub>+∆y

Соответственно, методы `ode_t` должны позволять рассчитывать значение производной 
в данной точке расчетной сетки или при данных значениях аргументов.

Методы базовых классов ode_t имеют одинаковый набор параметров:
`имя_метода(grid_index, point_vector)`
* где `grid_index` - индекс текущего узла расчетной сетки. Его можно использовать 
для дифференцирования параметра, профиль которого на данной расчетной сетке известен 
(например, гидростатический перепад dz/dt при известном профиле трубы);
* `point_vector` – вектор значений неизвестных в данной точке.

`solve_euler` работает с сущностями типа ode_t и ничего не знает о физическом смысле 
решаемой задачи. Сведение определенной задачи трубопроводного транспорта 
к базовому классу ode_t или pde_t – ответственность разработчика, использующего `pde_solvers`. 

Это сведение выполняется созданием класса, производного от базовых классов `ode_t` или `pde_t`. 
Производные классы переопределяют поведение функций базового класса  на основе математической модели задачи.


### Пример использования solve_euler

Используем solve_euler для расчет профиля давления в классической задаче PQ\
(см. [Лурье 2012](https://elib.gubkin.ru/content/19749 "Пособие в электронной нефтегазовой библиотеке"), раздел 4.2, Задача 1)

**Решаемое диффеенциальное уравнение реализовано следующим образом:**
```C++
/// @brief Уравнение трубы для задачи PQ
class Pipe_model_for_PQ_t : public ode_t<1>
{
public:
    using ode_t<1>::equation_coeffs_type;
    using ode_t<1>::right_party_type;
    using ode_t<1>::var_type;
protected:
    pipe_properties_t& pipe;
    oil_parameters_t& oil;
    double flow;

public:
    /// @brief Констуктор уравнения трубы
    /// @param pipe Ссылка на сущность трубы
    /// @param oil Ссылка на сущность нефти
    /// @param flow Объемный расход
    Pipe_model_for_PQ_t(pipe_properties_t& pipe, oil_parameters_t& oil, double flow)
        : pipe(pipe)
        , oil(oil)
        , flow(flow)
    {

    }

    /// @brief Возвращает известную уравнению сетку
    virtual const vector<double>& get_grid() const override {
        return pipe.profile.coordinates;
    }

    /// @brief Возвращает значение правой части ДУ
    /// @param grid_index Обсчитываемый индекс расчетной сетки
    /// @param point_vector Начальные условия
    /// @return Значение правой части ДУ в точке point_vector
    virtual right_party_type ode_right_party(
        size_t grid_index, const var_type& point_vector) const override
    {
        double rho = oil.density();
        double S_0 = pipe.wall.getArea();
        double v = flow / (S_0);
        double Re = v * pipe.wall.diameter / oil.viscosity();
        double lambda = pipe.resistance_function(Re, pipe.wall.relativeRoughness());
        double tau_w = lambda / 8 * rho * v * abs(v);
        /// Обработка индекса в случае расчетов на границах трубы
        /// Чтобы не выйти за массив высот, будем считать dz/dx в соседней точке
        grid_index = grid_index == 0 ? grid_index + 1 : grid_index;
        grid_index = grid_index == pipe.profile.heights.size() - 1 ? grid_index - 1 : grid_index;

        double height_derivative = (pipe.profile.heights[grid_index] - pipe.profile.heights[grid_index - 1]) /
            (pipe.profile.coordinates[grid_index] - pipe.profile.coordinates[grid_index - 1]);

        return { ((-4) / pipe.wall.diameter) * tau_w - rho * M_G * height_derivative };
    }

};
```
**Применение solve_euler к разработаному классу `Pipe_model_for_PQ_t`**

```C++
int main()
{
    // Создаем сущность трубы
    pipe_properties_t pipe;
    
    // Создаем сущность нефти
    oil_parameters_t oil;

    // Задаем объемнй расход нефти, [м3/с]
    double Q = 0.8;

    // Создаем расчетную модель трубы
    Pipe_model_for_PQ_t pipeModel(pipe, oil, Q);

    // Получаем указатель на начало слоя в буфере
    profile_wrapper<double, 1> start_layer(buffer.current());

    // Задаем конечное давление
    double Pout = 5e5;

    // Вектор для хранения профиля давления
    vector<double> profile(0, pipe.profile.getPointCount());

    /// Модифицированный метод Эйлера для модели pipeModel,
    /// расчет ведется справа-налево относительно сетки,
    /// начальное условие Pout, 
    /// результаты расчета запишутся в слой, на который указывает start_layer
    solve_euler_corrector<1>(pipeModel, -1,  Pout , &profile);
}
```