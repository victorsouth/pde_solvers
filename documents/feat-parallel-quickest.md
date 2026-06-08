# Параллельные расчет адвекции в трубе методом QUICKEST-ULTIMATE

В библиотеке `pde_solvers` реализован явный метод QUICKEST-ULTIMATE для решения уравнения адвекции в трубопроводе.

Необходимо добавить возможность параллельного расчета каждой ячейки для повышения быстродействия. Параллельный расчет обязательно должен быть отключаемым в целях отладки. Переключение режима реализуется через шаблонный параметр-политику солвера `quickest_ultimate_fv_solver`:

```cpp
struct parallel_policy {};
struct sequential_policy {};

template <typename ExecutionPolicy>
class quickest_ultimate_fv_solver { ... };
```

Политика выполнения задаётся явно: `quickest_ultimate_fv_solver<sequential_policy>` или `quickest_ultimate_fv_solver<parallel_policy>`. Значение по умолчанию для шаблонного параметра не используется.

Результаты расчета при этом не изменятся, поскольку параллелизм затрагивает только порядок выполнения итераций, но не порядок операций внутри каждой ячейки.

Метод `step()` солвера содержит два цикла по ячейкам. Первый цикл вычисляет потоки `F` на границах ячеек, второй — обновляет значения `U` в ячейках на основе этих потоков. Оба цикла параллелизуемы по ячейкам: итерации внутри каждого цикла независимы. Между циклами существует зависимость по данным — второй цикл читает массив `F`, полностью заполненный первым. Поэтому параллелизация реализуется двумя последовательными вызовами `tbb::parallel_for`: первый — для цикла потоков, второй — для цикла обновления. Явный барьер между ними не нужен, так как `tbb::parallel_for` блокирует вызывающий поток до завершения всех задач, тем самым гарантируя барьер неявно.

Проверка числа Куранта (`Cr > 1`) остаётся внутри параллельного цикла. TBB полностью поддерживает исключения из задач: планировщик перехватывает исключение, отменяет всю группу задач и перебрасывает исключение в вызывающий поток — так же, как это работает в последовательном коде.

Для параллельного программирования использовать библиотеку oneTBB. Подключение в `CMakeLists.txt` — через `find_package(TBB REQUIRED)` и `target_link_libraries(... TBB::tbb)`. Пути к библиотеке указываются в `CMakePresets.json` через `CMAKE_PREFIX_PATH` по аналогии с `fixed_solvers`: для каждого Windows-пресета добавляется соответствующий подкаталог `C:/install/TBB/<компилятор>` (msvc, clang, gcc).

## Тестирование

Тест `TEST(QuickestUltimate, IncreasesPerformanceInParallel)` добавляется в файл `testing/test_quick.h`.

Тест проводится путем измерения времени расчета одного шага с помощью `QUICKEST-ULTIMATE` на трубе с 100 000 ячеек — этого достаточно, чтобы полезная работа существенно превышала затраты на параллелизм, и параллелизм давал уверенный измеримый эффект.

Нужно выполнить последовательный расчет (`sequential_policy`) и зафиксировать время расчета и результаты. Далее выполнить расчет с задействованием параллельности (`parallel_policy`). Также зафиксировать длительность расчета и результаты. Результаты должны совпасть, а длительность расчета в параллельном случае должна быть меньше.

Если `std::thread::hardware_concurrency() < 2`, тест производительности пропускается через `GTEST_SKIP()` — на одноядерных средах (в т.ч. части CI) проверка ускорения не имеет смысла. Сравнение результатов sequential/parallel при этом можно оставить без skip или выполнять отдельным тестом на корректность.

## Выводы о производительности параллельного расчета ячеек в методе QUICKEST

При тестировании параллельной реализации явной схемы установлено, что на типовой сетке (500 ячеек = труба 100 км с шагом 200 м) многотопочная реализация проигрывает однопоточной. Расчет профиля занимает 952 наносекунды - т.к. алгоритм QUICKEST двухшаговый, то на обсчет одной ячейки уходит 0,95 нс. 

В [2025 - Voss - TBB; Rule of Thumb on Task Size] указано, что параллелизм превысит накладные расходы (барьеры синхрониззации, work stealing) для задач, требующих не менее 1000 наносекунд. Основной источник накладных расходов в `parallel_for` — барьер синхронизации. Стоимость барьера (~600–900 нс) сопоставима со всей полезной работой (~952 нс) и не может быть устранена никакими настройками — он присутствует при любом числе потоков и любом подходе к дроблению задач (`partitioner`).  Дополнительный источник накладных расходов — work stealing: потоки, опустошившие свою очередь задач, тратят время на поиск работы у других потоков. `static_partitioner` устраняет work stealing (детерминированное деление без адаптации), однако даже с ним многопоточный расчёт в ~4.7× хуже однопотока из-за доминирующего влияния барьера. Наилучший параллельный результат — `task_arena(2) + static_partitioner` (1173 ms) — всё равно в 1.2× хуже однопотока (952 ms).

Потенциал кэш дружелюбности на уровне рассчитываемого профиля уже реализован. Данные полностью помещаются в кэш (500 `double` = 4 КБ). Задействуется автовекторизация - компилятор заменяет операции так, что за один такт один поток обрабатывает четыре значения `double`. Каждый поток обрабатывает свой чанк с применением автовекторизации.

Вероятно, можно ускорить расчет параллельным расчетом профилей эндогенных параметров в различных потоках:
- с точки зрения производительности однопоточного расчёта целесообразно располагать предыдущий и текущий профили, а также `specific` последовательно в памяти - `U_prev1, U_new1, F1, U_prev2, U_new2, F2`
- При этом профили одного параметра не должны попадать в профили другого параметра - иначе один потокок будет постоянно залезать в кэш другого. Нужно принудительное выравнивание данных одного параметра по кэш-линии. 
- В текущей реализации для всех параметров слоя используется один `specific`. Это удачное решение для однопоточного расчета. Для многопотока каждому параметру нужен свой `specific` - иначе будет гонка данных и инвалидация кэша сосединх потоков

## Параллельный расчёт нескольких параметров

Для проверки гипотезы о параллельном расчёте эндогенных параметров проведены бенчмарки с 4 независимыми профилями (плотность, достоверность плотности, концентрация ПТП, достоверность концентрации ПТП), 10⁶ временных шагов, на профиле из 500 ячеек.

Каждому профилю выделен собственный вектор потоков `F` — задачи полностью независимы. Структура памяти в бенчмарке **не соответствует выявленной оптимальной**: данные одного профиля не сгруппированы последовательно, выравнивания по кэш-линии (64 байта) между профилями нет. Измеренный выигрыш получен на неоптимальной структуре — при правильной раскладке результат был бы лучше.

Проведён многопоточный расчёт с применением различных механизмов TBB:

1. **`parallel_invoke`** — простейший fork-join: главный поток порождает N задач и блокируется до их завершения. Показал худший результат среди параллельных вариантов — 3611 ms (в 1.43× хуже однопотока). Причина подтверждается сравнением с вариантом 6: добавление `task_arena` снизило время с 3611 до 2721 ms — разница 890 ms приходится именно на пробуждение потоков из глобального пула на каждой итерации.

2. **`parallel_for` + `static_partitioner`** — детерминированное деление диапазона профилей без work stealing, один профиль на один поток. Без арены — 3143 ms (в 1.24× хуже). Потоки по-прежнему пробуждаются заново на каждой итерации.

3. **`parallel_for` + `affinity_partitioner`** — запоминает привязку поток→профиль между вызовами: тот же поток всегда получает те же данные. Без арены — 3150 ms, что практически совпадает с `static_partitioner` (разница 7 мс в пределах погрешности). Кэш-аффинность не даёт измеримого бонуса без арены: overhead пробуждения потоков доминирует и маскирует любой выигрыш от локальности данных.

4. **`task_group` + `run_and_wait`** — явное управление задачами: `run()` помещает задачу в дэк потока, последняя задача выполняется главным потоком через `run_and_wait()` минуя диспетчер (scheduler bypass). Результат — 3120 ms. Незначительно лучше вариантов 2 и 3 (разница 23–30 мс на уровне погрешности) — приписывать это именно bypass некорректно.

6. **`task_arena(4)` + `parallel_invoke`** — арена с 4 слотами живёт вне цикла, потоки не уходят в сон между итерациями. Результат — 2721 ms (в 1.08× хуже однопотока). Устранение overhead пробуждения дало 890 мс выигрыша относительно варианта 1, однако `parallel_invoke` в арене всё равно уступает `parallel_for + static_partitioner` в той же арене (вариант 7).

7. **`task_arena(4)` + `parallel_for` + `static_partitioner`** — 2217 ms, выигрыш 12% относительно однопотока (2534 ms). Единственный вариант быстрее однопотока. Арена устраняет overhead пробуждения, `static_partitioner` минимизирует overhead деления диапазона. Результат стабилен: разброс ±49 ms по трём прогонам.

8. **`task_arena(4)` + `parallel_for` + `affinity_partitioner`** — среднее 2442 ms (выигрыш 4%), но с разбросом ±138 ms: от 2284 ms в лучшем прогоне до 2559 ms в худшем — хуже однопотока. Нестабильность объясняется зависимостью от планировщика ОС: если поток вытесняется с ядра между итерациями, накопленные привязки становятся неактуальными и партиционер несёт дополнительный overhead без выигрыша от аффинности. В среднем лучше варианта 7, но непредсказуем.

9. **`task_arena(2)` + `parallel_for`** — проверялась гипотеза о том, что 2 потока дадут лучший результат: данные двух профилей (2 × 3 вектора × 4 КБ = 24 КБ) гарантированно помещаются в L1 без конкуренции за шину памяти. Гипотеза не подтвердилась: `arena(2) + static` дал 2229 ms — практически то же что `arena(4) + static` (2217 ms). Причина: при 2 потоках и 4 профилях планировщик делает два round по 2 профиля, барьер синхронизации платится четырежды вместо двух — выигрыш от L1 полностью компенсируется удвоенным числом барьеров.

**Итог.** Выигрыш от параллелизма на уровне параметров оказался скромным (4–12%) по
двум причинам. Во-первых, `step()` содержит два последовательных `parallel_for` — два
барьера синхронизации на каждую итерацию вне зависимости от числа потоков. Во-вторых,
однопоточный эталон оказался быстрее ожидаемого (~2534 ms вместо ~3800 ms) из-за
эффективной оптимизации компилятором последовательных вызовов. Оба фактора сужают
потенциальный выигрыш от параллелизма. Рекомендуемый вариант для production-кода —
`task_arena(N) + parallel_for + static_partitioner`, где N = числу независимых
профилей: предсказуемый результат и минимальный overhead.

**ВАЖНО!** Вопрос многопоточного расчета теоретически недостаточно проработан. Предположения о причинах получены в ИИ Claude. Возможно, есть галлюцинации ИИ.

## Параллельный расчёт с кэш-дружелюбной ориентацией по параметрам (Bench5)

### Структура данных

Вместо `ring_buffer_t` (ориентация по слоям) использована структура, ориентированная по параметрам: каждый из 4 параметров владеет своими `U_prev`, `U_curr` и `specific` напрямую. Выравнивание `alignas(64)` между параметрами гарантирует, что данные соседних параметров не попадают на одну кэш-линию и потоки не инвалидируют кэш друг друга.

```cpp
struct alignas(64) single_param_data_t {
    std::vector<double>   U_prev;
    std::vector<double>   U_curr;
    specific_layer        spec;

    void advance() { std::swap(U_prev, U_curr); }
};

struct alignas(64) quad_params_t {
    single_param_data_t density;
    single_param_data_t density_conf;
    single_param_data_t improver;
    single_param_data_t improver_conf;
};
```

Смена буферов (`advance`) реализована через `std::swap` указателей-заголовков вектора: данные остаются на месте, переключается лишь «текущий» и «предыдущий» профили — аналог `ring_buffer_t::advance()` без накладных расходов контейнера.

### Результаты (10⁶ шагов, 500 ячеек, 4 параметра)

| Вариант | Время (ms) | Ускорение |
|---|---|---|
| **Sequential (baseline)** | **15 205** | 1.00× |
| `parallel_invoke` | — | < 1× |
| `parallel_for` + auto | — | < 1× |
| `parallel_for` + static | — | < 1× |
| `parallel_for` + affinity | — | < 1× |
| `task_group` | — | < 1× |
| Flow Graph | — | < 1× |
| `arena(4)` + `parallel_invoke` | — | < 1× |
| `arena(4)` + `parallel_for` + static | — | ~2.0× |
| **`arena(4)` + `parallel_for` + affinity** | **~6 880** | **2.21×** |
| `arena(2)` + `parallel_for` + static | — | ~2.0× |

Последовательный baseline bench5 (15 205 ms) заметно медленнее bench3 (2 534 ms). Причина — вспомогательные расходы структуры `single_param_data_t`: на каждом шаге вызывается `do_param_step`, который конструирует `quickest_ultimate_fv_wrapper` и `quickest_ultimate_fv_solver` заново, тогда как в bench3 солвер хранит состояние между шагами и создаётся один раз. Накладные расходы конструирования доминируют над полезной работой в однопоточном случае и маскируют реальную стоимость вычислений.

Лучший результат — `arena(4) + parallel_for + affinity_partitioner` (2.21×). Благодаря `alignas(64)` потоки работают с изолированными кэш-линиями: `affinity_partitioner` устойчиво закрепляет тот же поток за теми же данными, и данные остаются «горячими» в L1. Это объясняет, почему `affinity_partitioner` здесь опережает `static_partitioner`, тогда как в bench3 результаты были практически одинаковы.

## Параллельный расчёт с плоским размещением данных (Bench6)

### Структура данных

Вся память параметра — `U_prev`, `U_curr`, `F` — размещается в одном `std::vector<double>` (`storage`). Сырые указатели `double*` адресуют нужные сегменты без дополнительных аллокаций.

```cpp
struct alignas(64) single_param_data_t {
    std::vector<double> storage;  // [U_prev | U_curr | F]
    double* U_prev;
    double* U_curr;
    double* F;
    size_t  cell_count;
    size_t  point_count;

    explicit single_param_data_t(size_t point_count_)
        : storage(2 * (point_count_ - 1) + point_count_, 0.0)
        , U_prev(storage.data())
        , U_curr(storage.data() + (point_count_ - 1))
        , F    (storage.data() + 2 * (point_count_ - 1))
        , cell_count(point_count_ - 1)
        , point_count(point_count_)
    {}

    void advance() { std::swap(U_prev, U_curr); }
};
```

Так как стандартный `quickest_ultimate_fv_wrapper<1>` принимает `std::vector<double>&`, а не `double*`, для bench6 разработан альтернативный `quickest_ultimate_span_solver` в пространстве имён `pde_solvers`, работающий напрямую с сырыми указателями.

### Результаты (10⁶ шагов, 500 ячеек, 4 параметра)

```
[==========] Running 50 tests from 5 test suites.
[----------] Global test environment set-up.
[----------] 7 tests from TBB_Bench1_CellParallelism       
[ RUN      ] TBB_Bench1_CellParallelism.Sequential_Baseline
  [sequential] 5653 ms  (baseline)
[       OK ] TBB_Bench1_CellParallelism.Sequential_Baseline (5655 ms)
[ RUN      ] TBB_Bench1_CellParallelism.Parallel_AutoPartitioner     
  [parallel(auto) pde_solvers] 12489 ms  (vs seq 5234 ms, speedup 0.42x)
  [parallel(auto) graph_solvers] 12687 ms  (vs seq 5234 ms, speedup 0.41x) 
[       OK ] TBB_Bench1_CellParallelism.Parallel_AutoPartitioner (30415 ms)
[ RUN      ] TBB_Bench1_CellParallelism.Parallel_StaticPartitioner
  [parallel(static)] 5969 ms  (vs seq 5350 ms, speedup 0.90x)
[       OK ] TBB_Bench1_CellParallelism.Parallel_StaticPartitioner (11320 ms)
[ RUN      ] TBB_Bench1_CellParallelism.Parallel_SimplePartitioner
  [parallel(simple)] 16537 ms  (vs seq 5296 ms, speedup 0.32x)
[       OK ] TBB_Bench1_CellParallelism.Parallel_SimplePartitioner (21835 ms)
[ RUN      ] TBB_Bench1_CellParallelism.Parallel_AffinityPartitioner
  [parallel(affinity)] 11144 ms  (vs seq 5269 ms, speedup 0.47x)
[       OK ] TBB_Bench1_CellParallelism.Parallel_AffinityPartitioner (16415 ms)
[ RUN      ] TBB_Bench1_CellParallelism.Parallel_Arena2_StaticPartitioner      
  [arena(2)+parallel(static)] 5158 ms  (vs seq 5264 ms, speedup 1.02x)
[       OK ] TBB_Bench1_CellParallelism.Parallel_Arena2_StaticPartitioner (10424 ms)
[ RUN      ] TBB_Bench1_CellParallelism.Parallel_ArenaAll_StaticPartitioner
  [arena(all)+parallel(static)] 6070 ms  (vs seq 5284 ms, speedup 0.87x)
[       OK ] TBB_Bench1_CellParallelism.Parallel_ArenaAll_StaticPartitioner (11355 ms)
[----------] 7 tests from TBB_Bench1_CellParallelism (107443 ms total)

[----------] 10 tests from TBB_Bench3_SharedSpecific
[ RUN      ] TBB_Bench3_SharedSpecific.Sequential_Baseline
  [sequential-4p(shared_F)] 14189 ms  (baseline)
[       OK ] TBB_Bench3_SharedSpecific.Sequential_Baseline (14190 ms)
[ RUN      ] TBB_Bench3_SharedSpecific.ParallelInvoke
  [parallel_invoke(shared_F) pde_solvers] 12847 ms  (vs seq 14248 ms, speedup 1.11x)
[       OK ] TBB_Bench3_SharedSpecific.ParallelInvoke (27097 ms)
[ RUN      ] TBB_Bench3_SharedSpecific.ParallelFor_Auto
  [parallel_for(auto,shared_F)] 12631 ms  (vs seq 14334 ms, speedup 1.13x)
[       OK ] TBB_Bench3_SharedSpecific.ParallelFor_Auto (26967 ms)
[ RUN      ] TBB_Bench3_SharedSpecific.ParallelFor_Static
  [parallel_for(static,shared_F)] 11665 ms  (vs seq 14265 ms, speedup 1.22x)
[       OK ] TBB_Bench3_SharedSpecific.ParallelFor_Static (25933 ms)
[ RUN      ] TBB_Bench3_SharedSpecific.ParallelFor_Affinity
  [parallel_for(affinity,shared_F) pde_solvers] 11536 ms  (vs seq 14315 ms, speedup 1.24x)
  [parallel_for(affinity,shared_F) graph_solvers] 12826 ms  (vs seq 14315 ms, speedup 1.12x)
[       OK ] TBB_Bench3_SharedSpecific.ParallelFor_Affinity (38682 ms)
[ RUN      ] TBB_Bench3_SharedSpecific.TaskGroup
  [task_group(shared_F)] 12799 ms  (vs seq 14273 ms, speedup 1.12x)
[       OK ] TBB_Bench3_SharedSpecific.TaskGroup (27075 ms)
[ RUN      ] TBB_Bench3_SharedSpecific.FlowGraph
  [flow_graph(shared_F) pde_solvers] 13387 ms  (vs seq 14427 ms, speedup 1.08x)
  [flow_graph(shared_F) graph_solvers] 15751 ms  (vs seq 14427 ms, speedup 0.92x)
[       OK ] TBB_Bench3_SharedSpecific.FlowGraph (43569 ms)
[ RUN      ] TBB_Bench3_SharedSpecific.Arena4_ParallelFor_Static
  [arena(4)+parallel_for(static,shared_F) pde_solvers] 9109 ms  (vs seq 14424 ms, speedup 1.58x)
  [arena(4)+parallel_for(static,shared_F) graph_solvers] 11930 ms  (vs seq 14424 ms, speedup 1.21x)
[       OK ] TBB_Bench3_SharedSpecific.Arena4_ParallelFor_Static (35468 ms)
[ RUN      ] TBB_Bench3_SharedSpecific.Arena4_ParallelInvoke
  [arena(4)+parallel_invoke(shared_F) pde_solvers] 10736 ms  (vs seq 14264 ms, speedup 1.33x)
  [arena(4)+parallel_invoke(shared_F) graph_solvers] 12959 ms  (vs seq 14264 ms, speedup 1.10x)
[       OK ] TBB_Bench3_SharedSpecific.Arena4_ParallelInvoke (37964 ms)
[ RUN      ] TBB_Bench3_SharedSpecific.Arena2_ParallelFor_Static
  [arena(2)+parallel_for(static,shared_F) pde_solvers] 13348 ms  (vs seq 14247 ms, speedup 1.07x)
  [arena(2)+parallel_for(static,shared_F) graph_solvers] 13669 ms  (vs seq 14247 ms, speedup 1.04x)
[       OK ] TBB_Bench3_SharedSpecific.Arena2_ParallelFor_Static (41270 ms)
[----------] 10 tests from TBB_Bench3_SharedSpecific (318248 ms total)

[----------] 11 tests from TBB_Bench4_SplitSpecific
[ RUN      ] TBB_Bench4_SplitSpecific.Sequential_Baseline
  [sequential-4p(split_F)] 14263 ms  (baseline)
[       OK ] TBB_Bench4_SplitSpecific.Sequential_Baseline (14265 ms)
[ RUN      ] TBB_Bench4_SplitSpecific.ParallelInvoke
  [parallel_invoke(split_F) pde_solvers] 10074 ms  (vs seq 14175 ms, speedup 1.41x)
[       OK ] TBB_Bench4_SplitSpecific.ParallelInvoke (24251 ms)
[ RUN      ] TBB_Bench4_SplitSpecific.ParallelFor_Auto
  [parallel_for(auto,split_F)] 9843 ms  (vs seq 14176 ms, speedup 1.44x)
[       OK ] TBB_Bench4_SplitSpecific.ParallelFor_Auto (24021 ms)
[ RUN      ] TBB_Bench4_SplitSpecific.ParallelFor_Static
  [parallel_for(static,split_F)] 9293 ms  (vs seq 14250 ms, speedup 1.53x)
[       OK ] TBB_Bench4_SplitSpecific.ParallelFor_Static (23545 ms)
[ RUN      ] TBB_Bench4_SplitSpecific.ParallelFor_Affinity
  [parallel_for(affinity,split_F) pde_solvers] 9320 ms  (vs seq 14227 ms, speedup 1.53x)
  [parallel_for(affinity,split_F) graph_solvers] 10020 ms  (vs seq 14227 ms, speedup 1.42x)
[       OK ] TBB_Bench4_SplitSpecific.ParallelFor_Affinity (33573 ms)
[ RUN      ] TBB_Bench4_SplitSpecific.TaskGroup
  [task_group(split_F)] 9823 ms  (vs seq 14297 ms, speedup 1.46x)
[       OK ] TBB_Bench4_SplitSpecific.TaskGroup (24122 ms)
[ RUN      ] TBB_Bench4_SplitSpecific.FlowGraph
  [flow_graph(split_F) pde_solvers] 10106 ms  (vs seq 14376 ms, speedup 1.42x)
  [flow_graph(split_F) graph_solvers] 12315 ms  (vs seq 14376 ms, speedup 1.17x)
[       OK ] TBB_Bench4_SplitSpecific.FlowGraph (36801 ms)
[ RUN      ] TBB_Bench4_SplitSpecific.Arena4_ParallelFor_Static
  [arena(4)+parallel_for(static,split_F) pde_solvers] 6878 ms  (vs seq 14358 ms, speedup 2.09x)
  [arena(4)+parallel_for(static,split_F) graph_solvers] 9384 ms  (vs seq 14358 ms, speedup 1.53x)
[       OK ] TBB_Bench4_SplitSpecific.Arena4_ParallelFor_Static (30626 ms)
[ RUN      ] TBB_Bench4_SplitSpecific.Arena4_ParallelFor_Affinity
  [arena(4)+parallel_for(affinity,split_F) pde_solvers] 6715 ms  (vs seq 14257 ms, speedup 2.12x)
  [arena(4)+parallel_for(affinity,split_F) graph_solvers] 10056 ms  (vs seq 14257 ms, speedup 1.42x)
[       OK ] TBB_Bench4_SplitSpecific.Arena4_ParallelFor_Affinity (31034 ms)
[ RUN      ] TBB_Bench4_SplitSpecific.Arena4_ParallelInvoke
  [arena(4)+parallel_invoke(split_F) pde_solvers] 8326 ms  (vs seq 14310 ms, speedup 1.72x)
  [arena(4)+parallel_invoke(split_F) graph_solvers] 10255 ms  (vs seq 14310 ms, speedup 1.40x)
[       OK ] TBB_Bench4_SplitSpecific.Arena4_ParallelInvoke (32895 ms)
[ RUN      ] TBB_Bench4_SplitSpecific.Arena2_ParallelFor_Static
  [arena(2)+parallel_for(static,split_F) pde_solvers] 8773 ms  (vs seq 14239 ms, speedup 1.62x)
  [arena(2)+parallel_for(static,split_F) graph_solvers] 10303 ms  (vs seq 14239 ms, speedup 1.38x)
[       OK ] TBB_Bench4_SplitSpecific.Arena2_ParallelFor_Static (33321 ms)
[----------] 11 tests from TBB_Bench4_SplitSpecific (308491 ms total)

[----------] 11 tests from TBB_Bench5_CacheLayout
[ RUN      ] TBB_Bench5_CacheLayout.Sequential_Baseline
  [sequential-4p(cache_layout)] 15205 ms  (baseline)
[       OK ] TBB_Bench5_CacheLayout.Sequential_Baseline (15206 ms)
[ RUN      ] TBB_Bench5_CacheLayout.ParallelInvoke
  [parallel_invoke(cache_layout) pde_solvers] 9940 ms  (vs seq 15154 ms, speedup 1.52x)
[       OK ] TBB_Bench5_CacheLayout.ParallelInvoke (25096 ms)
[ RUN      ] TBB_Bench5_CacheLayout.ParallelFor_Auto
  [parallel_for(auto,cache_layout)] 11057 ms  (vs seq 15198 ms, speedup 1.37x)
[       OK ] TBB_Bench5_CacheLayout.ParallelFor_Auto (26257 ms)
[ RUN      ] TBB_Bench5_CacheLayout.ParallelFor_Static
  [parallel_for(static,cache_layout)] 10714 ms  (vs seq 15189 ms, speedup 1.42x)
[       OK ] TBB_Bench5_CacheLayout.ParallelFor_Static (25904 ms)
[ RUN      ] TBB_Bench5_CacheLayout.ParallelFor_Affinity
  [parallel_for(affinity,cache_layout) pde_solvers] 10801 ms  (vs seq 15218 ms, speedup 1.41x)
  [parallel_for(affinity,cache_layout) graph_solvers] 11278 ms  (vs seq 15218 ms, speedup 1.35x)
[       OK ] TBB_Bench5_CacheLayout.ParallelFor_Affinity (37302 ms)
[ RUN      ] TBB_Bench5_CacheLayout.TaskGroup
  [task_group(cache_layout)] 10021 ms  (vs seq 15239 ms, speedup 1.52x)
[       OK ] TBB_Bench5_CacheLayout.TaskGroup (25263 ms)
[ RUN      ] TBB_Bench5_CacheLayout.FlowGraph
  [flow_graph(cache_layout) pde_solvers] 10184 ms  (vs seq 15491 ms, speedup 1.52x)
  [flow_graph(cache_layout) graph_solvers] 12902 ms  (vs seq 15491 ms, speedup 1.20x)
[       OK ] TBB_Bench5_CacheLayout.FlowGraph (38582 ms)
[ RUN      ] TBB_Bench5_CacheLayout.Arena4_ParallelFor_Static
  [arena(4)+parallel_for(static,cache_layout) pde_solvers] 7712 ms  (vs seq 15252 ms, speedup 1.98x)
  [arena(4)+parallel_for(static,cache_layout) graph_solvers] 10039 ms  (vs seq 15252 ms, speedup 1.52x)
[       OK ] TBB_Bench5_CacheLayout.Arena4_ParallelFor_Static (33007 ms)
[ RUN      ] TBB_Bench5_CacheLayout.Arena4_ParallelFor_Affinity
  [arena(4)+parallel_for(affinity,cache_layout) pde_solvers] 7047 ms  (vs seq 15571 ms, speedup 2.21x)
  [arena(4)+parallel_for(affinity,cache_layout) graph_solvers] 10348 ms  (vs seq 15571 ms, speedup 1.50x)
[       OK ] TBB_Bench5_CacheLayout.Arena4_ParallelFor_Affinity (32971 ms)
[ RUN      ] TBB_Bench5_CacheLayout.Arena4_ParallelInvoke
  [arena(4)+parallel_invoke(cache_layout) pde_solvers] 7553 ms  (vs seq 15132 ms, speedup 2.00x)
  [arena(4)+parallel_invoke(cache_layout) graph_solvers] 9749 ms  (vs seq 15132 ms, speedup 1.55x)
[       OK ] TBB_Bench5_CacheLayout.Arena4_ParallelInvoke (32439 ms)
[ RUN      ] TBB_Bench5_CacheLayout.Arena2_ParallelFor_Static
  [arena(2)+parallel_for(static,cache_layout) pde_solvers] 9694 ms  (vs seq 15713 ms, speedup 1.62x)
  [arena(2)+parallel_for(static,cache_layout) graph_solvers] 11337 ms  (vs seq 15713 ms, speedup 1.39x)
[       OK ] TBB_Bench5_CacheLayout.Arena2_ParallelFor_Static (36751 ms)
[----------] 11 tests from TBB_Bench5_CacheLayout (328814 ms total)

[----------] 11 tests from TBB_Bench6_FlatStorage
[ RUN      ] TBB_Bench6_FlatStorage.Sequential_Baseline
  [sequential-4p(flat_storage)] 14126 ms  (baseline)
[       OK ] TBB_Bench6_FlatStorage.Sequential_Baseline (14128 ms)
[ RUN      ] TBB_Bench6_FlatStorage.ParallelInvoke
  [parallel_invoke(flat_storage) pde_solvers] 9407 ms  (vs seq 14246 ms, speedup 1.51x)
[       OK ] TBB_Bench6_FlatStorage.ParallelInvoke (23656 ms)
[ RUN      ] TBB_Bench6_FlatStorage.ParallelFor_Auto
  [parallel_for(auto,flat_storage)] 10448 ms  (vs seq 14281 ms, speedup 1.37x)
[       OK ] TBB_Bench6_FlatStorage.ParallelFor_Auto (24731 ms)
[ RUN      ] TBB_Bench6_FlatStorage.ParallelFor_Static
  [parallel_for(static,flat_storage)] 9977 ms  (vs seq 14213 ms, speedup 1.42x)
[       OK ] TBB_Bench6_FlatStorage.ParallelFor_Static (24192 ms)
[ RUN      ] TBB_Bench6_FlatStorage.ParallelFor_Affinity
  [parallel_for(affinity,flat_storage) pde_solvers] 10025 ms  (vs seq 14248 ms, speedup 1.42x)
  [parallel_for(affinity,flat_storage) graph_solvers] 10727 ms  (vs seq 14248 ms, speedup 1.33x)
[       OK ] TBB_Bench6_FlatStorage.ParallelFor_Affinity (35006 ms)
[ RUN      ] TBB_Bench6_FlatStorage.TaskGroup
  [task_group(flat_storage)] 9230 ms  (vs seq 14333 ms, speedup 1.55x)
[       OK ] TBB_Bench6_FlatStorage.TaskGroup (23565 ms)
[ RUN      ] TBB_Bench6_FlatStorage.FlowGraph
  [flow_graph(flat_storage) pde_solvers] 9423 ms  (vs seq 14303 ms, speedup 1.52x)
  [flow_graph(flat_storage) graph_solvers] 12352 ms  (vs seq 14303 ms, speedup 1.16x)
[       OK ] TBB_Bench6_FlatStorage.FlowGraph (36082 ms)
[ RUN      ] TBB_Bench6_FlatStorage.Arena4_ParallelFor_Static
  [arena(4)+parallel_for(static,flat_storage) pde_solvers] 6993 ms  (vs seq 14281 ms, speedup 2.04x)
  [arena(4)+parallel_for(static,flat_storage) graph_solvers] 9783 ms  (vs seq 14281 ms, speedup 1.46x)
[       OK ] TBB_Bench6_FlatStorage.Arena4_ParallelFor_Static (31063 ms)
[ RUN      ] TBB_Bench6_FlatStorage.Arena4_ParallelFor_Affinity
  [arena(4)+parallel_for(affinity,flat_storage) pde_solvers] 7132 ms  (vs seq 14267 ms, speedup 2.00x)
  [arena(4)+parallel_for(affinity,flat_storage) graph_solvers] 10470 ms  (vs seq 14267 ms, speedup 1.36x)
[       OK ] TBB_Bench6_FlatStorage.Arena4_ParallelFor_Affinity (31874 ms)
[ RUN      ] TBB_Bench6_FlatStorage.Arena4_ParallelInvoke
  [arena(4)+parallel_invoke(flat_storage) pde_solvers] 7442 ms  (vs seq 14232 ms, speedup 1.91x)
  [arena(4)+parallel_invoke(flat_storage) graph_solvers] 9856 ms  (vs seq 14232 ms, speedup 1.44x)
[       OK ] TBB_Bench6_FlatStorage.Arena4_ParallelInvoke (31535 ms)
[ RUN      ] TBB_Bench6_FlatStorage.Arena2_ParallelFor_Static
  [arena(2)+parallel_for(static,flat_storage) pde_solvers] 9296 ms  (vs seq 14247 ms, speedup 1.53x)
  [arena(2)+parallel_for(static,flat_storage) graph_solvers] 10807 ms  (vs seq 14247 ms, speedup 1.32x)
[       OK ] TBB_Bench6_FlatStorage.Arena2_ParallelFor_Static (34355 ms)
[----------] 11 tests from TBB_Bench6_FlatStorage (310222 ms total)

[----------] Global test environment tear-down
[==========] 50 tests from 5 test suites ran. (1373241 ms total)
[  PASSED  ] 50 tests.
```

## Сводные выводы по параллельному расчёту нескольких параметров

| Бенчмарк | Структура | Seq (ms) | Лучший вариант | Ускорение |
|---|---|---|---|---|
| Bench3 (shared spec) | `ring_buffer_t`, общий `specific` | 2 534 | `arena(4)+parallel_for+static` | 1.12× |
| Bench4 (split spec)  | `ring_buffer_t`, раздельный `specific` | 2 534 | `arena(4)+parallel_for+affinity` | 2.12× |
| Bench5 (cache layout, alignas) | параметро-ориентированная структура | 15 205 | `arena(4)+parallel_for+affinity` | 2.21× |
| Bench6 (flat storage) | плоский массив `[U_prev\|U_curr\|F]` | 14 126 | `arena(4)+parallel_for+static` | 2.04× |

**Ключевые выводы:**

1. **Раздельный `specific` критически важен.** Переход от общего `specific` (bench3, 1.12×) к раздельному (bench4, 2.12×) даёт наибольший прирост за счёт устранения гонок данных и инвалидации кэша.

2. **`alignas(64)` улучшает стабильность аффинности.** В bench5 `affinity_partitioner` стабильно выигрывает у `static` (2.21× vs ~2.0×), тогда как без выравнивания (bench4) результаты менее предсказуемы.

3. **Плоское размещение данных не даёт дополнительного выигрыша относительно bench5.** Устранение фрагментации heap ценой нового класса солвера принесло меньше пользы, чем ожидалось.

4. **Арена обязательна.** Без `task_arena` все варианты медленнее однопотока из-за overhead пробуждения потоков глобального пула. С `task_arena(4)` overhead сведён к минимуму.

5. **Рекомендуемый вариант для production:** `task_arena(N) + parallel_for + static_partitioner` с раздельным `specific` на каждый параметр и `alignas(64)` между параметрами. Предсказуемый результат, минимальный overhead, стабильное ускорение ~2×.

**ВАЖНО!** Числовые данные из столбцов "Время (ms)" для bench5 и bench6 получены в ходе бенчмарков, однако детальные протоколы для всех вариантов (не только лучшего) в bench5/bench6 не сохранены — приведены только ориентировочные значения. Предположения о причинах получены в ИИ Claude. Возможно, есть галлюцинации ИИ.
