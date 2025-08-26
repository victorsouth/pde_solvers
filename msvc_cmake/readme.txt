Подготовка Debug

Должны быть предсобраны и заинстоллены с помощью CMAKE:
- Eigen Debug
- Gtest Debug
- fixed_solvers Debug

Подробности см. fixed_solvers/msvc_cmake/readme.txt
(ЕСЛИ В БИБЛИОТЕКАХ ПОЯВИЛИСЬ ИЗМЕНЕНИЯ, КОТОРЫЕ НЕОБХОДИМЫ ДЛЯ РАБОТЫ ТЕКУЩЕЙ БИБЛИОТЕКИ, ТО ИХ НАДО ПЕРЕСОБИРАТЬ И ПЕРЕУСТАНАВЛИВАТЬ)

При всех собранных зависимостях используем типовую команду (запускать из pde_solvers/msvc_cmake)

cmake .. -G "Visual Studio 17 2022" -DCMAKE_PREFIX_PATH="C:/Program Files (x86)/googletest-distribution"  -DCMAKE_MSVC_RUNTIME_LIBRARY=MultiThreadedDebugDLL -DCMAKE_BUILD_TYPE=Debug -DPDE_SOLVERS_BUILD_TESTS=ON

В папке pde_solvers/msvc_cmake будет создано решение MSVC (.sln). Его открываем, работаем с ним как обычно.

Для установки с целью использования в другом зависимом проекте (graph_solvers и др.):

cmake --build . --config Debug -j 8
cmake --install . --config Debug



Подготовка Release

Должны быть предсобраны и заинстоллены с помощью CMAKE:
- Eigen Release
- Gtest Release
- fixed_solvers Release

Подробности см. fixed_solvers/msvc_cmake/readme.txt
(ЕСЛИ В БИБЛИОТЕКАХ ПОЯВИЛИСЬ ИЗМЕНЕНИЯ, КОТОРЫЕ НЕОБХОДИМЫ ДЛЯ РАБОТЫ ТЕКУЩЕЙ БИБЛИОТЕКИ, ТО ИХ НАДО ПЕРЕСОБИРАТЬ И ПЕРЕУСТАНАВЛИВАТЬ)

При всех собранных зависимостях используем типовую команду (запускать из pde_solvers/msvc_cmake)

cmake .. -G "Visual Studio 17 2022" -DCMAKE_PREFIX_PATH="C:/Program Files (x86)/googletest-distribution"  -DCMAKE_MSVC_RUNTIME_LIBRARY=MultiThreadedDLL -DCMAKE_BUILD_TYPE=Release -DPDE_SOLVERS_BUILD_TESTS=ON

В папке pde_solvers/msvc_cmake будет создано решение MSVC (.sln). Его открываем, работаем с ним как обычно.

Для установки с целью использования в другом зависимом проекте (graph_solvers и др.):

cmake --build . --config Release -j 8
cmake --install . --config Release

