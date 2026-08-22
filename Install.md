# Установка

## Первичная настройка

### Windows x64
- Установите Visual Studio 17 2022 (Desktop development with C++).
- Установите MSYS2 в `C:/msys64` (скачивание: https://www.msys2.org/) и добавьте компоненты UCRT64:

```powershell
C:\msys64\usr\bin\bash -lc "pacman -Syu --noconfirm"
C:\msys64\usr\bin\bash -lc "pacman -S --needed --noconfirm mingw-w64-ucrt-x86_64-toolchain mingw-w64-ucrt-x86_64-clang mingw-w64-ucrt-x86_64-gdb mingw-w64-ucrt-x86_64-lldb mingw-w64-ucrt-x86_64-lldb-dap mingw-w64-ucrt-x86_64-ninja"
C:\msys64\usr\bin\bash -lc "pacman -S --needed --noconfirm mingw-w64-ucrt-x86_64-gtest mingw-w64-ucrt-x86_64-eigen3"
```

- Установите vcpkg в `C:/vcpkg`:

```powershell
git clone https://github.com/microsoft/vcpkg C:\vcpkg
C:\vcpkg\bootstrap-vcpkg.bat
C:\vcpkg\vcpkg integrate install
```

Команда `vcpkg integrate install` опциональна и полезна для локальной работы в Visual Studio/MSBuild.

### Linux
- Установите инструменты: `git`, `cmake`, `ninja`, `gcc`, `clang`, `gdb`, `lldb`.
- Для отладки Clang и принтеров Eigen установите компоненты LLDB: `lldb`, Python-привязки для LLDB (например, `python3-module-lldb<version>`) и LLDB DAP-адаптер (`lldb-dap` или `lldb-vscode` из пакета LLVM tools).
- Установите vcpkg в `${HOME}/vcpkg`:

```bash
git clone https://github.com/microsoft/vcpkg "$HOME/vcpkg"
"$HOME/vcpkg/bootstrap-vcpkg.sh"
```

### macOS
- Установите Xcode Command Line Tools (`clang`, `lldb`, `git`).
- Установите инструменты: `cmake`, `ninja` (например, через Homebrew).
- Установите vcpkg в `${HOME}/vcpkg`:

```bash
git clone https://github.com/microsoft/vcpkg "$HOME/vcpkg"
"$HOME/vcpkg/bootstrap-vcpkg.sh"
```

## Поддерживаемые компиляторы

- Windows: MSVC, GCC, Clang
- Linux: GCC, Clang
- macOS: Clang

## Зависимости vcpkg

`pde_solvers` использует:
- GTest
- Eigen3
- Boost Headers

Установка для Windows (x64, только MSVC):

```powershell
C:\vcpkg\vcpkg install gtest:x64-windows eigen3:x64-windows boost:x64-windows
```

Для Windows GCC/Clang зависимости берутся из MSYS2 UCRT64 (`C:/msys64/ucrt64`) через пресеты CMake.

Установка для Linux:

```bash
$HOME/vcpkg/vcpkg install gtest:x64-linux eigen3:x64-linux boost:x64-linux
```

Установка для macOS:

```bash
$HOME/vcpkg/vcpkg install gtest:x64-osx eigen3:x64-osx boost:x64-osx
```

## Структура установки

Каждый пресет устанавливает артефакты в отдельный путь по компилятору и конфигурации:
- Windows: C:/install/<project>/<compiler>/<Debug|Release>
- Linux/macOS: $HOME/install/<project>/<compiler>/<Debug|Release>

Это позволяет держать Debug и Release рядом без переименования бинарников.
## Конфигурация / сборка / тесты

Используйте один из пресетов:
- `windows-msvc-debug`
- `windows-gcc-debug`
- `windows-clang-debug`
- `windows-msvc-release`
- `windows-gcc-release`
- `windows-clang-release`
- `linux-gcc-debug`
- `linux-clang-debug`
- `linux-gcc-release`
- `linux-clang-release`
- `macos-clang-debug`
- `macos-clang-release`

Примеры:

```powershell
cmake --preset windows-gcc-debug
cmake --build --preset windows-gcc-debug
ctest --preset windows-gcc-debug --output-on-failure
```

```bash
cmake --preset linux-clang-debug
cmake --build --preset linux-clang-debug
ctest --preset linux-clang-debug --output-on-failure
```

```powershell
cmake --preset windows-msvc-release
cmake --build --preset windows-msvc-release
cmake --install --preset windows-msvc-release
```

В установке используются пользовательские префиксы из пресетов, поэтому права администратора обычно не требуются.

## Отладка

- VS Code: используйте конфигурации запуска из `.vscode/launch.json`.
- Visual Studio 2022: откройте папку как CMake-проект и выберите `windows-msvc-debug`.

## Pretty-printers Eigen (общий корень установки)

Все отладчики используют единое расположение из `EIGEN_PRINTERS_ROOT`.
Путь по умолчанию:
- Windows: `C:/install/eigen-pretty-printers`
- Linux/macOS: `$HOME/install/eigen-pretty-printers`

Требуемая структура:
- `${EIGEN_PRINTERS_ROOT}/gdb/eigen.gdbinit`
- `${EIGEN_PRINTERS_ROOT}/gdb/printers.py`
- `${EIGEN_PRINTERS_ROOT}/lldb/eigen_printers.py`
- `${EIGEN_PRINTERS_ROOT}/natvis/Eigen.natvis`

Проверочный список:
- VS Code + `windows-msvc-debug`: `Eigen::Matrix` отображается через natvis.
- VS Code + `windows-gcc-debug`/`linux-gcc-debug`: `Eigen::Matrix` разворачивается через принтер GDB.
- VS Code + `windows-clang-debug`/`linux-clang-debug`/`macos-clang-debug`: доступны summaries/synthetic children в LLDB.