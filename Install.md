# Install

## One-time setup

### Windows x64
- Install Visual Studio 17 2022 (Desktop development with C++).
- Install MSYS2 to `C:/msys64` and install UCRT64 components:

```powershell
C:\msys64\usr\bin\bash -lc "pacman -Syu --noconfirm"
C:\msys64\usr\bin\bash -lc "pacman -S --needed --noconfirm mingw-w64-ucrt-x86_64-toolchain mingw-w64-ucrt-x86_64-clang mingw-w64-ucrt-x86_64-gdb mingw-w64-ucrt-x86_64-lldb mingw-w64-ucrt-x86_64-ninja"
C:\msys64\usr\bin\bash -lc "pacman -S --needed --noconfirm mingw-w64-ucrt-x86_64-gtest mingw-w64-ucrt-x86_64-eigen3"
```

- Install vcpkg to `C:/vcpkg`.

### Linux
- Install `gcc`, `clang`, `gdb`, `ninja`.
- Install LLDB components for Clang debugging and Eigen printers: `lldb`, Python bindings for LLDB (for example `python3-module-lldb<version>`), and LLDB DAP adapter (`lldb-dap` or `lldb-vscode` from LLVM tools package).
- Install vcpkg to `${HOME}/vcpkg`.

### macOS
- Install Xcode Command Line Tools (`clang`, `lldb`).
- Install `ninja` (for example via Homebrew).
- Install vcpkg to `${HOME}/vcpkg`.

## Supported compiler matrix

- Windows: MSVC, GCC, Clang
- Linux: GCC, Clang
- macOS: Clang

## vcpkg dependencies

`pde_solvers` uses:
- GTest
- Eigen3
- Boost is not used in this project.

Install for Windows (x64, MSVC only):

```powershell
C:\vcpkg\vcpkg install gtest:x64-windows eigen3:x64-windows
```

For Windows GCC/Clang, dependencies are taken from MSYS2 UCRT64 (`C:/msys64/ucrt64`) via CMake presets.

Install for Linux:

```bash
$HOME/vcpkg/vcpkg install gtest:x64-linux eigen3:x64-linux
```

Install for macOS:

```bash
$HOME/vcpkg/vcpkg install gtest:x64-osx eigen3:x64-osx
```

## Install layout

Each preset installs to a compiler/config specific directory:
- Windows: C:/install/<project>/<compiler>/<Debug|Release>
- Linux/macOS: $HOME/install/<project>/<compiler>/<Debug|Release>

This keeps Debug and Release artifacts side-by-side without renaming binaries.
## Configure/build/test

Use one of presets:
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

Example:

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

Install uses user-scoped prefixes from presets, so admin privileges are typically not required.

## Debug

- VS Code: use launch configurations from `.vscode/launch.json`.
- Visual Studio 2022: open folder as CMake project and pick `windows-msvc-debug`.

## Eigen pretty-printers (shared install root)

All debuggers use one shared location from `EIGEN_PRINTERS_ROOT`.
Default location:
- Windows: `C:/install/eigen-pretty-printers`
- Linux/macOS: `$HOME/install/eigen-pretty-printers`

Required layout:
- `${EIGEN_PRINTERS_ROOT}/gdb/eigen.gdbinit`
- `${EIGEN_PRINTERS_ROOT}/gdb/printers.py`
- `${EIGEN_PRINTERS_ROOT}/lldb/eigen_printers.py`
- `${EIGEN_PRINTERS_ROOT}/natvis/Eigen.natvis`

Set `EIGEN_PRINTERS_ROOT`:

```powershell
# Windows (PowerShell, current user)
[Environment]::SetEnvironmentVariable("EIGEN_PRINTERS_ROOT", "C:/install/eigen-pretty-printers", "User")
```

```bash
# Linux/macOS
export EIGEN_PRINTERS_ROOT="$HOME/install/eigen-pretty-printers"
```

Validation checklist:
- VS Code + `windows-msvc-debug`: `Eigen::Matrix` displayed via natvis.
- VS Code + `windows-gcc-debug`/`linux-gcc-debug`: `Eigen::Matrix` expanded with GDB printer.
- VS Code + `windows-clang-debug`/`linux-clang-debug`/`macos-clang-debug`: LLDB summaries/synthetic children are available.
- Visual Studio 2022 (`windows-msvc-debug`): natvis is loaded automatically from CMake (`CMAKE_VS_NATVIS_PATH`).