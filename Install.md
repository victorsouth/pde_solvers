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
cmake --build --preset windows-msvc-release-install
```

## Debug

- VS Code: use launch configurations from `.vscode/launch.json`.
- Visual Studio 2022: open folder as CMake project and pick `windows-msvc-debug`.
