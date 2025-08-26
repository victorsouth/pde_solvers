���������� Debug

������ ���� ����������� � ������������ � ������� CMAKE:
- Eigen Debug
- Gtest Debug
- fixed_solvers Debug

����������� ��. fixed_solvers/msvc_cmake/readme.txt
(���� � ����������� ��������� ���������, ������� ���������� ��� ������ ������� ����������, �� �� ���� ������������ � �����������������)

��� ���� ��������� ������������ ���������� ������� ������� (��������� �� pde_solvers/msvc_cmake)

cmake .. -G "Visual Studio 17 2022" -DCMAKE_PREFIX_PATH="C:/Program Files (x86)/googletest-distribution"  -DCMAKE_MSVC_RUNTIME_LIBRARY=MultiThreadedDebugDLL -DCMAKE_BUILD_TYPE=Debug -DPDE_SOLVERS_BUILD_TESTS=ON

� ����� pde_solvers/msvc_cmake ����� ������� ������� MSVC (.sln). ��� ���������, �������� � ��� ��� ������.

��� ��������� � ����� ������������� � ������ ��������� ������� (graph_solvers � ��.):

cmake --build . --config Debug -j 8
cmake --install . --config Debug



���������� Release

������ ���� ����������� � ������������ � ������� CMAKE:
- Eigen Release
- Gtest Release
- fixed_solvers Release

����������� ��. fixed_solvers/msvc_cmake/readme.txt
(���� � ����������� ��������� ���������, ������� ���������� ��� ������ ������� ����������, �� �� ���� ������������ � �����������������)

��� ���� ��������� ������������ ���������� ������� ������� (��������� �� pde_solvers/msvc_cmake)

cmake .. -G "Visual Studio 17 2022" -DCMAKE_PREFIX_PATH="C:/Program Files (x86)/googletest-distribution"  -DCMAKE_MSVC_RUNTIME_LIBRARY=MultiThreadedDLL -DCMAKE_BUILD_TYPE=Release -DPDE_SOLVERS_BUILD_TESTS=ON

� ����� pde_solvers/msvc_cmake ����� ������� ������� MSVC (.sln). ��� ���������, �������� � ��� ��� ������.

��� ��������� � ����� ������������� � ������ ��������� ������� (graph_solvers � ��.):

cmake --build . --config Release -j 8
cmake --install . --config Release

