���������� Debug
������ ���� ���������� � ����������� Gtest Debug

cmake .. -G "Visual Studio 17 2022" -DCMAKE_PREFIX_PATH="C:/Program Files (x86)/googletest-distribution"  -DCMAKE_MSVC_RUNTIME_LIBRARY=MultiThreadedDebugDLL -DCMAKE_BUILD_TYPE=Debug -DPDE_SOLVERS_BUILD_TESTS=ON