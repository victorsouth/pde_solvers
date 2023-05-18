# header only library
# for install, dependencies and correct include path
project(pde_solvers)
cmake_minimum_required(VERSION 3.10)

set(CMAKE_CXX_STANDARD 17)#23?
set(CMAKE_CXX_STANDARD_REQUIRED True)

find_package(fixed_solvers)

set(HEADERS
    pde_solvers/pde_solvers.h  pde_solvers/pipe.h
    )
set(HEADERS_CORE
    pde_solvers/core/custom_buffer.h  pde_solvers/core/differential_equation.h  pde_solvers/core/profile_structures.h
    )
set(HEADERS_PIPE
    pde_solvers/pipe/oil.h  pde_solvers/pipe/pipe_advection_pde.h  pde_solvers/pipe/pipe_hydraulic_computations.h  pde_solvers/pipe/pipe_hydraulic_pde.h  pde_solvers/pipe/pipe_hydraulic_struct.h
    )
set(HEADERS_SOLVERS
    pde_solvers/solvers/godunov_solver.h  pde_solvers/solvers/moc_solver.h  pde_solvers/solvers/ode_solver.h
)

if("3.19.0" VERSION_LESS ${CMAKE_VERSION})
    add_library(${PROJECT_NAME} INTERFACE ${HEADERS} ${HEADERS_LS})
else()
    add_library(${PROJECT_NAME} INTERFACE)
endif()
target_link_libraries(${PROJECT_NAME} INTERFACE fixed_solvers::fixed_solvers)
target_include_directories(${PROJECT_NAME}
    INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:include/pde_solvers>
)
add_library(pde_solvers::pde_solvers ALIAS pde_solvers)

option(PDE_SOLVERS_INSTALL "" ON)

if(PDE_SOLVERS_INSTALL)
    include(CMakePackageConfigHelpers)
    include(GNUInstallDirs)
    install(TARGETS pde_solvers
        DESTINATION ${CMAKE_INSTALL_LIBDIR}
        EXPORT pde_solversTargets
    )
install(EXPORT pde_solversTargets NAMESPACE pde_solvers::
    FILE pde_solversTargets.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/pde_solvers
)
# generate the config file that is includes the exports
configure_package_config_file(${CMAKE_SOURCE_DIR}/Config.cmake.in
    "${CMAKE_CURRENT_BINARY_DIR}/pde_solversConfig.cmake"
    INSTALL_DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/pde_solvers"
    NO_SET_AND_CHECK_MACRO
    NO_CHECK_REQUIRED_COMPONENTS_MACRO
)
install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/pde_solversConfig.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/pde_solvers
)

install(FILES ${HEADERS} DESTINATION include/pde_solvers/)
install(FILES ${HEADERS_CORE} DESTINATION include/pde_solvers/core)
install(FILES ${HEADERS_PIPE} DESTINATION include/pde_solvers/pipe)
install(FILES ${HEADERS_SOLVERS} DESTINATION include/pde_solvers/solvers)
endif()

option(PDE_SOLVERS_BUILD_TESTS "" OF)

if(PDE_SOLVERS_BUILD_TESTS)

find_package(Threads)
find_package(GTest REQUIRED)

add_executable(pde_tests testing/test_main.cpp)
target_link_libraries(pde_tests pde_solvers::pde_solvers GTest::gtest)

endif()
