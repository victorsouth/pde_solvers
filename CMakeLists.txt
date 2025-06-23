cmake_minimum_required(VERSION 3.10)
# header only library
# for install, dependencies and correct include path
project(pde_solvers)

set(CMAKE_CXX_STANDARD 17)#23?
set(CMAKE_CXX_STANDARD_REQUIRED True)

if(NOT TARGET fixed_solvers::fixed_solvers)
find_package(fixed_solvers)
endif()

set(HEADERS
    pde_solvers/pde_solvers.h  pde_solvers/pipe.h pde_solvers/timeseries.h
    )
set(HEADERS_CORE
    pde_solvers/core/differential_equation.h  pde_solvers/core/profile_structures.h  pde_solvers/core/ring_buffer.h
    )
set(HEADERS_PIPE
    pde_solvers/pipe/oil.h
    pde_solvers/pipe/pipe_advection_solver.h
    pde_solvers/pipe/pipe_hydraulic_pde.h
    pde_solvers/pipe/pipe_advection_pde.h
    pde_solvers/pipe/pipe_hydraulic_computations.h
    pde_solvers/pipe/pipe_hydraulic_struct.h
    pde_solvers/pipe/pipe_profile_utils.h
    pde_solvers/pipe/heat_transfer.h          
    pde_solvers/pipe/pipe_dynamic_soil_multizone.h  
    pde_solvers/pipe/pipe_heat_computations.h
    pde_solvers/pipe/pipe_heat_pde.h
    pde_solvers/pipe/pipe_heat_struct.h
    pde_solvers/pipe/pipe_heat_util.h
    pde_solvers/pipe/pipe_dynamic_soil.h      
)
set(HEADERS_SOLVERS
    pde_solvers/solvers/diffusion_solver.h
    pde_solvers/solvers/godunov_solver.h
    pde_solvers/solvers/moc_solver.h
    pde_solvers/solvers/ode_solver.h
    pde_solvers/solvers/quick_solver.h
)
file(GLOB HEADERS_TASKS pde_solvers/tasks/* )

set(HEADERS_TIME
pde_solvers/timeseries/csv_readers.h  pde_solvers/timeseries/timeseries_helpers.h  pde_solvers/timeseries/vector_timeseries.h
pde_solvers/timeseries/synthetic_timeseries.h
)


if("3.19.0" VERSION_LESS ${CMAKE_VERSION})
    add_library(${PROJECT_NAME} INTERFACE ${HEADERS} ${HEADERS_CORE} ${HEADERS_PIPE} ${HEADERS_SOLVERS} ${HEADERS_TIME} ${HEADERS_TASKS})
else()
    add_library(${PROJECT_NAME} INTERFACE)
endif()
target_link_libraries(${PROJECT_NAME} INTERFACE fixed_solvers::fixed_solvers)
target_include_directories(${PROJECT_NAME}
    INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/pde_solvers>
    $<INSTALL_INTERFACE:include>
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
configure_package_config_file(${CMAKE_CURRENT_SOURCE_DIR}/Config.cmake.in
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
install(FILES ${HEADERS_TIME} DESTINATION include/pde_solvers/timeseries)
install(FILES ${HEADERS_TASKS} DESTINATION include/pde_solvers/tasks)
endif()

option(PDE_SOLVERS_BUILD_TESTS "" ON)

if(PDE_SOLVERS_BUILD_TESTS)

    find_package(Threads)
    find_package(GTest REQUIRED)
    set(TESTS_HEADERS
        testing/test_advection_moc_solver.h
        testing/test_diffusion.h
        testing/test_moc.h
        testing/test_quick.h
        testing/test_static_pipe_solver.h
        testing/test_timeseries.h
    )
    add_executable(pde_tests testing/test_main.cpp ${TESTS_HEADERS})
    target_link_libraries(pde_tests pde_solvers::pde_solvers GTest::gtest)


    set(RESEARCH_HEADERS
        research/2023-12-diffusion-of-advection/diffusion_of_advection.h
        research/2024-02-quick-with-quasistationary-model/quick_with_quasistationary_model.h
        research/2024-08-quasistationary-with-real-data/quasistationary_with_real_data.h
        research/2024-10-ident-quasistatic-isothermal/ident_quasistatic_isothermal.h
        research/2025-01-ident-quasistatic-nonisothermal/ident_quasistatic_nonisothermal.h
        research/2025-01-ident-quasistatic-nonisothermal/nonisothermal_quasistationary_with_real_data.h
        research/2025-04-calc-mass-on-isothermal-quasistatic/calc_mass_on_isothermal_quasistatic.h
    )
    add_executable(pde_research research/research_main.cpp ${RESEARCH_HEADERS})
    target_link_libraries(pde_research pde_solvers::pde_solvers GTest::gtest)



endif()
