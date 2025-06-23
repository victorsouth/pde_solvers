
#define GTEST_BREAK_ON_FAILURE 1
#define GTEST_CATCH_EXCEPTIONS 0
#define GTEST_HAS_SEH 0
#define _VARIADIC_MAX 10 /* for gtest */
#include "gtest/gtest.h"

#include <fixed/fixed.h>
#include <pde_solvers/pde_solvers.h>

#include <iostream>
#include <fstream>
#include <filesystem>


#include <time.h>
#include <algorithm>



inline std::string prepare_research_folder()
{
    auto test_info = ::testing::UnitTest::GetInstance()->current_test_info();
    std::string research_name = std::string(test_info->test_case_name());
    std::string case_name = std::string(test_info->name());

    std::string path = std::string("../research_out/") + 
        research_name + "/" + case_name + "/";
    std::filesystem::create_directories(path);
    return path;
}

using namespace pde_solvers;


#include "../research/2023-12-diffusion-of-advection/diffusion_of_advection.h"
#include "../research/2024-02-quick-with-quasistationary-model/quick_with_quasistationary_model.h"
#include "../research/2024-08-quasistationary-with-real-data/quasistationary_with_real_data.h"
#include "../research/2024-10-ident-quasistatic-isothermal/ident_quasistatic_isothermal.h"
#include "../research/2025-04-calc-mass-on-isothermal-quasistatic/calc_mass_on_isothermal_quasistatic.h"
#include "../research/2025-01-ident-quasistatic-nonisothermal/nonisothermal_quasistationary_with_real_data.h"
#include "../research/2025-01-ident-quasistatic-nonisothermal/ident_quasistatic_nonisothermal.h"


int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
#if defined(_WIN32) && !defined(__MINGW32__)
    std::wcout.imbue(std::locale("rus_rus.866"));
#endif
int res= RUN_ALL_TESTS();
    return res;
}
