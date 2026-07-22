#include <pde_solvers/pde_solvers.h>

#define GTEST_BREAK_ON_FAILURE 1
#define GTEST_CATCH_EXCEPTIONS 0
#define GTEST_HAS_SEH 0
#define _VARIADIC_MAX 10 /* for gtest */
#include "gtest/gtest.h"


#include <iostream>
#include <fstream>
#include <filesystem>


#include <time.h>
#include <algorithm>

/// @brief Путь к каталогу с исследованиями
inline std::filesystem::path research_source_dir()
{
    return std::filesystem::path(__FILE__).parent_path();
}

/// @brief Путь к каталогу с результатами исследований
inline std::filesystem::path research_out_dir()
{
    return std::filesystem::path(research_source_dir().string() + "_out");
}

inline std::string prepare_research_folder()
{
    auto test_info = ::testing::UnitTest::GetInstance()->current_test_info();
    std::string research_name = std::string(test_info->test_case_name());
    std::string case_name = std::string(test_info->name());

    std::filesystem::path p = research_out_dir() / research_name / case_name;
    std::filesystem::create_directories(p);
    return p.string() + "/";
}

inline std::string prepare_research_folder_for_qsm_model2()
{
    std::string path = prepare_research_folder();
    for (const auto& entry : std::filesystem::directory_iterator(path)) {
        std::filesystem::remove_all(entry.path());
    }
    return path;
}

/// @brief Возвращает путь к папке, где лежат данные по данной трубе
/// В папке лежат как параметры самой трубы (профиль), так и временные ряды (.csv)
/// @param pipe_name Стандартное название трубы (cold_lu, hot_lu, condensate_lu)
inline std::string get_pipe_data_path(std::string pipe_name = "")
{
    std::filesystem::path p = research_out_dir() / "data";
    if (!pipe_name.empty()) {
        p = p / pipe_name;
    }
    return p.string() + "/";
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
