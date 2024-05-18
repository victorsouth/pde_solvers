#pragma once
//#include <cmath>
#include <pde_solvers/pde_solvers.h>

/// @brief Генерация значейний высотных отметок профиля
/// @param coords Сетка координат генерируемого профиля
inline vector<double> build_heights(const vector<double>& coords)
{
	size_t point_cnt = coords.size();
	vector<double> res;
	for (size_t index = 0; index < point_cnt; index++)
		// 50 * cos(x/10) + 10 * sin(x + 150)
		res.push_back(50 * cos(coords[index] / 10) + 10 * sin(coords[index]) + 150);
	return res;
}

/// @brief Генерация профиля трубопровода с непостоянным шагом по детерменированной функции
/// @param length Длина генерируемого профиля
/// @param step Примерный шаг сетки по координате
inline vector<vector<double>> build_profile(const double length, const double step)
{
	vector<vector<double>> coord_heights(2);
	coord_heights[0] = { 0 };
	size_t point_count = static_cast<size_t>(length / step);

	for (size_t i = 1; i <= point_count; i++)
	{
		coord_heights[0].push_back(coord_heights[0][i - 1] + step + 10 * (i % 2));
	}

	if (coord_heights[0].back() != length)
		coord_heights[0].back() = length;

	coord_heights[1] = build_heights(coord_heights[0]);

	return coord_heights;
}

/// @brief Вывод профиля в файл
/// @param coord_heights вектор, первый элемент которого сетка координат, а второй - соответствующие высотные отметки
/// @param filename путь и имя файла
inline void write_profile(const vector<vector<double>>& coord_heights,const string filename)
{
	std::ofstream output_file;
	size_t profCount = coord_heights[0].size();
	output_file.open(filename + ".csv");
	output_file << "plug,km,height" << std::endl;
	for (size_t i = 0; i < profCount; i++)
	{
		output_file << "0," << coord_heights[0][i] << ',' << coord_heights[1][i] << endl;
	}
}

/// @brief Генерация профиля с переменным шагом по координате и примнение его 
/// в качестве исходного профиля для создания профиля с постоянным шагом
/// @param L Длина трубопровода
/// @param dx Примерный шаг исходного профиля
/// @param desired_dx Желаемый шаг нового профиля
/// @param write_flg Флаг записи профилей в файл
/// @param file_name Название файла для вывода
/// @return Возвращает пару - исходный профиль с переменным шагом и новый с постоянным шагом
inline std::pair<vector<vector<double>>, PipeProfile> different_source_prof(const double L, const double dx, const double desired_dx, 
	const bool write_flg = false, const string file_name = "file_name")
{
	// Генерируем исходный профиль по детерменированной функции с переменным шагом
	vector<vector<double>> source_profile = build_profile(L, dx);

	// Создаём новый профиль с постоянным шагом
	PipeProfile new_prof = pipe_profile_uniform::get_uniform_profile(source_profile, desired_dx);

	// Записываем в файлы при необходимости
	if (write_flg)
	{
		write_profile(source_profile, file_name);
		write_profile({ new_prof.coordinates, new_prof.heights }, file_name + "_new");
	}
	

	return { source_profile, new_prof };
}

/// @brief Создание профиля с шагом равным желаемому
TEST(CreatingProfile, CreateWithEqStep)
{
	double L = 1000;
	double dx = 50;

	double desired_dx = 200;

	auto [source_prof, new_prof] = different_source_prof(L, dx, desired_dx);

	// Проверка на равенство желаемого шага и реального в новом профиле
	ASSERT_NEAR(desired_dx, new_prof.coordinates[1] - new_prof.coordinates[0], 1e-8);
}

/// @brief Создание профиля с шагом большим чем желаемый
TEST(CreatingProfile, CreateWithDiffStep)
{
	double L = 1000;
	double dx = 50;

	double desired_dx = 300;

	auto [source_prof, new_prof] = different_source_prof(L, dx, desired_dx);

	// Проверка на неравенство желаемого шага и реального в новом профиле
	ASSERT_LT(0, new_prof.coordinates[1] - new_prof.coordinates[0] - desired_dx);
}

/// @brief Создание профиля, когда желаемый шаг больше длины исходного профиля
TEST(CreatingProfile, CreateFromShortPipe)
{
	double L = 200;
	double dx = 50;

	double desired_dx = 300;

	auto [source_prof, new_prof] = different_source_prof(L, dx, desired_dx);

	// Проверка на неравенство желаемого шага и реального в новом профиле
	ASSERT_NEAR(desired_dx, new_prof.getLength(), 1e-8);
}

/// @brief Создание профиля из различных видов исходного профиля
/// Плотный профиль - профиль, шаг сетки которого меньше желаемого профиля
/// Разреженный профиль - профиль, шаг сетки которого больше желаемого профиля
/// Третий случай - когда шаг исходного профиля и желаемый шаг совпадают
TEST(CreatingProfile, CreateFromDifferentSourceProfile)
{
	string path = prepare_test_folder();

	double L = 1000;
	double dx = 50;
	double desired_dx;
	
	// Создание профиля в случае, когда исходный профиль плотный
	desired_dx = 100;
	different_source_prof(L, dx, desired_dx, true, path + "a_from_tight_prof");
	
	// Создание профиля в случае, когда желаемый шаг совпадает с шагом исходного профиля
	desired_dx = 50;
	different_source_prof(L, dx, desired_dx, true, path + "b_from_norm_prof");

	// Создание профиля в случае, когда исходный профиль разреженный
	desired_dx = 10;
	different_source_prof(L, dx, desired_dx, true, path + "c_from_rare_prof");

}

/// @brief Создание профиля из различных видов исходного профиля
/// Плотный профиль - профиль, шаг сетки которого меньше желаемого профиля
/// Разреженный профиль - профиль, шаг сетки которого больше желаемого профиля
/// Третий случай - когда шаг исходного профиля и желаемый шаг совпадают
TEST(CreatingProfile, UseCaseFromFile)
{
	// Указываем имя файла и желаемый шаг новой сетки
	string file_name = "coord_heights.csv";
	double desired_dx = 100;

	// Создаём новый профиль с постоянным шагом
	PipeProfile new_prof = pipe_profile_uniform::get_uniform_profile_from_csv(desired_dx, file_name);
}

