#pragma once

/// @brief Генерация значений высотных отметок профиля
/// по определённой гармонической функции: 50 * cos(x/10) + 10 * sin(x + 150)
/// @param coords Сетка координат генерируемого профиля
inline vector<double> build_heights_on_harmonic_function(const vector<double>& coords)
{
	size_t point_cnt = coords.size();
	vector<double> res;
	for (size_t index = 0; index < point_cnt; index++)
		res.push_back(50 * cos(coords[index] / 10) + 10 * sin(coords[index]) + 150);
	return res;
}

/// @brief Генерация профиля трубопровода с непостоянным шагом по детерменированной функции
/// @param length Длина генерируемого профиля
/// @param step Примерный шаг сетки по координате
inline PipeProfile build_profile_with_diff_step(const double length, const double step)
{
	PipeProfile source_prof;
	source_prof.coordinates = { 0 };
	size_t point_count = static_cast<size_t>(length / step);

	for (size_t i = 1; i <= point_count; i++)
	{
		// К каждому нечётному шагу прибавляется 10
		source_prof.coordinates.push_back(source_prof.coordinates[i - 1] + step + step * 0.05 * (i % 2));
	}

	if (source_prof.coordinates.back() != length)
		source_prof.coordinates.back() = length;
	// Построение профиля высоток
	source_prof.heights = build_heights_on_harmonic_function(source_prof.coordinates);
	// Несущая способность
	source_prof.capacity = vector<double>(source_prof.getPointCount(), 10e6);

	return source_prof;
}

/// @brief Вывод профиля в файл
/// @param coord_heights вектор, первый элемент которого сетка координат, а второй - соответствующие высотные отметки
/// @param filename путь и имя файла
inline void write_profile(const PipeProfile& profile,const string filename)
{
	std::ofstream output_file;
	size_t profCount = profile.getPointCount();
	output_file.open(filename + ".csv");
	output_file << "plug,km,height" << std::endl;
	for (size_t i = 0; i < profCount; i++)
	{
		output_file << "0," << profile.coordinates[i] << ',' << profile.heights[i] << endl;
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
inline std::pair<PipeProfile, PipeProfile> create_uniform_prof_from_generated_prof(const double L, const double dx, const double desired_dx, const string file_name = "")
{
	// Генерируем исходный профиль по детерменированной функции с переменным шагом
	PipeProfile source_profile = build_profile_with_diff_step(L, dx);

	// Создаём новый профиль с постоянным шагом
	PipeProfile new_prof = create_uniform_profile(source_profile, desired_dx);

	// Записываем в файлы при необходимости
	if (!file_name.empty())
	{
		write_profile(source_profile, file_name);
		write_profile(new_prof, file_name + "_new");
	}
	
	return { source_profile, new_prof };
}

/// @brief Проверка случая, когда шаг итогового профиля совпадает с желаемым
TEST(UniformProfile, CheckResStepEqualDesStep)
{
	// Длина исходного профиля
	double L = 1000;
	// Шаг исходного профиля
	double dx = 50;

	// Шаг итогового профиля совпадает с желаемым, если длина трубы делится на желаемый шаг нацело
	double desired_dx = 200;

	auto [source_prof, new_prof] = create_uniform_prof_from_generated_prof(L, dx, desired_dx);

	// Проверка на равенство желаемого шага и реального в новом профиле
	ASSERT_NEAR(desired_dx, new_prof.coordinates[1] - new_prof.coordinates[0], 1e-8);
}

/// @brief Проверка случая, когда шаг итогового профиля больше желаемого
TEST(UniformProfile, CheckResStepGreaterDesStep)
{
	// Длина исходного профиля
	double L = 1000;
	// Шаг исходного профиля
	double dx = 50;

	// Шаг итогового профиля не совпадает с желаемым, если длина трубы не делится на желаемый шаг нацело
	double desired_dx = 300;

	auto [source_prof, new_prof] = create_uniform_prof_from_generated_prof(L, dx, desired_dx);

	// Проверка на неравенство желаемого шага и реального в новом профиле
	ASSERT_LT(0, new_prof.coordinates[1] - new_prof.coordinates[0] - desired_dx);
}

/// @brief Создание профиля, когда желаемый шаг больше длины исходного профиля
TEST(UniformProfile, HandlesShortSourceProf)
{
	// Длина исходного профиля
	double L = 200;
	// Шаг исходного профиля
	double dx = 50;
	//Желаемый шаг
	double desired_dx = 300;

	auto [source_prof, new_prof] = create_uniform_prof_from_generated_prof(L, dx, desired_dx);

	// Проверка на неравенство желаемого шага и реального в новом профиле
	ASSERT_NEAR(desired_dx, new_prof.getLength(), 1e-8);
}

/// @brief Сравнение итоговых профилей при различных видах исходного профиля
/// Плотный профиль - профиль, шаг сетки которого меньше желаемого профиля
/// Разреженный профиль - профиль, шаг сетки которого больше желаемого профиля
/// Третий случай - когда шаг исходного профиля и желаемый шаг совпадают
TEST(UniformProfile, CompareResProfFromDiffSourceProf)
{
	string path = prepare_test_folder();

	// Длина исходного профиля
	double L = 1000;
	// Шаг исходного профиля
	double dx = 50;
	//Желаемый шаг
	double desired_dx;
	
	// Создание профиля в случае, когда исходный профиль плотный
	desired_dx = 100;
	create_uniform_prof_from_generated_prof(L, dx, desired_dx, path + "a_from_tight_prof");
	
	// Создание профиля в случае, когда желаемый шаг совпадает с шагом исходного профиля
	desired_dx = 50;
	create_uniform_prof_from_generated_prof(L, dx, desired_dx, path + "b_from_norm_prof");

	// Создание профиля в случае, когда исходный профиль разреженный
	desired_dx = 10;
	create_uniform_prof_from_generated_prof(L, dx, desired_dx, path + "c_from_rare_prof");

}

/// @brief Пример созддания профиля, когда исходный профиль считывается из файла
TEST(UniformProfile, UseCaseSourceProfFromFile)
{
	// Указываем имя файла и желаемый шаг новой сетки
	string file_name = "coord_heights.csv";
	//Желаемый шаг
	double desired_dx = 100;

	// Создаём новый профиль с постоянным шагом
	PipeProfile new_prof = pipe_profile_uniform::get_uniform_profile_from_csv(desired_dx, file_name);
}

