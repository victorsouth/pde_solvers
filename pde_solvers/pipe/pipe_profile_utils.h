
#include <fstream>

namespace pde_solvers {
;

using std::max;
using std::string;


/// @brief Чтение координат и соответствующих высоток из файла csv
/// первая строка с названием колонок, следующие строки в формате km;m
/// @param filename Путь к файлу
/// @return Вектор векторов - координаты и высотки
inline vector<vector<double>> read_coordinates_and_heights_file(const std::string filename)
{
    std::ifstream file(filename);
    std::string line;
    char delimiter = ';';
    vector<double> coords;
    vector<double> heights;

    // Первую строку c названиями колонок пропускаем
    std::getline(file, line);
    while (std::getline(file, line))
    {
        std::stringstream stream(line);
        std::string coord;
        std::string height;

        std::getline(stream, coord, delimiter);
        std::getline(stream, height, delimiter);
        coords.push_back(stod(coord) * 1000);
        heights.push_back(stod(height));
    }

    file.close();
    return { coords, heights };
}


/// @brief Класс для создание профиля с желаемым постоянным шагом по координате
/// из исходного профиля, который в общем случае имеет непостоянный шаг сетки
class pipe_profile_uniform
{
protected: // обработка координаток
    /// @brief Создание сетки нового профиля с заданным шагом по координате
    /// @param segment_len Длина шага по координате
    /// @param point_cnt Количество точек профиля - на 1 больше чем количество сегментов
    /// @param offset Координата начала профиля
    /// @return Новая координатная сетка с заданным шагом
    static vector<double> generate_uniform_grid(const double segment_len, const size_t point_cnt, const double offset = 0)
    {
        vector<double> grid(point_cnt);
        for (size_t index = 0; index < point_cnt; index++)
        {
            grid[index] = offset + segment_len * index;
        }

        return grid;
    }

    static vector<double> generate_uniform_grid(double initial_coordinate,
        double pipe_length, double desired_uniform_segment)
    {
        // В большинстве случаев длина сегмента естественным образом вырастет по отношению к желаемому
        // В случае короткой трубы, меньшей desired_uniform_segment, обеспечивается как минимум один 
        // сегмент длиной desired_uniform_segment
        size_t segment_count = max<size_t>(1, static_cast<size_t>(pipe_length / desired_uniform_segment));
        desired_uniform_segment = max<double>(desired_uniform_segment, pipe_length / segment_count);

        return generate_uniform_grid(desired_uniform_segment, segment_count + 1, initial_coordinate);
    }

protected:  // обработка коротких труб
    /// @brief Обработка случая короткой трубы исходного профиля по отношению к желаемому шагу по координате
    /// @param new_coordinates Координатная сетка нового профиля
    /// @param exact_coordinates Координатная сетка исходного профиля
    /// @param _exact_parameters Профиль из исходных данных
    static void extend_values(const vector<double>& new_coordinates,
        const vector<double>& exact_coordinates, vector<double>* _exact_parameters)
    {
        vector<double>& exact_parameters = *_exact_parameters;

        // если НАЧАЛЬНАЯ точка new_coordinates меньше НАЧАЛЬНОЙ точки exact_coordinates, 
        // экстраполировать exact_coordinates ВЛЕВО
        if (exact_coordinates.front() - new_coordinates.front() > 1e-8) {
            //exact_coordinates.insert(exact_coordinates.begin(), new_coordinates.front() + offset);
            exact_parameters.insert(exact_parameters.begin(), exact_parameters.front());
        }

        // если КОНЕЧНАЯ точка new_coordinates меньше конечной. точки exact_coordinates, 
        // экстраполировать exact_coordinates ВПРАВО
        if (new_coordinates.back() - exact_coordinates.back() > 1e-8) {
            //exact_coordinates.push_back(new_coordinates.back() + offset);
            exact_parameters.push_back(exact_parameters.back());
        }
    }

    /// @brief Обработка случая короткой трубы исходного профиля по отношению к желаемому шагу по координате
    /// @param new_coordinates Координатная сетка нового профиля
    /// @param exact_coordinates Координатная сетка исходного профиля
    static void extend_arguments(const vector<double>& new_coordinates,
        vector<double>* _exact_coordinates)
    {
        vector<double>& exact_coordinates = *_exact_coordinates;

        // если НАЧАЛЬНАЯ точка new_coordinates меньше НАЧАЛЬНОЙ точки exact_coordinates, 
        // экстраполировать exact_coordinates ВЛЕВО
        if (exact_coordinates.front() - new_coordinates.front() > 1e-8) {
            exact_coordinates.insert(exact_coordinates.begin(), new_coordinates.front());
        }

        // если КОНЕЧНАЯ точка new_coordinates больше КОНЕЧНОЙ точки exact_coordinates, 
        // экстраполировать exact_coordinates ВПРАВО
        if (new_coordinates.back() - exact_coordinates.back() > 1e-8) {
            exact_coordinates.push_back(new_coordinates.back());
        }
    }
private: // обработка высоток и несущей
    /// @brief Увеличение количества точек в случае разреженного исходного профиля
    /// @param source_prof Исходный профиль
    /// @param max_segment Длина желаемого сегмента
    /// @return Более плотный профиль в случае разреженного исходного профиля
    /// Исходный профиль в обратном случае
    static pipe_profile_t subdivide_irregular_profile(pipe_profile_t source_prof, const double max_segment)
    {
        // Если все dl меньше чем max_segment, то профиль останется неизменным
        // Там, где dl больше чем max_segment, профиль будет дополнен 
        // с шагом не меньшим чем max_segment
        pipe_profile_t new_prof;
        for (size_t segment = 0; segment < source_prof.get_point_count() - 1; segment++)
        {
            double dl = source_prof.coordinates[segment + 1] - source_prof.coordinates[segment];
            double dh = source_prof.heights[segment + 1] - source_prof.heights[segment];
            double dcapacity = source_prof.capacity[segment + 1] - source_prof.capacity[segment];

            size_t divide_cnt = static_cast<size_t>(ceil(dl / max_segment) + 1e-8);

            // Заполнение промежутков новыми точками
            for (size_t offset = 0; offset < divide_cnt; offset++)
            {
                new_prof.coordinates.push_back(source_prof.coordinates[segment] + dl * offset / divide_cnt);
                new_prof.heights.push_back(source_prof.heights[segment] + dh * offset / divide_cnt);
                new_prof.capacity.push_back(source_prof.capacity[segment] + dcapacity * offset / divide_cnt);
            }
        }

        new_prof.coordinates.push_back(source_prof.coordinates.back());
        new_prof.heights.push_back(source_prof.heights.back());
        new_prof.capacity.push_back(source_prof.capacity.back());

        return new_prof;
    }

    /// @brief Определение границ областей притяжения точек нового профиля
    /// @param uniform_coordinates Координатная сетка нового профиля
    /// @return Координаты границ областей притяжения точек нового профиля
    static vector<double> generate_influence_segments(const vector<double>& uniform_coordinates)
    {
        if (uniform_coordinates.size() <= 1) {
            throw std::logic_error("generate_influence_segments(): uniform_coordinates.size() <= 1 ");
        }
        else if (uniform_coordinates.size() == 2) {
            return vector<double>{ uniform_coordinates.front() }; // здесь 0 сегментов
        }

        double segment_length = uniform_coordinates[1] - uniform_coordinates[0];

        // Не хотим менять крайние точки
        // Вторая сначала и предпоследняя с конца точка имеет область притяжения 1.5 сегмента
        // Все остальные имеют область притяжения по полсегмента в каждую сторону, т.е. 1.0 сегмент

        vector<double> result;
        result.push_back(uniform_coordinates.front());
        for (size_t index = 1; index < uniform_coordinates.size() - 2; ++index) {
            result.push_back(uniform_coordinates[index] + 0.5 * segment_length);
        }
        result.push_back(uniform_coordinates.back());
        return result;
    }

    /// @brief Определние индексов границ областей притяжения на исходном профиле
    /// @param exact_coordinates Исходный профиль
    /// @param sample_coordinates Границы областей притяжения 
    /// @return Индексы границ областей притяжения на исходном профиле
    static vector<size_t> get_decimated_coordinates(const vector<double>& exact_coordinates,
        const vector<double>& sample_coordinates)
    {
        if (exact_coordinates.size() == 0)
            throw std::logic_error("get_decimated_coordinates(): empty coordinates");

        size_t sample_index = 0;

        vector<size_t> indices;
        indices.reserve(sample_coordinates.size());
        for (size_t index = 0; index < exact_coordinates.size(); ++index) {
            while (exact_coordinates[index] > sample_coordinates[sample_index] - 1e-8)
            {
                indices.push_back(index);
                sample_index++;
                if (sample_index >= sample_coordinates.size())
                    return indices;
            }
        }

        return indices;
    }

    /// @brief Определение высотных отметок соответствующим координатам нового профиля
    /// @tparam Function Принцип выбора соответствующей высотной отметни на области притяжения
    /// @param indices Индексы границ областей притяжения на исходном профиле
    /// @param values Высотные отметки исходного профиля
    /// @param function Принцип выбора соответствующей высотной отметни на области притяжения
    /// @return Вектор соответствующих высоток для сетки нового профиля
    template <typename Function>
    static vector<double> execute_function_on_profile_segment(const vector<size_t>& indices,
        const vector<double>& values, Function& function)
    {
        if (indices.empty()) {
            throw std::logic_error("execute_function_on_profile_segment() empty indices");
        }
        if (indices.size() == 1) {
            return vector < double > {values.front(), values.back()};
        }

        vector<double> result(indices.size() + 1);
        result.front() = values.front();
        result.back() = values.back();

        // первая и последняя точка из values всегда попадают в результат
        // их надо исключить из подмножеств
        // последняя точка исключена из диапазона ввиду соглашений на вызов функтора (STL-like)
        // первую точку надо исключить специальной проверкой

        for (size_t segment = 0; segment < indices.size() - 1; ++segment) {
            size_t index_from = indices[segment];
            if (index_from == 0)
                index_from = 1; // исключаем нулевую точку из нулевого сегмента

            size_t index_to = indices[segment + 1];
            // функция вызывается в стандарте STL, второй итератор указывает на элемент, следующий за последним
            result[segment + 1] = *function(values.begin() + index_from, values.begin() + index_to);
        }

        return result;
    }
    
    /// @brief Формирует профили высоток и несущей для заднного равномерного профиля координат
    /// Подразумевает, что труба source_profile не короче трубы uniform_coordinates
    /// @return пара ["равномерные" высотки; "равномерная" несущая]
    static pair<vector<double>, vector<double>> create_uniform_height_and_capacity(pipe_profile_t source_profile, 
        const vector<double>& uniform_coordinates) 
    {
        // Увеличение плотности в случае разреженного профиля
        double segment_length = uniform_coordinates[1] - uniform_coordinates[0];
        source_profile = subdivide_irregular_profile(source_profile, segment_length / 2);

        // Определение границ областей притяжения точек нового профиля
        vector<double> influence_segments = generate_influence_segments(uniform_coordinates);

        // получения индексов начада областей притяжения на исходном профиле
        vector<size_t> influence_segment_indices =
            get_decimated_coordinates(source_profile.coordinates, influence_segments);

        // Определение высотных отметок на новом профиле по максимальному значению в соответствующей области притяжения
        vector<double> uniform_heights = execute_function_on_profile_segment(influence_segment_indices,
            source_profile.heights, std::max_element<vector<double>::const_iterator>);

        vector<double> uniform_capacity = execute_function_on_profile_segment(influence_segment_indices,
            source_profile.capacity, std::min_element<vector<double>::const_iterator>);

        // Проверка на соответствие количества точек координат количеству высотных отметок
        if (uniform_coordinates.size() != uniform_heights.size()) {
            throw std::logic_error("pipeline_profile_t::create_uniform_profile(): coordinates.size() != heights.size()");
        }

        return std::make_pair(std::move(uniform_heights), std::move(uniform_capacity));
    }

public:

    /// @brief Создание профиля с постоянным шагом по координате
    /// @param source_profile Исходный профиль с непостоянным 
    /// в общем случае шагом по координате
    /// @param desired_uniform_segment Желаемый шаг по координате для новой сетки
    /// @return Профиль с постоянным близким к желаемому шагом по координате
    static pipe_profile_t create_uniform_profile(const pipe_profile_t& source_profile, double desired_uniform_segment)
    {
        // Создание итогого профиля
        pipe_profile_t uniform_profile;
        // Создание сетки с постоянным шагом по координате
        uniform_profile.coordinates = generate_uniform_grid(source_profile.coordinates.front(),
            source_profile.get_length(), desired_uniform_segment);

        // Обработка случая короткой трубы. 
        // Может получиться, если source_profile короче, чем desired_uniform_segment
        // В этом случае удлиняем исходную трубу
        pipe_profile_t extended_source_profile = source_profile;
        extend_values(uniform_profile.coordinates, extended_source_profile.coordinates, 
            &extended_source_profile.heights); 
        extend_values(uniform_profile.coordinates, extended_source_profile.coordinates,
            &extended_source_profile.capacity);
        extend_arguments(uniform_profile.coordinates,
            &extended_source_profile.coordinates);

        // Подготовка профиля трассы и несущей способности 
        // на основе (возможно удлиненной) исходной трубы и подготовленного равномерного профиля
        std::tie(uniform_profile.heights, uniform_profile.capacity) =
            create_uniform_height_and_capacity(extended_source_profile, uniform_profile.coordinates);

        return uniform_profile;

    }

    /// @brief Создание профиля с постоянным шагом
    /// @param coord_heights Вектор двух векторов - координат и соответствующих им высоток
    /// @param desired_segment Желаемый постоянный шаг по координате новой сетки профиля 
    /// @return Профиль с постоянным шагом по координате
    static pipe_profile_t get_uniform_profile(const vector<vector<double>>& coord_heights,
        double desired_segment, double capacity_value = 10e6)
    {
        pipe_profile_t source_prof;
        source_prof.coordinates = coord_heights[0];
        source_prof.heights = coord_heights[1];

        // Заглушка для несущей
        source_prof.capacity = vector<double>(source_prof.get_point_count(), capacity_value);

        return create_uniform_profile(source_prof, desired_segment);
    }

    /// @brief Создание профиля с постоянным шагом по данным профиля из файла csv
    /// @param desired_segment Желаемый постоянный шаг новой сетки профиля
    /// @param filename Название файла с координатами и соответствующими высотками
    /// @return Профиль с постоянным шагом по координате
    static pipe_profile_t get_uniform_profile_from_csv(const double desired_segment, const string& filename)
    {
        vector<vector<double>> km_heights = read_coordinates_and_heights_file(filename);
        return get_uniform_profile(km_heights, desired_segment);
    }
};

/// @brief Создание профиля с заданным постоянным шагом по координате на основе любого исходного профиля
inline pipe_profile_t create_uniform_profile(const pipe_profile_t& source_profile, double desired_uniform_segment)
{
    return pipe_profile_uniform::create_uniform_profile(source_profile, desired_uniform_segment);
}


/// @brief Вывод профиля в файл
/// @param coord_heights вектор, первый элемент которого сетка координат, а второй - соответствующие высотные отметки
/// @param filename путь и имя файла
inline void write_profile(const pipe_profile_t& profile, const string filename)
{
    std::ofstream output_file;
    size_t profCount = profile.get_point_count();
    output_file.open(filename + ".csv");
    output_file << "plug,km,height" << std::endl;
    for (size_t i = 0; i < profCount; i++)
    {
        output_file << "0," << profile.coordinates[i] << ',' << profile.heights[i] << std::endl;
    }
}



}
