#pragma once

/// @brief Контейнер слоев с удобным доступом при численном расчете задач на ДУЧП и им подобным
/// Организует циклическую смену буферов
/// Наиболее типичное использование - организация двух слоев 
/// как текущий/предыдущий с простым переключением
/// @tparam LayerType Тип слоя
template <typename LayerType>
class layer_container_t {
    /// @brief Буфер слоев
    vector<LayerType> layers;
    /// @brief Индекс текущего слоя
    size_t current_layer{ 0 };
protected:
    /// @brief Рассчитывает индекс слоя в layers на основе 
    /// смещения offset от текущего индекса current_layer
    size_t advanced_layer_index(int offset) const {
        return (current_layer + layers.size() + offset) % layers.size();
    }
    /// @brief Эквивалентен advanced_layer_index(-1)
    size_t previous_layer_index() const {
        return (current_layer + layers.size() - 1) % layers.size();
    }
public:
    /// @brief Конструктор с инициализацией буфера слоев по переданному слою layer
    /// @param layer_count Количество слоев в буфере
    /// @param layer Слой для инициализации
    layer_container_t(size_t layer_count, const LayerType& layer)
        : layers(layer_count, layer)
    {
    }
    /// @brief Конструктор с инициализацией буфера по размерности каждого слоя
    /// @param layer_count Количество
    /// @param profile_length Передается в конструктор LayerType 
    /// TODO: уточнить, перепроверить
    layer_container_t(size_t layer_count, size_t profile_length)
        : layers(layer_count, profile_length)
    {
    }
    /// @brief Возвращает внутренний буфер слоев 
    const vector<LayerType>& get_layers() const {
        return layers;
    }
    /// @brief Частный геттер, возвращающий из слоев ячейки
    /// TODO: зачем он нужен, такой частный?
    vector<vector<double>> get_layers_cell_values() const {
        vector<vector<double>> result;
        transform(layers.begin(), layers.end(),
            back_inserter(result),
            [](const LayerType& layer) { return layer.cell.value; });

        return result;
    }
    /// @brief Смещает внутренний
    /// @param offset 
    void advance(int offset) {
        current_layer = advanced_layer_index(offset);
    }
    /// @brief Ссылка на текущий слой 
    LayerType& current() { return layers[current_layer]; }
    /// @brief Константная ссылка на текущий слой
    const LayerType& current() const { return layers[current_layer]; }
    /// @brief Ссылка на предыдущий слой 
    LayerType& previous() { return layers[previous_layer_index()]; }
    /// @brief Константная ссылка на предыдущий слой 
    const LayerType& previous() const { return layers[previous_layer_index()]; }
};
