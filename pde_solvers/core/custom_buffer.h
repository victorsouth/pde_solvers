#pragma once

/// @brief Специфический слой
/// @tparam LayerType 
template <typename LayerType>
class custom_buffer_t
{
    vector<LayerType> layers;
    size_t current_layer{ 0 };
protected:
    size_t advanced_layer_index(int offset) const {
        return (current_layer + layers.size() + offset) % layers.size();
    }
    size_t previous_layer_index() const {
        return (current_layer + layers.size() - 1) % layers.size();
    }
public:
    custom_buffer_t(size_t layer_count, const LayerType& layer)
        : layers(layer_count, layer)
    {
    }
    custom_buffer_t(size_t layer_count, size_t profile_length)
        : layers(layer_count, profile_length)
    {
    }

    const vector<LayerType>& get_layers() const {
        return layers;
    }

    vector<vector<double>> get_layers_cell_values() const {
        vector<vector<double>> result;
        transform(layers.begin(), layers.end(),
            back_inserter(result),
            [](const LayerType& layer) { return layer.cell.value; });

        return result;
    }


    void advance(int offset) {
        current_layer = advanced_layer_index(offset);
    }
    LayerType& current() { return layers[current_layer]; }
    const LayerType& current() const { return layers[current_layer]; }

    LayerType& previous() { return layers[previous_layer_index()]; }
    const LayerType& previous() const { return layers[previous_layer_index()]; }
};
