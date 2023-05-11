%% Парсит сырые данные при выводе
% Возвращает: data = Map<ArrayList<Matlab Vector>> 
function [time, data, unitof, group2vars, lims, gridof] = parse_csv_raw(raw_data)
% Материалы по Java Map, ArrayList
% https://vertex-academy.com/tutorials/ru/list-java-primer/

% читаем профили для разных моментов времени и группируем их по переменным

n = size(raw_data, 1);
data = containers.Map; 
time = zeros(n, 1);
lims = containers.Map;

for i = 1:n
    [t, name, ~, ~, values] = parse_row(raw_data, i);
    if ~data.isKey(name)
        data(name) = java.util.ArrayList;
        lims(name) = [+inf; -inf];
    end
   
    ylimits = lims(name);
    ymin = min(ylimits(1), min(values));
    ymax = max(ylimits(2), max(values));
    lims(name) = [ymin; ymax];
    
    data(name).add(values);
    T = data(name).size;
    time(T) = t;
    
end

% реализация через Java.Hashtable
% https://www.cs.utah.edu/~germain/PPS/Topics/Matlab/maps_via_java.html
%data = java.util.Hashtable; 
% for i = 1:n
%     [name, ~, ~, values] = parse_row(raw_data, i);
%     data_form_name = data.get(name); 
%     if (isempty(data_form_name)) 
%         data.put(name, java.util.ArrayList);
%     end
%     x = data.get(name);
%     x.add(values);
% end

% готовим размерности и группы
%names = unique(raw_data(:, 1)); % по идее, можно просто из data достать
layer_length = data.Count;

time = time(1:(n / layer_length));
group2vars = containers.Map('KeyType', 'int64', 'ValueType', 'any');  
unitof = containers.Map;  
gridof = containers.Map;  
for i = 1:layer_length
    [~, name, group, unit, ~, grid_type] = parse_row(raw_data, i);
    if ~group2vars.isKey(group)
        group2vars(group) = {name};
    else
        group_names = group2vars(group);
        group_names{length(group_names) + 1} = name;
        group2vars(group) = group_names;
    end
    
    gridof(name) = grid_type;
    unitof(name) = unit;
end


end

%%
function [t, name, group, unit, values, grid_type] = parse_row(Data, row)

t = Data{row, 1};
name = Data{row, 2};
if strcmp(name, 'cells') || strcmp(name, 'points')
    grid_type = name;
    name = Data{row, 3};
    group = int32(Data{row, 4});
    unit = Data{row, 5};
    values = Data(row, 6:end);
    values = cell2mat(values);
else
    grid_type = 'points';
    group = int32(Data{row, 3});
    unit = Data{row, 4};
    values = Data(row, 5:end);
    values = cell2mat(values);
end
    


end
