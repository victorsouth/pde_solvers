%% Достает профиль из данных (data), получает оттуда количество точек N, 
% вычисляет шаг расчетной сетки и формирует саму сетку (вектор x)
function [x, N] = get_grid(Length, data, gridof)
names = keys(data);



profiles = data(names{1}); % массив профилей для первой переменной группы
profile = profiles.get(0);

gridtype = gridof(names{1});
if gridtype == 'cells'
    % N - количество точек, их на одну больше для профиля ячеек, ...
    N = length(profile) + 1; 
else
    % ... а если это профиль точек, то N - сразу то, что надо
    N = length(profile);
end


dx = Length/(N - 1); % Шаг по координате
x = 0:dx:Length;


end

