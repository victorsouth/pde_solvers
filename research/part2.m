clc
clear all
currentDirectory = pwd;
% Используйте fileparts для разделения пути
[upperPath, ~, ~] = fileparts(currentDirectory);

fname = 'output Cr=0.5.csv';
fname_1 = 'output Cr=1.csv';
relativePath_2 = fullfile('testing_out', 'QUICKEST.UseCaseStepDensity');
relativePath_3 = fullfile('testing_out', 'QUICKEST_ULTIMATE.UseCaseStepDensity');
fullpath_1 = fullfile(upperPath, relativePath_3, fname_1);
fullpath_2 = fullfile(upperPath, relativePath_2, fname);
fullpath_3 = fullfile(upperPath, relativePath_3, fname);
% Загрузка данных из CSV файлов
data_2 = readtable(fullpath_2);
data_3 = readtable(fullpath_3);


% Извлечение значений моментов времени
t = table2array(data_2(:, 1));
% Находим номер строки в момент врмени начиная от 39350
y = table2array(data_2(:, end));

% Построение графика
plot(t, y, 'x', 'MarkerSize', 3, 'Color', [0, 0, 1.0]);

hold on;

% Извлечение значений моментов времени
t = table2array(data_3(:, 1));
% Находим номер строки в момент врмени начиная от 39350
y = table2array(data_3(:, end));

% Построение графика
plot(t, y, 'LineWidth', 2, "LineStyle","-.", Color='g');

hold on;


data_1 = readtable(fullpath_1);

% Извлечение значений моментов времени
t = table2array(data_1(:, 1));
% Находим номер строки в момент врмени начиная от 39350
y = table2array(data_1(:, end));

% Построение графика

plot(t, y, LineWidth = 1, LineStyle = "-", Color='0 0 0');
hold on;

legend('QUICKEST','QUICKEST ULTIMATE', 'Аналитическое решение')

grid on

xlim([284000 298000])

ylim([848 864])

xlabel('t')

ylabel('rho')
