clc
clear all
currentDirectory = pwd;
% Используйте fileparts для разделения пути
[upperPath, ~, ~] = fileparts(currentDirectory);

fname = 'output Cr=0.5.csv';
relativePath_1 = fullfile('testing_out', 'QUICK.UseCaseStepDensity');
relativePath_2 = fullfile('testing_out', 'QUICKEST.UseCaseStepDensity');
relativePath_3 = fullfile('testing_out', 'QUICKEST_ULTIMATE.UseCaseStepDensity');
relativePath_4 = fullfile('testing_out', 'UpstreamDifferencing.UseCaseStepDensity');
relativePath_5 = fullfile('testing_out', 'MOC_Solver.MOC_Compare_With_QUICK');
fullpath_1 = fullfile(upperPath, relativePath_1, fname);
fullpath_2 = fullfile(upperPath, relativePath_2, fname);
fullpath_3 = fullfile(upperPath, relativePath_3, fname);
fullpath_4 = fullfile(upperPath, relativePath_4, fname);
fullpath_5 = fullfile(upperPath, relativePath_5, fname);
% Загрузка данных из CSV файлов
data_1 = readtable(fullpath_1);
data_2 = readtable(fullpath_2);
data_3 = readtable(fullpath_3);
data_4 = readtable(fullpath_4);
data_5 = readtable(fullpath_5);
% Значения для оси x
x = linspace(100, 700000, 7000);
% Извлечение значений моментов времени
t = table2array(data_1(:, 1));
% Находим номер строки в момент врмени начиная от 39350
index = find(t >= 39350 & t <= 39370);
% Извлечение значений для оси y (последний столбец - 20)
y = table2array(data_1(index, 6:end));

% Построение графика

plot(x, y, LineWidth = 1, LineStyle = "--", Color='b');
plot(x)
hold on;

% Значения для оси x
x = linspace(100, 700000, 7000);
% Извлечение значений моментов времени
t = table2array(data_2(:, 1));
% Находим номер строки в момент врмени начиная от 39350
index = find(t >= 39350 & t <= 39370);
% Извлечение значений для оси y (последний столбец - 20)
y = table2array(data_2(index, 6:end));

% Построение графика
plot(x, y, 'LineWidth', 1, "LineStyle","-.", Color='c');

hold on;

% Значения для оси x
x = linspace(100, 700000, 7000);
% Извлечение значений моментов времени
t = table2array(data_3(:, 1));
% Находим номер строки в момент врмени начиная от 39350
index = find(t >= 39350 & t <= 39370);
% Извлечение значений для оси y (последний столбец - 20)
y = table2array(data_3(index, 6:end));

% Построение графика
plot(x, y, 'LineWidth', 2, 'Color', "r");

hold on;

% Значения для оси x
x = linspace(100, 700000, 7000);
% Извлечение значений моментов времени
t = table2array(data_4(:, 1));
% Находим номер строки в момент врмени начиная от 39350
index = find(t >= 39350 & t <= 39370);
% Извлечение значений для оси y (последний столбец - 20)
y = table2array(data_4(index, 6:end));

% Построение графика
plot(x, y, LineStyle=":", LineWidth=2, Color='g');

hold on

% Значения для оси x
x = linspace(100, 700000, 7000);
% Извлечение значений моментов времени
t = table2array(data_5(:, 1));
% Находим номер строки в момент врмени начиная от 39350
index = find(t >= 39350 & t <= 39370);
% Извлечение значений для оси y (последний столбец - 20)
y = table2array(data_5(index, 6:end-1));

% Построение графика
plot(x, y, 'LineWidth', 2, Color='o');

hold on;

%xline(40000-2250, LineStyle="--", LineWidth=3, Color="r")

%xline(40000+2250, LineStyle="--", LineWidth=3, Color="r")

legend('QUICK','QUICKEST','QUICKEST ULTIMATE','UpstreamDifferencing', 'Метод характеристик')

grid on

xlim([88000 102000])

ylim([848 864])
