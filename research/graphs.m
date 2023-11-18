clc
clear all
currentDirectory = pwd;
% Используйте fileparts для разделения пути
[upperPath, ~, ~] = fileparts(currentDirectory);

fname = 'output Cr=0.9.csv';
relativePath_1 = fullfile('testing_out', 'QUICK.UseCaseStepDensity');
relativePath_2 = fullfile('testing_out', 'QUICKEST.UseCaseStepDensity');
relativePath_3 = fullfile('testing_out', 'QUICKEST_ULTIMATE.UseCaseStepDensity');
relativePath_4 = fullfile('testing_out', 'UpstreamDifferencing.UseCaseStepDensity');
fullpath_1 = fullfile(upperPath, relativePath_1, fname);
fullpath_2 = fullfile(upperPath, relativePath_2, fname);
fullpath_3 = fullfile(upperPath, relativePath_3, fname);
fullpath_4 = fullfile(upperPath, relativePath_4, fname);
% Загрузка данных из CSV файлов
data_1 = readtable(fullpath_1);
data_2 = readtable(fullpath_2);
data_3 = readtable(fullpath_3);
data_4 = readtable(fullpath_4);
% Значения для оси x
x = linspace(1000, 50000, 50);
% Извлечение значений моментов времени
t = table2array(data_1(:, 1));
% Находим номер строки в момент врмени начиная от 30000
index = find(t >= 30000 & t <= 31000);
% Извлечение значений для оси y (последний столбец - 20)
y = table2array(data_1(index, 6:end));

% Построение графика
plot(x, y, 'LineWidth', 1, 'Color', "#EDB120", "LineStyle","--");

hold on;

% Значения для оси x
x = linspace(1000, 50000, 50);
% Извлечение значений моментов времени
t = table2array(data_2(:, 1));
% Находим номер строки в момент врмени начиная от 30000
index = find(t >= 30000 & t <= 31000);
% Извлечение значений для оси y (последний столбец - 20)
y = table2array(data_2(index, 6:end));

% Построение графика
plot(x, y, 'LineWidth', 1, 'Color', "k", "LineStyle","-.");

hold on;

% Значения для оси x
x = linspace(1000, 50000, 50);
% Извлечение значений моментов времени
t = table2array(data_3(:, 1));
% Находим номер строки в момент врмени начиная от 30000
index = find(t >= 30000 & t <= 31000);
% Извлечение значений для оси y (последний столбец - 20)
y = table2array(data_3(index, 6:end));

% Построение графика
plot(x, y, 'LineWidth', 2, 'Color', "b");

hold on;

% Значения для оси x
x = linspace(1000, 50000, 50);
% Извлечение значений моментов времени
t = table2array(data_4(:, 1));
% Находим номер строки в момент врмени начиная от 30000
index = find(t >= 30000 & t <= 31000);
% Извлечение значений для оси y (последний столбец - 20)
y = table2array(data_4(index, 6:end));

% Построение графика
plot(x, y, LineStyle=":", LineWidth=2, Color="#77AC30");



hold on

xline(40000-2250, LineStyle="--", LineWidth=3, Color="r")

xline(40000+2250, LineStyle="--", LineWidth=3, Color="r")

legend('QUICK','QUICKEST','QUICKEST ULTIMATE','UpstreamDifferencing', 'Границы физ. диффузии')

grid on

xlim([35000 50000])

ylim([848 864])
