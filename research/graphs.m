clc
clear all

fname = 'output Cr=0.9.csv';
relativePath_1 = fullfile('testing_out', 'QUICK.UseCaseStepDensity');
relativePath_2 = fullfile('testing_out', 'QUICKEST.UseCaseStepDensity');
relativePath_3 = fullfile('testing_out', 'QUICKEST_ULTIMATE.UseCaseStepDensity');
relativePath_4 = fullfile('testing_out', 'UpstreamDifferencing.UseCaseStepDensity');
fullpath_1 = fullfile(pwd, relativePath_1);
fullpath_2 = fullfile(pwd, relativePath_2);
fullpath_3 = fullfile(pwd, relativePath_3);
fullpath_4 = fullfile(pwd, relativePath_4);
% Загрузка данных из CSV файлов
data_1 = readtable(fullpath_1);
data_2 = readtable(fullpath_2);
data_3 = readtable(fullpath_3);
data_4 = readtable(fullpath_4);
% Извлечение значений для оси x (первый столбец)
x = table2array(data_1(:, 1));

% Извлечение значений для оси y (последний столбец - 20)
y = table2array(data_1(:, end-20));

% Построение графика
plot(x, y, 'LineWidth', 1, 'Color', "#EDB120", "LineStyle","--");

hold on;

% Извлечение значений для оси x (первый столбец)
x = table2array(data_2(:, 1));

% Извлечение значений для оси y (последний столбец - 20)
y = table2array(data_2(:, end-20));

% Построение графика
plot(x, y, 'LineWidth', 1, 'Color', "k", "LineStyle","-.");

hold on;

% Извлечение значений для оси x (первый столбец)
x = table2array(data_3(:, 1));

% Извлечение значений для оси y (последний столбец - 20)
y = table2array(data_3(:, end-20));

% Построение графика
plot(x, y, 'LineWidth', 2, 'Color', "b");

hold on;

% Извлечение значений для оси x (первый столбец)
x = table2array(data_4(:, 1));

% Извлечение значений для оси y (последний столбец - 20)
y = table2array(data_4(:, end-20));

% Построение графика
plot(x, y, LineStyle=":", LineWidth=2, Color="#77AC30");



hold on

xline(22783-2250, LineStyle="--", LineWidth=3, Color="r")

xline(22783+2250, LineStyle="--", LineWidth=3, Color="r")

legend('QUICK','QUICKEST','QUICKEST ULTIMATE','UpstreamDifferencing', 'Границы физ. диффузии')

grid on

ylim([848 864])
