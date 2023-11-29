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

% Извлечение значений моментов времени
t = table2array(data_1(:, 1));
% Находим номер строки в момент врмени начиная от 39350
y = table2array(data_1(:, end));

% Построение графика

plot(t, y, LineWidth = 1, LineStyle = "--", Color='b');
hold on;

% Извлечение значений моментов времени
t = table2array(data_2(:, 1));
% Находим номер строки в момент врмени начиная от 39350
y = table2array(data_2(:, end));

% Построение графика
plot(t, y, 'LineWidth', 1, "LineStyle","-.", Color='b');

hold on;

% Извлечение значений моментов времени
t = table2array(data_3(:, 1));
% Находим номер строки в момент врмени начиная от 39350
y = table2array(data_3(:, end));

% Построение графика
plot(t, y, 'LineWidth', 1, 'Color', "m");

hold on;

% Извлечение значений моментов времени
t = table2array(data_4(:, 1));
% Находим номер строки в момент врмени начиная от 39350
y = table2array(data_4(:, end));

% Построение графика
plot(t, y, LineStyle=":", LineWidth=2, Color='g');

hold on

% Извлечение значений моментов времени
t = table2array(data_5(:, 1));
% Находим номер строки в момент врмени начиная от 39350
y = table2array(data_5(:, end));

% Построение графика
plot(t, y, 'LineWidth', 2, Color='yellow');

hold on;


fname_1 = 'time2.txt';
fname_2 = 'CL2.txt';
relativePath_6 = fullfile('msvc');
fullpath_6_1 = fullfile(upperPath, relativePath_6, fname_1);
fullpath_6_2 = fullfile(upperPath, relativePath_6, fname_2);
t = load(fullpath_6_1);
y = load(fullpath_6_2);

plot(t, y, LineWidth = 2, LineStyle = "--", Color='r');
hold on;
%xline(40000-2250, LineStyle="--", LineWidth=3, Color="r")

%xline(40000+2250, LineStyle="--", LineWidth=3, Color="r")

legend('QUICK','QUICKEST','QUICKEST ULTIMATE','UpstreamDifferencing', 'Метод характеристик', 'Физическая диффузия')

grid on

xlim([280000 300000])

ylim([848 864])
