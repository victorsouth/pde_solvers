clc
clear all
currentDirectory = pwd;
% Используйте fileparts для разделения пути
[upperPath, ~, ~] = fileparts(currentDirectory);
fname = 'output Cr=0.5.csv';
fname_1 = 'output cr=1.csv';
relativePath_1 = fullfile('testing_out', 'UpstreamDifferencing.UseCaseStepDensity');
relativePath_2 = fullfile('testing_out', 'MOC_Solver.MOC_Compare_With_QUICK');
relativePath_3 = fullfile('testing_out', 'QUICKEST_ULTIMATE.UseCaseStepDensity');
fullpath_1 = fullfile(upperPath, relativePath_1, fname);
fullpath_2 = fullfile(upperPath, relativePath_2, fname);
fullpath_3 = fullfile(upperPath, relativePath_3, fname_1);
% Загрузка данных из CSV файлов
data_1 = readtable(fullpath_1);
data_2 = readtable(fullpath_2);

% Извлечение значений моментов времени
t = table2array(data_1(:, 1));
% Находим номер строки в момент врмени начиная от 39350
y = table2array(data_1(:, end));

% Построение графика

plot(t, y, 'x', 'MarkerSize', 3, 'Color', [0, 0, 1.0]);
hold on;

% Извлечение значений моментов времени
t = table2array(data_2(:, 1));
% Находим номер строки в момент врмени начиная от 39350
y = table2array(data_2(:, end));

% Построение графика

plot(t, y, 'LineWidth', 2, "LineStyle","-.", Color='g');
hold on;

fname_1 = 'time2.txt';
fname_2 = 'CL2.txt';
relativePath_3 = fullfile('msvc');
fullpath_3_1 = fullfile(upperPath, relativePath_3, fname_1);
fullpath_3_2 = fullfile(upperPath, relativePath_3, fname_2);
t = load(fullpath_3_1);
y = load(fullpath_3_2);



plot(t, y, LineWidth = 2, LineStyle = "--", Color='r');
hold on;

data_3 = readtable(fullpath_3);

% Извлечение значений моментов времени
t = table2array(data_3(:, 1));
% Находим номер строки в момент врмени начиная от 39350
y = table2array(data_3(:, end));

% Построение графика

plot(t, y, LineWidth = 1, LineStyle = "-", Color='0 0 0');
hold on;

%xline(t(27), LineStyle=":", LineWidth=3, Color="r")

%xline(t(58), LineStyle=":", LineWidth=3, Color="r")

legend('UpstreamDifferencing', 'Метод характеристик', 'Физическая диффузия', 'Аналитическое решение')

grid on

xlim([284000 298000])

ylim([848 864])

xlabel('t')

ylabel('rho')