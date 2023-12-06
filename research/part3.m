clc
clear all
currentDirectory = pwd;
% Используйте fileparts для разделения пути
[upperPath, ~, ~] = fileparts(currentDirectory);

n=0.05;

relativePath = fullfile('testing_out', 'QUICKEST_ULTIMATE.UseCaseStepDensity');

for i = 1 : 19
    n_str = num2str(n);
    fname{i} = strcat('output Cr=',n_str,'.csv');
    graphname{i} = strcat('Cr = ',n_str);
    n = n + 0.05;
end

for i = 1 : 19
    fullpath{i} = fullfile(upperPath, relativePath, fname{i});
end

for i = 1 : 19
    data = readtable(fullpath{i});
% Извлечение значений моментов времени
    t = table2array(data(:, 1));
    y = table2array(data(:, end));
% Построение графика
    plot(t, y, 'LineWidth', 3, 'LineStyle', '--');
    hold on;
end
data_1 = 0;

fname_1 = 'time2.txt';
fname_2 = 'CL2.txt';
relativePath_3 = fullfile('msvc');
fullpath_3_1 = fullfile(upperPath, relativePath_3, fname_1);
fullpath_3_2 = fullfile(upperPath, relativePath_3, fname_2);
t = load(fullpath_3_1);
y = load(fullpath_3_2);

plot(t, y, LineWidth = 2, LineStyle = "-.", Color='r');
hold on;
graphname{20} = 'Физ. диффузия';

legend(graphname)

grid on

xlim([284000 298000])

ylim([848 864])

xlabel('t')

ylabel('rho')