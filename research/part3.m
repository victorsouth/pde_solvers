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
    t = t/3600;
    y = table2array(data(:, end));
% Построение графика
    plot(t, y, 'LineWidth', 3, 'LineStyle', '--');
    hold on;
end
data = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
v = 2.4096;
diam_vnutr = 0.514;
nu = 15 * 10^-6;
Re = v*diam_vnutr/nu;
sherokh = 15*10^-5;
lambda = 0.11 * (sherokh + 68/Re)^0.25;
K_CnaV = 3.211 * sqrt(lambda) * diam_vnutr;
L = 6.58 * sqrt(K_CnaV) * sqrt(700*10^3);
t_diff=L/v;
fname_1 = 'time2.txt';
fname_2 = 'CL2.txt';
relativePath_3 = fullfile('msvc');
fullpath_3_1 = fullfile(upperPath, relativePath_3, fname_1);
fullpath_3_2 = fullfile(upperPath, relativePath_3, fname_2);
t = load(fullpath_3_1);
t = t/3600;
y = load(fullpath_3_2);
index = find(y >= (max(y)+min(y))/2,1);
plot(t, y, LineWidth = 2, LineStyle = "-.", Color='r');
hold on;
graphname{20} = 'Физ. диффузия';

xline(t(index)+t_diff/7200, LineStyle=":", LineWidth=3, Color="r")

xline(t(index)-t_diff/7200, LineStyle=":", LineWidth=3, Color="r")
graphname{21} = 'Границы физ. диффузии';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


legend(graphname)

grid on

xlim([80.5 81])

ylim([848 864])

xlabel('t, ч', 'FontSize',18)

ylabel('\rho, кг/м^3', 'FontSize',18)