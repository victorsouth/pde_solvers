clc
clear all
currentDirectory = pwd;

[upperPath, ~, ~] = fileparts(currentDirectory);
[upperPath, ~, ~] = fileparts(upperPath);

n=0.05;

relativePath = fullfile('research_out', 'DiffusionOfAdvection','CompareQuickestDiffusion');

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

fname = 'output physical.csv';
fullpath = fullfile(upperPath, relativePath, fname);
data = readtable(fullpath);
% Извлечение значений моментов времени
t = table2array(data(:, 1));
t = t/3600;
y = table2array(data(:, end));
% Построение графика
plot(t, y, 'LineWidth', 3, 'LineStyle', '--');
hold on;

graphname{20} = 'Физ. диффузия';



filename_pd = 'physical_diffusion.txt';
fullpath_3 = fullfile(upperPath, relativePath, filename_pd);
t = load(fullpath_3);
t = t/3600;
xline(t(1), LineStyle=":", LineWidth=3, Color="r")

xline(t(2), LineStyle=":", LineWidth=3, Color="r")
graphname{21} = 'Границы физ. диффузии';


legend(graphname)

grid on

xlim([80.5 81])

ylim([848 864])

xlabel('t, ч', 'FontSize',18)

ylabel('\rho, кг/м^3', 'FontSize',18)

hold off

% Настройка размера графика
fig = gcf;
fig.Position = [0, 0, 1920, 1080]; % [x, y, width, height]

% Изменение цвета фона
fig.Color = [1, 1, 1]; % Белый фон

% Добавление легенды и увеличение ее масштаба
legend('show');
lgd = legend;
lgd.FontSize = 14; % Размер шрифта
fullpath = fullfile(upperPath, relativePath, 'Part 3.png');

saveas(gcf, fullpath)