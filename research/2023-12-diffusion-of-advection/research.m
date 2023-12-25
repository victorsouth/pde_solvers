function main()
    clc
    clear
    currentDirectory = pwd;
    [upperPath, ~, ~] = fileparts(currentDirectory);
    [upperPath, ~, ~] = fileparts(upperPath);

    n = 0.05;
    relativePath = fullfile('research_out', 'DiffusionOfAdvection', 'CompareQuickestDiffusion');

    [fname, graphname, fullpath] = generateFileNames(upperPath, relativePath, n);

    plotGraphs(fullpath);

    plotPhysicalDiffusion(upperPath, relativePath);

    customizePlot(upperPath, relativePath, graphname); % Передаем upperPath и relativePath
end

function [fname, graphname, fullpath] = generateFileNames(upperPath, relativePath, n)
    fname = cell(1, 19);
    graphname = cell(1, 21);
    fullpath = cell(1, 19);

    for i = 1:19
        n_str = num2str(n);
        fname{i} = strcat('output Cr=', n_str, '.csv');
        graphname{i} = strcat('Cr = ', n_str);
        n = n + 0.05;
    end

    for i = 1:19
        fullpath{i} = fullfile(upperPath, relativePath, fname{i});
    end

    fname{20} = 'output physical.csv';
    fullpath{20} = fullfile(upperPath, relativePath, fname{20});
    graphname{20} = 'Физ. диффузия';

    fname{21} = 'physical_diffusion.txt';
    fullpath{21} = fullfile(upperPath, relativePath, fname{21});
    graphname{21} = 'Границы физ. диффузии';
end

function plotGraphs(fullpath)
    figure;

    for i = 1:19
        data = readtable(fullpath{i});
        plotData(data);
    end

    fname = 'output physical.csv';
    data = readtable(fullpath{20});
    plotData(data);


    grid on;
    xlim([80.5 81]);
    ylim([848 864]);
    xlabel('t, ч', 'FontSize', 18);
    ylabel('\rho, кг/м^3', 'FontSize', 18);
    hold off;
end

function plotData(data)
    t = table2array(data(:, 1));
    t = t / 3600;
    y = table2array(data(:, end));
    plot(t, y, 'LineWidth', 3, 'LineStyle', '--');
    hold on;
end

function plotPhysicalDiffusion(upperPath, relativePath)
    filename_pd = 'physical_diffusion.txt';
    fullpath_3 = fullfile(upperPath, relativePath, filename_pd);
    t = load(fullpath_3);
    t = t / 3600;
    xline(t(1), 'LineStyle', ':', 'LineWidth', 3, 'Color', 'r');
    xline(t(2), 'LineStyle', ':', 'LineWidth', 3, 'Color', 'r'); 

end

function customizePlot(upperPath, relativePath, graphname)
    fig = gcf;
    fig.Position = [0, 0, 1920, 1080];
    fig.Color = [1, 1, 1];
    legend(graphname);
    legend('show');
    lgd = legend;
    lgd.FontSize = 14;

    fullpath = fullfile(upperPath, relativePath, 'Part 3.png');
    saveas(gcf, fullpath);
end