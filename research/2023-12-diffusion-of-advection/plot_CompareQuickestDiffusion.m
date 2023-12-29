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

    %количество файлов - 2 файла, output physical и output Cr=1
    numfiles = length(dir(fullfile(upperPath, relativePath, '*.csv'))) - 2; 
    fname = cell(1, numfiles);
    graphname = cell(1, numfiles + 2);
    fullpath = cell(1, numfiles);

    for i = 1:numfiles
        n_str = num2str(n);
        fname{i} = strcat('output Cr=', n_str, '.csv');
        graphname{i} = strcat('Cr = ', n_str);
        n = n + 0.05;
    end

    for i = 1:numfiles
        fullpath{i} = fullfile(upperPath, relativePath, fname{i});
    end

    fname{numfiles+1} = 'output physical.csv';
    fullpath{numfiles+1} = fullfile(upperPath, relativePath, fname{numfiles+1});
    graphname{numfiles+1} = 'Физ. диффузия';
    graphname{numfiles+2} = 'Границы физ. диффузии';
end

function plotGraphs(fullpath)
    figure;
    numfiles = length(fullfile(fullpath)); 
    for i = 1:numfiles
        data = readtable(fullpath{i});
        plotData(data);
    end

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