function main()
    clc
    clear
    currentDirectory = pwd;
    [upperPath, ~, ~] = fileparts(currentDirectory);
    [upperPath, ~, ~] = fileparts(upperPath);

    relativePath = fullfile('research_out', 'DiffusionOfAdvection', 'CompareQuickestAndQuickestUltimateDiffusion');

    [graphname, fullpath] = readFiles(upperPath, relativePath);

    plotGraphs(fullpath);

    customizePlot(upperPath, relativePath, graphname); % Передаем upperPath и relativePath
end

function [graphname, fullpath] = readFiles(upperPath, relativePath)
    numfiles = length(dir(fullfile(upperPath, relativePath, '*.csv'))); 
    directory = dir(fullfile(upperPath, relativePath, '*.csv')); 
    fname = {directory.name};
    graphname = cell(1, numfiles);
    fullpath = cell(1, numfiles);
   
    for i = 1:numfiles
        fullpath{i} = fullfile(upperPath, relativePath, fname{i});
    end

    graphname = strrep(fname,'output','');
    graphname = strrep(graphname,' Cr=0.5.csv','');
    for i = 1:numfiles
        if contains(graphname{i}, ' Cr=1.csv')
            graphname{i} = 'Аналитическое решение';
        end
    end
end

function plotGraphs(fullpath)
    figure;
    numfiles = length(fullfile(fullpath)); 
    for i = 1:numfiles
        plotData(fullpath{i});
    end

    grid on;
    xlim([79.5 82]);
    ylim([848 864]);
    xlabel('t, ч', 'FontSize', 18);
    ylabel('\rho, кг/м^3', 'FontSize', 18);
    hold off;
end

function plotData(data)
    table = readtable(data);
    t = table2array(table(:, 1));
    t = t / 3600;
    y = table2array(table(:, end));
    if contains(data, ' Cr=1.csv')
        plot(t, y, LineWidth = 1, LineStyle = "-", Color='0 0 0');
    else
        plot(t, y, 'LineWidth', 3);
    end
    hold on;
end


function customizePlot(upperPath, relativePath, graphname)
    fig = gcf;
    fig.Position = [0, 0, 1920, 1080];
    fig.Color = [1, 1, 1];
    legend(graphname);
    legend('show');
    lgd = legend;
    lgd.FontSize = 14;

    fullpath = fullfile(upperPath, relativePath, 'Part 2.png');
    saveas(gcf, fullpath);
end