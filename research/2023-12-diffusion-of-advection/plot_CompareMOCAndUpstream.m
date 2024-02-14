function main()
    clc
    clear
    currentDirectory = pwd;
    [upperPath, ~, ~] = fileparts(currentDirectory);
    [upperPath, ~, ~] = fileparts(upperPath);

    relativePath = fullfile('research_out', 'DiffusionOfAdvection', 'CompareUpstreamDiffAndMOC');

    [graphname, fullpath] = readFiles(upperPath, relativePath);

    plotGraphs(fullpath);

    plotPhysicalDiffusion(upperPath, relativePath);

    customizePlot(upperPath, relativePath, graphname); % Передаем upperPath и relativePath
end

function [graphname, fullpath] = readFiles(upperPath, relativePath)
    numfiles = length(dir(fullfile(upperPath, relativePath, '*.csv'))); 
    directory = dir(fullfile(upperPath, relativePath, '*.csv')); 
    fname = {directory.name};
    fullpath = cell(1, numfiles);
    graphname = cell(1, numfiles);
    for i = 1:numfiles
        fullpath{i} = fullfile(upperPath, relativePath, fname{i});
    end

    graphname = strrep(fname,'output','');
    graphname = strrep(graphname,' Cr=0.5.csv','');
    for i = 1:numfiles
        if contains(graphname{i}, ' Cr=1.csv')
            graphname{i} = 'Аналитическое решение';
        elseif contains(graphname{i}, ' physical')
            graphname{i} = 'Физическая диффузия';
        elseif contains(graphname{i}, ' MOC')
            graphname{i} = 'Метод характеристик';
        elseif contains(graphname{i}, ' upstream')
            graphname{i} = 'Метод передней разности';
        end
    end
    graphname{numfiles+1} = 'Границы физической диффузии';
end

function plotGraphs(fullpath)
    figure;
    numfiles = length(fullfile(fullpath)); 
    for i = 1:numfiles
        plotData(fullpath{i});
    end

    grid on;
    xlim([80, 81.5]);
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
    elseif contains(data, 'output MOC Cr=0.5')
        plot(t, y, 'LineWidth', 5, 'LineStyle','--');
    elseif contains(data, 'output upstream Cr=0.5')
        plot(t, y, 'LineWidth', 3, 'LineStyle','-.');
    else
        plot(t, y, 'LineWidth', 3);
    end
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

    fullpath = fullfile(upperPath, relativePath, 'Part 1.png');
    saveas(gcf, fullpath);
end