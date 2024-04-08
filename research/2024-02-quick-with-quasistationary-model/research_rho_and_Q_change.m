function main()
    clc
    clear
    import matlab.io.*;
    currentDirectory = pwd;
    [upperPath, ~, ~] = fileparts(currentDirectory);
    [upperPath, ~, ~] = fileparts(upperPath);   
    relativePath_1 = fullfile('research_out', 'QSM_models', 'MocWithQuasiStationaryModel', 'WorkingWithTimeSeries');
    
    % Загрузка и обработка данных из первой папки
    processDirectoryData(upperPath, relativePath_1);

    relativePath_2 = fullfile('research_out', 'QSM_models', 'QuickWithQuasiStationaryModel', 'WorkingWithTimeSeries');

    % Загрузка и обработка данных из второй папки
    processDirectoryData(upperPath, relativePath_2);
end

function processDirectoryData(upperPath, relativePath)
    % Загрузка данных из файлов CSV
    [data, data2, data3] = loadDataFromFiles(upperPath, relativePath);
    
    % Нахождение минимальных и максимальных значений
    [minValue, maxValue] = calculateMinMax(data);
    [minValue2, maxValue2] = calculateMinMax2(data2);
    [minValue3, maxValue3] = calculateMinMax3(data3);
    [minValue4, maxValue4] = calculateMinMax4(data2);

    % Отображение данных перед началом цикла и получение гифки
    plotData(data, data2, data3, km, minValue, maxValue, minValue2, maxValue2, minValue3, maxValue3, minValue4, maxValue4);
    createGif(data, data2, data3, km, minValue, maxValue, minValue2, maxValue2, minValue3, maxValue3, minValue4, maxValue4);
end

function [data, km] = loadData(filename)
    % Загрузка данных из файла CSV
    opts = detectImportOptions(filename);
    opts.DataLines = [1, Inf]; % Загрузить все строки
    opts.Delimiter = ';'; % Указываем разделитель поля
    opts.VariableNamesLine = 1; % Первая строка содержит имена переменных
    
    % Указываем типы переменных
    numColumns = 21; % Укажите фактическое количество столбцов в вашем файле CSV
    variableTypes = cell(1, numColumns);
    variableTypes{1} = 'datetime'; % Первая переменная - дата и время
    for i = 2:numColumns
        variableTypes{i} = 'double'; % Все остальные переменные - числовые
    end
    opts.VariableTypes = variableTypes;

    % Чтение таблицы
    data = readtable(filename, opts);

    % Преобразование времени в числовой формат
    data.t = datenum(data.t);
    data.t = (data.t - data.t(1)) * 24 * 3600; % Преобразование в секунды от начала

    % Получение значений km из данных (пример)
    km = 0:0.1:200;
end

function [data, data2, data3] = loadDataFromFiles(upperPath, relativePath)
    % Загрузка данных из файлов CSV
    data = loadData(fullfile(upperPath, relativePath, 'output pressure_delta.csv'));
    data2 = loadData(fullfile(upperPath, relativePath, 'output pressure.csv'));
    data3 = loadData(fullfile(upperPath, relativePath, 'density.csv.csv'));
end
function [minValue, maxValue] = calculateMinMax(data)
    % Нахождение минимального и максимального значения для data
    minValue = min(data(:, 2:end-1), [], 'all') - max(data(:, 2:end-1), [], 'all') * 0.1;
    maxValue = max(data(:, 2:end-1), [], 'all') + max(data(:, 2:end-1), [], 'all') * 2;
end

function [minValue2, maxValue2] = calculateMinMax2(data)
    % Нахождение минимального и максимального значения для data2
    minValue2 = min(data(:, 2:end-1), [], 'all') - 0.1e6;
    maxValue2 = max(data(:, 2:end-1), [], 'all') + 0.1e6;
end

function [minValue3, maxValue3] = calculateMinMax3(data)
    % Нахождение минимального и максимального значения для data3
    minValue3 = min(data(:, 2:end-1), [], 'all') - 10;
    maxValue3 = max(data(:, 2:end-1), [], 'all') + 10;
end

function [minValue4, maxValue4] = calculateMinMax4(data)
    % Нахождение минимального и максимального значения для data3
    minValue4 = min(data(:, end-1), [], 'all') - max(data(:, end-1), [], 'all')*0.1;
    maxValue4 = max(data(:, end-1), [], 'all') + max(data(:, end-1), [], 'all')*0.1;
end

function plotData(data, data2, data3, km, minValue, maxValue, minValue2, maxValue2, minValue3, maxValue3, minValue4, maxValue4)
    figure;
    % Первый подграфик
    subplot(5, 1, 2);
    plot(km, data(1, 2:end-1), 'Color', 'b', LineWidth=2);
    hold on;
    xlabel('Труба, км');
    ylabel('Разница давлений, Па');
    title('Профиль разницы давлений');
    xlim([0, 100]);
    ylim([minValue, maxValue]);

    % Второй подграфик
    subplot(5, 1, 1);
    plot(km, data2(1, 2:end-1), 'Color', 'b');
    hold on;
    xlabel('Труба, км');
    ylabel('Давление, Па');
    title('Профиль давления');
    xlim([0, 100]);
    ylim([minValue2, maxValue2]);

    % Третий подграфик
    subplot(5, 1, 3);
    plot(km, data3(1, 2:end-1), 'Color', 'b', LineWidth=2);
    xlabel('Труба, км');
    ylabel('Плотность, кг/м3');
    title('Профиль плотности');
    xlim([0, 100]);
    ylim([minValue3, maxValue3]);
    
    % Четвертый подграфик
    subplot(5, 1, 4);
    % Подготовка данных
    t = data2(1:end, 1);
    t = t/3600;
    plot(t, data2(1:end, end-1), 'Color', 'b', LineWidth=2);
    hold on;
    plot(t(1), data2(1, end-1), "Marker",".","LineStyle","none",MarkerSize=20, Color='r');
    hold on;
    text(t(1), data2(1, end-1), ['t = ' num2str(data2(1,1)) ', с'], 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
    hold off;
    xlabel('Время, ч');
    ylabel('Давление, Па');
    title(['Времяной ряд давления на выходе']);
    xlim([0, 42]);
    ylim([minValue4,  maxValue4]);
end

function createGif(data, data2, data3, km, minValue, maxValue, minValue2, maxValue2, minValue3, maxValue3, minValue4, maxValue4)
    % Получение кадра для первого кадра
    frame = getframe(gcf);
    [im, map] = rgb2ind(frame.cdata, 256, 'nodither');
    im(1, 1, 1, size(data, 1)) = 0;

    % Анимация
    for i = 2:size(data(:,1))
        % Отображение данных
        subplot(5, 1, 2);
        plot(km, data(i, 2:end-1), 'Color', 'r', LineWidth=2);
        hold on;
        plot(km, data(1, 2:end-1), 'Color', 'b', LineWidth=2);
        xlim([0, 100]);
        ylim([minValue, maxValue]);
        xlabel('Труба, км');
        ylabel('Разница давлений, Па');
        title('Профиль разницы давлений');
        hold off;

        subplot(5, 1, 1);
        plot(km, data2(i, 2:end-1), 'Color', 'r');
        hold on;
        plot(km, data2(1, 2:end-1), 'Color', 'b');
        xlim([0, 100]);
        ylim([minValue2, maxValue2]);
        xlabel('Труба, км');
        ylabel('Давление, Па');
        title('Профиль давления');
        hold off;

        subplot(5, 1, 3);
        plot(km, data3(i, 2:end-1), 'Color', 'b', LineWidth=2);
        xlabel('Труба, км');
        ylabel('Плотность, кг/м3');
        title('Профиль плотности');
        xlim([0, 100]);
        ylim([minValue3, maxValue3]);
        
        subplot(5, 1, 4);
        % Подготовка данных
        t = data2(1:end, 1);
        t = t/3600;
        plot(t, data2(1:end, end-1), 'Color', 'b', LineWidth=2);
        hold on
        plot(t(i), data2(i, end-1), "Marker",".","LineStyle","none",MarkerSize=20, Color='r')
        hold on;
        text(t(i), data2(i, end-1), ['t = ' num2str(data2(i,1)) ', с'], 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
        hold off;
        xlabel('Время, ч');
        ylabel('Давление, Па');
        title(['Времяной ряд давления на выходе']);
        xlim([0, 42]);
        ylim([minValue4,maxValue4]);
    end

    % Сохранение гифки в файл
    filename = 'меняем плотность и расход_3.gif';
    imwrite(im, map, filename, 'DelayTime', 0.02, 'LoopCount', inf);
    disp(['Гифка сохранена в файл: ' filename]);
end