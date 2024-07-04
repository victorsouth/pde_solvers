function main()
    clc
    clear
    import matlab.io.*;
    currentDirectory = pwd;
    [upperPath, ~, ~] = fileparts(currentDirectory);
    [upperPath, ~, ~] = fileparts(upperPath);
    relativePath = fullfile('research_out', 'QSM_models', 'QuasiStationaryModel');

    % Загрузка и обработка данных из первой папки
    processDirectoryData(upperPath, relativePath, 'MocWithQuasiStationaryModel');
    %processDirectoryData(upperPath, relativePath, 'ChangeFlowMocWithQuasiStationaryModel');
    % Загрузка и обработка данных из второй папки
    %processDirectoryData(upperPath, relativePath, 'QuickWithQuasiStationaryModel');
    %processDirectoryData(upperPath, relativePath, 'IdealQuickWithQuasiStationaryModel');
    %processDirectoryData(upperPath, relativePath, 'OptionalStepMocWithQuasiStationaryModel');
    %processDirectoryData(upperPath, relativePath, 'IdealMocWithQuasiStationaryModel');
    %processDirectoryData(upperPath, relativePath, 'ChangeFlowQuickWithQuasiStationaryModel');
    %processDirectoryData(upperPath, relativePath, 'ChangeFlowOptionalStepMocWithQuasiStationaryModel');
end

function processDirectoryData(upperPath, relativePath, name)
    % Загрузка данных из файлов CSV
    relativePath = fullfile(relativePath,name);
    [data, data2, data3] = loadDataFromFiles(upperPath, relativePath);
    % Длина трубы от 0 до 200 км
    km = linspace(0, 200, size(data, 2) - 1);
    % Нахождение минимальных и максимальных значений
    [minValue, maxValue] = calculateMinMax(data);
    [minValue2, maxValue2] = calculateMinMax2(data2);
    [minValue3, maxValue3] = calculateMinMax3(data3);
    [minValue4, maxValue4] = calculateMinMax4(data2);
    %minValue4 = 2.43e+06;
    %maxValue4 = 2.44e+06;
    % Отображение данных перед началом цикла и получение гифки
    %plotData(data, data2, data3, km, minValue, maxValue, minValue2, maxValue2, minValue3, maxValue3, minValue4, maxValue4, name);
    %createGif(data, data2, data3, km, minValue, maxValue, minValue2, maxValue2, minValue3, maxValue3, minValue4, maxValue4, name);
    plotScaledData(data, data3, km, minValue3, maxValue3, name);
end

function [data, data2, data3] = loadDataFromFiles(upperPath, relativePath)
    % Полные пути к файлам
    filePath1 = fullfile(upperPath, relativePath, 'output pressure_delta.csv');
    filePath2 = fullfile(upperPath, relativePath, 'output pressure.csv');
    filePath3 = fullfile(upperPath, relativePath, 'output density.csv');


    % Загрузка данных из файлов CSV
    data = readtable(filePath1);
    data2 = readtable(filePath2);
    data3 = readtable(filePath3);
end

function [minValue, maxValue] = calculateMinMax(data)
    % Нахождение минимального и максимального значения для data
    numericData = data{:, 2:end-1}; % Извлечение числовых данных из таблицы
    minValue = min(numericData, [], 'all') - max(numericData, [], 'all') * 0.1;
    maxValue = max(numericData, [], 'all') + max(numericData, [], 'all') * 0.1;
end

function [minValue2, maxValue2] = calculateMinMax2(data)
    % Нахождение минимального и максимального значения для data2
    numericData = data{:, 2:end-1}; % Извлечение числовых данных из таблицы
    minValue2 = min(numericData, [], 'all') - 0.1e6;
    maxValue2 = max(numericData, [], 'all') + 0.1e6;
end

function [minValue3, maxValue3] = calculateMinMax3(data)
    % Нахождение минимального и максимального значения для data3
    numericData = data{:, 2:end-1}; % Извлечение числовых данных из таблицы
    minValue3 = min(numericData, [], 'all') - 10;
    maxValue3 = max(numericData, [], 'all') + 10;
end

function [minValue4, maxValue4] = calculateMinMax4(data)
    % Извлечение последнего столбца числовых данных из таблицы
    numericData = data{:, end};

    % Нахождение минимального и максимального значения только в последнем столбце
    minValue4 = min(numericData)-max(numericData)*0.01;
    maxValue4 = max(numericData)+max(numericData)*0.01;
end

function plotData(data, data2, data3, km, minValue, maxValue, minValue2, maxValue2, minValue3, maxValue3, minValue4, maxValue4, name)
    figure;
    % Первый подграфик
    subplot(4, 1, 2);
    plot(km, table2array(data(1, 2:end)), 'Color', 'b', 'LineWidth', 2);
    hold on;
    xlabel('Труба, км');
    ylabel('Разница давлений, Па');
    title('Профиль разницы давлений');
    xlim([0, 200]);
    ylim([minValue, maxValue]);

    % Второй подграфик
    subplot(4, 1, 1);
    plot(km, table2array(data2(1, 2:end), 'Color', 'b'));
    hold on;
    xlabel('Труба, км');
    ylabel('Давление, Па');
    title('Профиль давления');
    xlim([0, 200]);
    ylim([minValue2, maxValue2]);

    % Третий подграфик
    subplot(4, 1, 3);
    if strcmp(name, 'QuickWithQuasiStationaryModel') || strcmp(name, 'IdealQuickWithQuasiStationaryModel')|| strcmp(name, 'ChangeFlowQuickWithQuasiStationaryModel')
         % Создание нового массива km_interp с размером 2000
        km_interp = conv(km, [0.5, 0.5], 'valid');
        stairs(km_interp, table2array(data3(1, 2:end)), 'b', 'LineWidth', 2);

    else
        plot(km, table2array(data3(1, 2:end)), 'Color', 'b', LineWidth=2);
    end
    xlabel('Труба, км');
    ylabel('Плотность, кг/м3');
    title('Профиль плотности');
    xlim([0, 200]);
    ylim([minValue3, maxValue3]);

    % Четвертый подграфик
    % Четвертый подграфик
    subplot(4, 1, 4);
    % Подготовка данных
    % Преобразование времени в числовой формат
    t = table2array(data2(:, 1));  % Преобразование временных меток в массив
    t = (t - t(1));  % Преобразование времени в секунды от начала
    plot(t, table2array(data2(1:end, end)), 'Color', 'b', 'LineWidth', 2);
    hold on;
    plot(t(1), table2array(data2(1, end)), "Marker", ".", "LineStyle", "none", "MarkerSize", 20, "Color", 'r');
    hold on;
    text(t(1), table2array(data2(1, end)), ['t = ' string(data2{1,1}) ', с'], 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
    hold off;
    xlabel('Время, ч');
    ylabel('Давление, Па');
    title('Временной ряд давления на выходе');
    %xticks(0:12:84);  % устанавливаем метки с интервалом в 12 часов
    ylim([minValue4, maxValue4]);
    figure_size = [960, 540, 960, 540];
    set(gcf, 'Position', figure_size);
end

function createGif(data, data2, data3, km, minValue, maxValue, minValue2, maxValue2, minValue3, maxValue3, minValue4, maxValue4, name)
    % Получение кадра для первого кадра
    frame = getframe(gcf);
    [im, map] = rgb2ind(frame.cdata, 256, 'nodither');
    im(1, 1, 1, size(data, 1)) = 0;

    % Анимация
    for i = 2:(size(data(:,1))-2)
        % Отображение данных
        subplot(4, 1, 2);
        plot(km, table2array(data(i, 2:end)), 'Color', 'r', LineWidth=2);
        hold on;
        plot(km, table2array(data(1, 2:end)), 'Color', 'b', LineWidth=2);
        xlim([0, 200]);
        ylim([minValue, maxValue]);
        xlabel('Труба, км');
        ylabel('Разница давлений, Па');
        title('Профиль разницы давлений');
        hold off;

        subplot(4, 1, 1);
        plot(km, table2array(data2(i, 2:end)), 'Color', 'r');
        hold on;
        plot(km, table2array(data2(1, 2:end)), 'Color', 'b');
        xlim([0, 200]);
        ylim([minValue2, maxValue2]);
        xlabel('Труба, км');
        ylabel('Давление, Па');
        title('Профиль давления');
        hold off;

        subplot(4, 1, 3);
        if strcmp(name, 'QuickWithQuasiStationaryModel') || strcmp(name, 'IdealQuickWithQuasiStationaryModel')|| strcmp(name, 'ChangeFlowQuickWithQuasiStationaryModel')
            % Создание нового массива km_interp с размером 2000
            km_interp = conv(km, [0.5, 0.5], 'valid');
            stairs(km_interp, table2array(data3(i, 2:end)), 'b', 'LineWidth', 2);

        else
            plot(km, table2array(data3(i, 2:end)), 'Color', 'b', LineWidth=2);
        end
        xlabel('Труба, км');
        ylabel('Плотность, кг/м3');
        title('Профиль плотности');
        xlim([0, 200]);
        ylim([minValue3, maxValue3]);

        subplot(4, 1, 4);
        % Подготовка данных
        % Четвертый подграфик
        subplot(4, 1, 4);
        % Подготовка данных
        % Преобразование времени в числовой формат
        t = table2array(data2(:, 1));  % Преобразование временных меток в массив
        t = (t - t(1));  % Преобразование времени в секунды от начала
        plot(t, table2array(data2(1:end, end)), 'Color', 'b', 'LineWidth', 2);
        hold on;
        plot(t(i), table2array(data2(i, end)), "Marker", ".", "LineStyle", "none", "MarkerSize", 20, "Color", 'r');
        hold on;
        text(t(i), table2array(data2(i, end)), ['t = ' string(data2{i,1}) ', с'], 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
        hold off;
        xlabel('Время, ч');
        ylabel('Давление, Па');
        title('Временной ряд давления на выходе');
        ylim([minValue4, maxValue4]);
        figure_size = [960, 540, 960, 540];
        set(gcf, 'Position', figure_size);
        % Получение кадра
        frame = getframe(gcf);
        % Добавляем кадр к гифке
        im(:, :, 1, i) = rgb2ind(frame.cdata, map, 'nodither');
    end

    % Сохранение гифки в файл
    filename = strcat(name, '.gif');
    imwrite(im, map, filename, 'DelayTime', 0.02, 'LoopCount', inf);
    disp(['Гифка сохранена в файл: ' filename]);
end

function plotScaledData(data, data3, km, minValue3, maxValue3, name)
figure;
%i = 2423;
i = 1694;
%i = 1962;
%i = 2021;
%i = 2694;
        % Отображение данных
        subplot(2, 1, 1);
        plot(km, table2array(data(i, 2:end)), 'Color', 'r', LineWidth=2);
        hold on;
        plot(km, table2array(data(1, 2:end)), 'Color', 'b', LineWidth=2);
        xlim([192, 198]);
        %xlim([165, 175]);
        ylim([4*10000, 4.11*10000]);
        xlabel('Труба, км');
        ylabel('Разница давлений, Па');
        title('Профиль разницы давлений');
        grid on;
        hold off;


        subplot(2, 1, 2);
        if strcmp(name, 'QuickWithQuasiStationaryModel') || strcmp(name, 'IdealQuickWithQuasiStationaryModel')|| strcmp(name, 'ChangeFlowQuickWithQuasiStationaryModel')
            % Создание нового массива km_interp с размером 2000
            km_interp = conv(km, [0.5, 0.5], 'valid');
            stairs(km_interp, table2array(data3(i, 2:end)), 'b', 'LineWidth', 2);
        else
            plot(km, table2array(data3(i, 2:end)), 'Color', 'b', LineWidth=2);
        end
        xlabel('Труба, км');
        ylabel('Плотность, кг/м3');
        title('Профиль плотности');
        xlim([192, 198]);
        %xlim([165, 175]);
        ylim([minValue3, maxValue3]);
        grid on;
end

