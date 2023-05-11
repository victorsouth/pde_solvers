%%
function plot_groups(filename, path, clusters, plot_position, dt, L, time_mult)

[~,name,~]=fileparts(filename);
path = [path '/' name];
if not(isfolder(path))
    mkdir(path);
else
    delete([path '/*.png']);
end


% Зачитка csv файла
raw_data = readcell(filename);
[time, data, unitof, group2vars, lims, gridof] = parse_csv_raw(raw_data);
time = time * time_mult;


% Расчет координаты по csv
x = get_grid(L, data, gridof);

% отображение графиков

plot_data_at_grid(path, x, time, clusters, plot_position, dt, data, unitof, group2vars, lims);


end

%%
function plot_data_at_grid(path, x, time, clusters, plot_position, dt, data, unitof, group2vars, lims)

x_cell = x(1:end-1) + (x(2)-x(1))*0.5;
T = get_duration(data);

% Цикл по группам переменных, заданных в CSV

for c = clusters
    cluster = c{1};
    
    prev_time = -inf;
    % итерация по времени
    for t = 0:T - 1
        current_time = time(t+1);
        if current_time - prev_time < dt
            continue;
        end
        prev_time = current_time;
        time_str = ['t = ' num2str(current_time, '%f')];
        
        % формируем vars_str - строка с переменными и их размерностями
        vars_str = '';
        for group_index = 1:length(cluster)
            group = cluster(group_index);

            group_names = group2vars(group);
            for i = 1:length(group_names)
                unit = unitof(group_names{i});
                vars_str = [vars_str group_names{i} ', ' unit ' | '];
            end        
        end
        
        title_str = [vars_str time_str];
       

        % для данного момента времени строим график по всем группам по всем
        % переменным
        figure('Position', plot_position);  hold on;  grid on;   

        % пройтись по всем группам кластера
        for group_index = 1:length(cluster)
            group = cluster(group_index);

            group_names = group2vars(group);
            if length(cluster) > 1
                subplot(length(cluster), 1, group_index);
                hold on;  grid on; 
            end
            
            if group_index == 1
                title(title_str);
            end
            
            % для текущего момента времени перебираем все переменные
            for i = 1:length(group_names)
                profiles = data(group_names{i}); % массив профилей для переменной name
                profile = profiles.get(t);
                if length(x) == length(profile)
                    % это точки, т.е. границы ячеек
                    plot(x, profile, ...
                        'DisplayName', [group_names{i} ', ' unitof(group_names{i})]);
                elseif length(x) == length(profile) + 1
                    % это точки, т.е. границы ячеек
                    plot(x_cell, profile, ...
                        'DisplayName', [group_names{i} ', ' unitof(group_names{i})]);
                else
                    throw ('Profile pointcount incosisten with grid')
                end
                
            end
            legend;
            xlabel('x, м');
            
            [ymin, ymax] = get_group_limits(group_names, lims);
            if ymin < ymax
                ylim([ymin ymax]);
            end
            
        end
        

        title_str = replace(title_str, {'/', '|', '*', '^', '\'},'_');

        saveas(gcf, [path '/' title_str '.png'])
        close all;
    end            

    
end



end




