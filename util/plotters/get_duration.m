%% Достает длину
function T = get_duration(data)
names = keys(data);
profiles = data(names{1}); % массив профилей для первой переменной группы
T = profiles.size();

end