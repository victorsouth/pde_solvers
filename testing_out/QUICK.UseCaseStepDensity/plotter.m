%%
function plotter()
clc
close all

L = 50000; % Длина линейного участка

dt = 0;

files = ls('*.csv');

for i = 1:length(files)
    filename =  files(i, :);

    plot_groups(filename, './', ...
         {[2]}, [50 50 1000 700], dt, L, 1);
end


end


