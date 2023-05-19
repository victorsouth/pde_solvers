%%
function plotter()
clc
close all

L = 50000; % Длина линейного участка

dt = 0;

files = dir('*.csv');

for i = length(files):-1:1
    filename =  files(i).name;

    disp(filename)
    plot_groups(filename, './', ...
          {[2]}, [50 50 1000 700], dt, L, 1);
end


end


