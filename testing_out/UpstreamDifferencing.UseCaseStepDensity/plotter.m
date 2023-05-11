%%
function plotter()
clc
close all

L = 50000; % Длина линейного участка

dt = 0;

filename =  'output.csv';

plot_groups(filename, './', ...
     {[2]}, [50 50 1000 700], dt, L, 1);


end


