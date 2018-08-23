%% Clean up
close all;
clear variables;
clc;
format;


%% Load the systeminfo.out file
fileID  = fopen('../systeminfo.out', 'r');
data    = textscan(fileID, '%s %f');



%% Load .xyz file
fileID  = fopen('../positions.xyz', 'r');
data    = textscan(fileID, '%d %f %f %f', 'HeaderLines', 2);

% Organize cell array into vectors.
type = data{1};
x    = data{2};
y    = data{3};
z    = data{4};


%% Plot the particle positions
hFig = figure();
hPlot = scatter3(x, y, z, type .* 500, 'k.');