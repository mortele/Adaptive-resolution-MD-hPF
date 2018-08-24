%% Clean up
close all;
clear variables;
clc;
format;


%% Load the systeminfo.out file
fileID  = fopen('../systeminfo.out', 'r');
data    = textscan(fileID, '%s %f');
data    = cellfun(@transpose, data, 'UniformOutput', false);
data    = {data{1}, num2cell(data{2})};
system  = cell2struct(data{2}, data{1}, 2);


%% Load .xyz file
fileID  = fopen('../positions.xyz', 'r');
data    = textscan( fileID, ...
                    '%d %f %f %f', ...              % Data format in file
                    system.number_of_particles, ... % Number of lines to read
                    'HeaderLines', 2);              % Number of lines to skip

% Organize cell array into vectors.
type = data{1};
x    = data{2};
y    = data{3};
z    = data{4};


%% Plot the particle positions
hFig = figure();
hPlot = scatter3(x, y, z, type .* 500, 'k.');