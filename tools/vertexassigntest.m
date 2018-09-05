%% Clean up
close all;
clear variables;
clc;
format;


%% Parameters
l = 5.27597;


%% 
nodes =  [  [0,0,0] 
            [1,0,0] 
            [0,1,0]
            [0,0,1] 
            [1,1,0] 
            [1,0,1] 
            [0,1,1] 
            [1,1,1]  ] .* l;
distances = zeros(1,8);
fractions = zeros(1,8);

p = rand(1,3) .* l;
L = 0;

for i=1:8
    distances(i) = norm(nodes(i,:) - p);
    L = L + distances(i)^2;
end
distances;
% L - (8*(p(1)^2 + p(2)^2 + p(3)^2) + 12*l^2 - 8*l * (p(1) + p(2) + p(3)))

% Computing vertex fractions
fractions(1) = (l-p(1)) * (l-p(2)) * (l-p(3)) / l^3;  % (0,0,0)
fractions(2) = ( p(1) ) * (l-p(2)) * (l-p(3)) / l^3;  % (1,0,0)
fractions(3) = (l-p(1)) * ( p(2) ) * (l-p(3)) / l^3;  % (0,1,0)
fractions(4) = (l-p(1)) * (l-p(2)) * ( p(3) ) / l^3;  % (0,0,1)

fractions(5) = ( p(1) ) * ( p(2) ) * (l-p(3)) / l^3;  % (1,1,0)
fractions(6) = ( p(1) ) * (l-p(2)) * ( p(3) ) / l^3;  % (1,0,1)
fractions(7) = (l-p(1)) * ( p(2) ) * ( p(3) ) / l^3;  % (0,1,1)
fractions(8) = ( p(1) ) * ( p(2) ) * ( p(3) ) / l^3;  % (1,1,1)

p
sum(fractions)