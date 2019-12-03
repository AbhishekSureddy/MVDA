%clc;
clear all;
%% loading the data
load("ncadata.mat");
%% Applying the Fast-NCA algorithm
% variables
Z = measabs;
P = pureabs;
Astruct = [1,1,0;1,0,1;0,1,1;1,0,1;1,1,0;1,0,1;0,1,1];
disp("Structural matrix, Astruct = ")
disp(Astruct);
p = 3 ;
[A,P] = fastNCA(Z, Astruct, p );
disp('The Mixing Matrix estimated from FastNCA is: ');
disp(A);
corr_coef = corr(P' , pureabs');
disp('The Correlation between the estimated and true pure component spectra for all three species using FastNCA is = ')
disp(corr_coef);