clc;
clear all;
%% loading the data
load('Inorfull.mat');
Z = DATA(1:5:130,:);  % taking the first sample of the replicates from each mixture
C = CONC(1:5:130,:);
%% finding the scores matrix "T"
rmse = loocv(C,Z,2);



