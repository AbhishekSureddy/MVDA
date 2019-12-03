clc;
clear all;
%% loading the data
load('Inorfull.mat');
Z = DATA(1:5:130,:);  % taking the first sample of the replicates from each mixture
C = CONC(1:5:130,:);

%% visualizing the plots of standard deviations
for i = 1:26
std_wave_length(i,:) = std(DATA([5*i-4:5*i],:));
end
avg_std = mean(std_wave_length);
std_mixture = std(DATA([1:5],:),1,2);
% scaling the data matrix by dividing with respective error standard
% deviations for different wavelengths
Z = Z./avg_std;
rmse = loocv(C,Z,2);
