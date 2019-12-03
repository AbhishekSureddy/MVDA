clc;
clear all;
%% loading the data
load('Inorfull.mat');
wave_lengths = 300:2:650;
%% ploting the data and visualizing the max absorbance
figure;
plot(wave_lengths,PureCo,'linewidth',2,'color','r');
hold on;
plot(wave_lengths,PureCr,'linewidth',2,'color','g');
plot(wave_lengths,PureNi,'linewidth',2,'color','b');
legend('Co','Cr','Ni');
xlabel("wave lengths");
ylabel("Absorption");
title("max absorption wave lengths for Co, Cr, Ni ")
[val,lam_Co_max_ind] = max(PureCo);
[val,lam_Cr_max_ind] = max(PureCr);
[val,lam_Ni_max_ind] = max(PureNi);
disp("The max value of absorbance in Co, Cr, Ni occurs at : ")
disp(wave_lengths(lam_Co_max_ind))
disp(wave_lengths(lam_Cr_max_ind))
disp(wave_lengths(lam_Ni_max_ind))

%% Dealing with data
Z = DATA(1:5:130,[lam_Co_max_ind,lam_Cr_max_ind,lam_Ni_max_ind])  % taking the first sample of the replicates from each mixture
C = CONC(1:5:130,:);
rmse = loocv(C,Z,1);
