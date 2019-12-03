% Problem 2
clc;
clear all;
%% loading the data
load('Inorfull.mat');
wave_lengths = 300:2:650;

%% Prob 2a-- visualizing the wavelengths corresponding to max absorbance of pure species and applying OLS calibration model
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
%disp("==================================Prob 2a ======================================")
disp(" ")
disp("The max value of absorbance in Co, Cr, Ni occurs at : ")
disp(string(wave_lengths(lam_Co_max_ind))+" nm")
disp(string(wave_lengths(lam_Cr_max_ind))+" nm")
disp(string(wave_lengths(lam_Ni_max_ind))+" nm")

% Dealing with data
% taking the first sample of the replicates from each mixture
Z = DATA(1:5:130,[lam_Co_max_ind,lam_Cr_max_ind,lam_Ni_max_ind]);
C = CONC(1:5:130,:);
rmse = loocv(C,Z,1);        % applying LOOCV

disp("RMSE for 26 samples using LOOCV = ")
disp(rmse)
disp("The maximum value of RMSE obtained = "+ string(max(rmse)))


% ============================ End of prob2a ==================================================
%% prob 2b -- PCR
% taking the first sample of the replicates from each mixture
Z = DATA(1:5:130,:);
C = CONC(1:5:130,:);
% applying LOOCV, a function that is defined by me
rmse = loocv(C,Z,2);
%disp(" ")
%disp("==================================Prob 2b ======================================")
%disp(" ")
disp("RMSE for 26 samples using LOOCV obtained for choice of PCs from 1 to 5 = ")
disp(rmse)
disp("The maximum value of RMSE obtained obtained for choice of PCs from 1 to 5 = ")
disp(max(rmse))
disp("The mean value of RMSE obtained obtained for choice of PCs from 1 to 5 = ")
disp(mean(rmse))

%============================End of prob2b ==================================================
%% prob 2c -- mlpca
avg_std = mean(stdDATA);
std_mixture = std(DATA([1:5],:),1,2);
% scaling the data matrix by dividing with respective error standard
% deviations for different wavelengths
Z = Z./avg_std;
rmse = loocv(C,Z,2);              % Applying LOOCV
%disp(" ")
%disp("==================================Prob 2c ======================================")
%disp(" ")
figure;
plot(stdDATA(1,:));
title("plot showing standard deviations along wavelengths for first mixture")
disp("RMSE for 26 samples using LOOCV obtained for choice of PCs from 1 to 5 = ")
disp(rmse)
disp("The maximum value of RMSE obtained obtained for choice of PCs from 1 to 5 = ")
disp(max(rmse))
disp("The mean value of RMSE obtained obtained for choice of PCs from 1 to 5 = ")
disp(mean(rmse))

% End of prob2c

%% prob 2d -- ipca
Z = DATA(1:5:130,:);
flag = 1;
iter = 0;
sumsing = 0;        
[nvar,nsamples]=size(Z');
scale_factor = sqrt(nsamples)*ones(1,nvar);
nfact = 3;  % total number of independent variables
while(flag==1)
    
    iter = iter + 1;
    Z_s = Z./scale_factor;
    % Estimate number of PCs to retain
    [u s v] = svd(Z_s); % use svd not svds because we need 173 columns of v
    sdiag = diag(s);
    sumsingnew = sum(sdiag(nfact+1:end));
    Amat_s = v(:,[nfact+1:end])';
    Amat   = Amat_s./scale_factor;
    if(abs(sumsingnew-sumsing)< 0.001)
        flag = 0;
    else
        est_std = std_est(Z,Amat);
        scale_factor = sqrt(nsamples)*est_std';
        %scale_factor = est_std';
        sumsing = sumsingnew;
    end
end

Z = Z./est_std';
rmse = loocv(C,Z,2);
% disp(" ")
% disp("==================================Prob 2d ======================================")
% disp(" ")
disp("RMSE for 26 samples using LOOCV obtained for choice of PCs from 1 to 5 = ")
disp(rmse)
disp("The maximum value of RMSE obtained for choice of PCs from 1 to 5 = ")
disp(max(rmse))
disp("The mean value of RMSE obtained for choice of PCs from 1 to 5 = ")
disp(mean(rmse))
disp("The eigen values are ")
disp(sdiag.^2)




