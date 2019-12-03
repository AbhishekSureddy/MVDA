clc;
clear all;
%% loading the data
load('Inorfull.mat');
Z = DATA(1:5:130,:);  % taking the first sample of the replicates from each mixture
C = CONC(1:5:130,:);

%% ipca

flag = 1;
iter = 0;
sumsing = 0;
[nvar,nsamples]=size(Z');
scale_factor = sqrt(nsamples)*ones(1,nvar);
nfact = 3;  % total number of independent variables
while(flag==1)
    
    iter = iter + 1
    Z_s = Z./scale_factor;
    % Estimate number of PCs to retain
    [u s v] = svd(Z_s); % use svd not svds because we need 173 columns of v
    sdiag = diag(s);
    sumsingnew = sum(sdiag(nfact+1:end));
    Amat_s = v(:,[nfact+1:nvar])';
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

