load("arx.mat");
%% first order (a - part)
X_ols = [ymeas(1:999)',umeas(1:999)'];
Y_ols = ymeas(2:end)';
params = pinv(X_ols)*Y_ols;
disp("OLS");
disp(params);

%% TLS - (b - part)

Z = [X_ols,Y_ols];
[u s v] = svd(Z,0);
params = v(:,3)/v(3,3);
disp("TLS");
disp(params)

%% order 10 (d - part)
m =5;
Zu = [];
Zy = [];
for i = 1:m
   tempu = umeas(i:1000-m+i-1)';
   tempy = ymeas(i:1000-m+i-1)';
   Zu = [Zu,tempu];
   Zy = [Zy,tempy];
   
end
ymat = ymeas(m+1:1000)';
Z = [ymat,Zy,Zu];
Zs = Z - mean(Z);
[u s v] = svd(Zs,0);
params = v(:,end)/v(1,end);
disp(params);

%% DIPCA
%% 1b --- IPCA,Estimation of error variances in case of  independent variables
flag = 1;
iter = 0;
nfact = 4;
sumsing = 0;
[nvar nsamples] = size(Zs');
scale_factor = sqrt(nsamples)*ones(1,nvar);

while(flag)
    iter = iter + 1;
    Zss = zeros(size(Z));
    Zss = Zs./scale_factor;% Equivalent to dividing by cholesky matrix(diagonal here)
    % Estimate number of PCs to retain
    [u s v] = svd(Zss,0);
    sdiag = diag(s);
    %sum of singular values of constraint eigen vectors
    sumsingnew = sum(sdiag(nfact+1:end));
    Amat_s = v(:,[nfact+1:nvar])';
    Amat_b   = Amat_s./scale_factor;
    if(abs(sumsingnew-sumsing)< 0.001)
        flag = 0;
    else
        est_std = std_est(Zs,Amat_b);
        scale_factor = sqrt(nsamples)*est_std';
        %scale_factor = est_std';
        sumsing = sumsingnew;
    end
end
disp('The estimated variances are: ');
disp(est_std);
disp('The eigenvalues are: ');
disp(sdiag.^2);








