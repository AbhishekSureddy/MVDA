% Assignment --3
% Problem -1 -- Identification of model using PCA
clc;
clear all;
% loading the data
load('flowdata2.mat');

%% Part - (a) -- Applying PCA to Fmeas data matrix 
% Assuming initial intercept term to be zero
[u,s,v] = svds(Fmeas);
s_diag = diag(s);

% measured value
% Assuming "three" linear relationships
v_ind = 2;
Ameas = v(: , (v_ind+1):end)'; 
Am_dep = Ameas(:,[1,2,4]); % dependent variable coefficients
Am_ind = Ameas(:,[3,5]);   % independent variable coefficients

est_reg_coeff = -inv(Am_dep)*Am_ind;    % dependent variables interms of independent variables
                                        % here F3 and F5 are independent
% for true A
At_dep = Atrue(:,[1,2,4]); % dependent variable coefficients
At_ind = Atrue(:,[3,5]);   % independent variable coefficients

Z_dt = -inv(At_dep)*At_ind; % dependent variables interms of independent

disp('The eigenvalues are: ');
disp(s_diag.^2);

disp("maximum absolute difference = ")
disp(max(max(abs(est_reg_coeff-Z_dt))))

%% 1b --- IPCA,Estimation of error variances in case of 2 independent variables
flag = 1;
iter = 0;
nfact = 2;
sumsing = 0;
Z = Fmeas;
[nvar nsamples] = size(Z');
scale_factor = sqrt(nsamples)*ones(1,nvar);

while(flag)
    iter = iter + 1;
    Zs = zeros(size(Z));
    Zs = Z./scale_factor;% Equivalent to dividing by cholesky matrix(diagonal here)
    % Estimate number of PCs to retain
    [u s v] = svds(Zs);
    sdiag = diag(s);
    %sum of singular values of constraint eigen vectors
    sumsingnew = sum(sdiag(nfact+1:end));
    Amat_s = v(:,[nfact+1:nvar])';
    Amat_b   = Amat_s./scale_factor;
    if(abs(sumsingnew-sumsing)< 0.001)
        flag = 0;
    else
        est_std = std_est(Fmeas,Amat_b);
        scale_factor = sqrt(nsamples)*est_std';
        %scale_factor = est_std';
        sumsing = sumsingnew;
    end
end

est_reg_coef_ipca =  -inv(Amat_b(:,[1,2,4]))*Amat_b(:,[3,5]);
disp("For n = 2, independent variables or for 3 constraint equations");
disp("Estimated regression coefficients by IPCA are :")
disp(est_reg_coef_ipca)
disp('The estimated variances are: ');
disp(est_std);
disp('The eigenvalues are: ');
disp(sdiag.^2);

%% 1c --- IPCA,Estimation of error variances in case of 1 independent variables

flag = 1;
iter = 0;
nfact = 1;
sumsing = 0;
Z = Fmeas;
[nvar nsamples] = size(Z');
scale_factor = sqrt(nsamples)*ones(1,nvar);

while(flag)
    iter = iter + 1;
    Zs = zeros(size(Z));
    Zs = Z./scale_factor;% Equivalent to dividing by cholesky matrix(diagonal here)
    % Estimate number of PCs to retain
    [u s v] = svds(Zs);
    sdiag = diag(s);
    %sum of singular values of constraint eigen vectors
    sumsingnew = sum(sdiag(nfact+1:end));
    Amat_s = v(:,[nfact+1:nvar])';
    Amat   = Amat_s./scale_factor;
    if(abs(sumsingnew-sumsing)< 0.001)
        flag = 0;
    else
        est_std = std_est(Fmeas,Amat);
        scale_factor = sqrt(nsamples)*est_std';
        %scale_factor = est_std';
        sumsing = sumsingnew;
    end
end
disp("For n = 1 independent variable or for 4 constraint equations");
disp('The estimated variances are: ');
disp(est_std);
disp('The eigenvalues are: ');
disp(sdiag.^2);

%% prob -- 1d -- Finding the best and worst possible sets
var = [1 2 3 4 5];
Combinations = combnk(var,2);    % listing all the combinations of ind variables (taken 2 at a time) 
Det_matrix = [];
Max_Diff = []; 
for i = 1:length(Combinations)
    A1 = Amat_b;  % matrix calculated in part b of the question
    A2 = Atrue; % temporary matrix
    At_ind3 = [A2(:,[Combinations(i,1),Combinations(i,2)])];
    Ac_ind3 = [A1(:,[Combinations(i,1),Combinations(i,2)])];
    A1(:,[Combinations(i,1),Combinations(i,2)]) = [];
    A2(:,[Combinations(i,1),Combinations(i,2)]) = [];
    %calculating the determinant of the 3x3 matrix formed by coefficients of
    %dependent variables
    Det_matrix = [Det_matrix , det(A1)];
    estd_regress_coeff3 = -inv(A1)*Ac_ind3;
    true_regress_coeff3 = -inv(A2)*At_ind3;
    maxdiff2 = max(max(abs(true_regress_coeff3 - estd_regress_coeff3))); % max_diff calculation
    Max_Diff = [Max_Diff , maxdiff2]; 
 
end

min_max_diff = min(Max_Diff);        % minimum of max_diff values
n = 2;
B = zeros(n,1);% temp matrrix
index = zeros(n,1);
D1 = Det_matrix; 
for i=1:n
    [B(i),index(i)] = min(abs(D1));% storing the indices of minimum determinant 
    D1(index(i)) = 1000;% large number
end
disp('The Combinations of independent variables are: ');
disp("F"+string(Combinations));
disp('The determinant values for the matrix with dependent variable coefficients: ')
disp(Det_matrix);
disp('The max_diff for different combinations of independent variables considered: ');
disp(Max_Diff);
disp('The worst possible sets of independent variables are: ');
disp("F"+string(Combinations(index(1,1),:)));
disp("F"+string(Combinations(index(2,1),:)));
disp('The best set of independent variables are: ');
disp("F"+string(Combinations(find(Max_Diff == min_max_diff),:)));

