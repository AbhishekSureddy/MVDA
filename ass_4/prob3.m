clc;
clear all;
% loading the data
load('Inorfull.mat');
A = CONC;
P = [PureCo;PureCr;PureNi];
Z = DATA;

% modifying Z
ind = find(Z<0); % indices of negative elements
Z(ind) = 0;      % clipping the negative elements to zero

%% prob a

% data 
Z1 = Z(1:5:130,:);
[u,s,v] = svds(Z1);
% using the absolute values of loadings (v) and scores 
% as mentioned in the question.
A_init = abs(Z1*v(:,1:3));
P_init = abs(v(:,1:3))';

[A1,P1] = nmf(Z1,A_init,P_init,1e-10,100,10000);

cor_coef = corr(P1',P');
disp("The correlation coefficients are ")
disp(cor_coef)

%% prob b
Z_avg = [];
for i=1:5:130
    temp = [mean(DATA(i:i+4,:))];
    %disp(temp)
    Z_avg=[Z_avg;temp];
end
Z_avg(Z_avg < 0) = 0;
[u,s,v] = svds(Z_avg);
% using the absolute values of loadings (v) and scores 
% as mentioned in the question.
A_init_avg = abs(Z_avg*v(:,1:3));
P_init_avg = abs(v(:,1:3))';

[A_avg,P_avg] = nmf(Z_avg,A_init_avg,P_init_avg,1e-10,100,10000);

cor_coef = corr(P_avg',P');
disp("The correlation coefficients are ")
disp(cor_coef)
% we are able to resolve the ambiguity