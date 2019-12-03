load("vpdata.mat");
% constants
A = 14.0568;
B = 2825.42;
C = 230.44;
sd_t = 0.18;
sd_p = 2;

% training first 70 samples
T_tr = (temp(1:70)-mean(temp(1:70)))/std(temp(1:70));
P_tr = psat(1:70);
% validation 
T_va = (temp(71:end)-mean(temp(1:70)))/std(temp(1:70));
P_va = psat(71:end);
w = 50;% width
% Kernel matrix creation
K =[];
for i=1:length(T_tr)
    k = T_tr-T_tr(i);
    kk = exp(-k.^2/w);
    K =[K,kk];
end

% eig vals and eig vecs

[vec,lamda] = eig(K);
lamda = diag(lamda);
disp(lamda);
lamda = lamda.^(-0.5);
n = 2;
v1 = vec(:,70-n:end);
lam1 = lamda(70-n:end);
lam1 = diag(lam1);
T1 = (lam1)*v1'*K;
Beta = inv(T1*T1')*T1*P_tr;

% prediction
pred = [];
for i = 1:length(T_va)
   k_p = exp(-(T_tr-T_va(i)).^2/w);
   t = (lam1)*v1'*k_p;
   pred_val = Beta'*t;
   pred = [pred;pred_val];
    
end

% accuracy

mse = 1/30*sum((pred-P_va).^2);
disp(['mse :' num2str(mse)]);
pred = [];

T_pred = [(55-mean(temp(1:70)))/std(temp(1:70));(100-mean(temp(1:70)))/std(temp(1:70))] ;
for i = 1:length(T_pred)
   k_p = exp(-(T_tr-T_pred(i)).^2/w);
   t = (lam1)*v1'*k_p;
   pred_val = Beta'*t;
   pred = [pred;pred_val];
    
end
disp(pred)






