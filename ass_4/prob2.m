%% Question 2
clc;
clear all ;
%% loading the data
load("yeastdata.mat");

Z = microarraydata;
% applying NCA
zind = find(~Astruct);
Astruct_mod = rand(size(Astruct));
Astruct_mod(zind) = 0;
[A_nca,P_nca,iter,ss] = gnca_fast(Z,Astruct_mod);

% for normalizing each column of A_nca
scale = [];
m = size(Astruct,2)
for i = 1:m
    scale = [scale,norm(A_nca(:,i))];
end
% normalized A_nca
A  =  A_nca ./ scale ;
P  =  P_nca .* scale';
var_P = var(P');
TFs   = [];
for i = 1:11
    [val,ind] = max(var_P);
    TFs       = [TFs,tfa(ind)];
    var_P(ind) = 0;
end
disp('The actual TFs implicated in cell cycle regulation are: ')
disp('Ace2, Fkh1, Fkh2, Mbp1, Mcm1, Ndd1, Skn7, Stb1, Swi4, Swi5, Swi6');
disp('The eleven TFs identified using the data are :');
disp(TFs);





