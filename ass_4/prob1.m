%% Question 1
clc;
clear all;
%% loading the data
load("ncadata.mat");
%% prob 1b
Z = measabs;
[u s v] = svds(Z);
% assuming the total number of species = 3
m =3;
scores = u(:,1:m)'*Z;
Z_est  = u(:,1:m)*scores;% denoised data
P = scores;
A = u(:,1:m);
% finding M (Rotational matrix) 
M = [];
Astruct = [1,1,0;1,0,1;0,1,1;1,0,1;1,1,0;1,0,1;0,1,1];

for i = 1:m
    zind  = find(~Astruct(:,i));
    nzind = find(Astruct(:,i));
    Az    = A(zind,:);
    [u,s,v] = svd(Az);
    mi    = v(:,m)*1/(v(i,m));
    M = [M,mi];
end
disp("Rotation matrix = ")
disp(M);
disp("The mixing matrix obtained with PCA-Rotation is :")
disp(A*M);
disp("Estimated pure component spectra = inv(M)*P")
disp("Correlation between estimated and true pure component spectra =")
disp(corr((inv(M)*P)',pureabs'))

%% prob 1c ----- Applying NCA to estimate pure-species spectra
zind = find(~Astruct);
Astruct = abs(rand(size(Astruct)));
Astruct(zind) = 0;
[A_nca,P_nca,iter,ss] = gnca_fast(Z,Astruct);
disp("A_nca = ")
disp(A_nca)
disp("correlation between estimated and true pure component spectra using nca is ")
disp(corr(P_nca',pureabs'))



