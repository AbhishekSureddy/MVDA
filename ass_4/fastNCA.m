function [ A P ] = fastNCA(Z, Astruct, p )
%fastNCA Summary of this function goes here
%   Detailed explanation goes here
%     Z : Data matrix
%     Astruct : structural matrix
%     p : rank
[u,s,v] = svd(Z,0);
S1      = s(1:p,1:p);
W_star = u(:,1:p);%Z*v(:,1:p)*inv(S1);
Amix =[];
[ns,nvar] = size(Astruct);
for k = 1:nvar
    % step2 : rearrange rows of W_star
    [Wc Wr] = rearrange(W_star, Astruct, k);
    % step3 : consider Wr
    % step4 : S such that Q_r*S = 0
    [u1,s1,v1] = svd(Wr,0);
    S  = v1(:,end);
    % step5 : compute Wc*S  = a1*q1'*S + Ac*Qr*S = a1*q_t'
    %         q_t' = q1'*S  (1xn row vector)
    %         implies Wc*S has rank = 1
    [u2,s2,v2] = svds(Wc*S);
    a1 = u2(:,1);
    az = zeros(ns-length(a1),1);
    a1 =[a1;az];
    Amix =  [Amix,a1];
end 
A = reconstitute(Amix,Astruct);
P = pinv(A)*Z;
end

