function [std] = std_est( Z , Amat )
%   Detailed explanation goes here
r = Z*Amat';
[m,n]=size(Amat);
if(n > m*(m+1)/2)
    disp("NOt possible to estimate the error variance")
    return
else 
    y = cov(r,1);  % covariance of residues
    y = y(:);    % vec(A*cov(e)*A')
    A = Amat;
    G = [];
    for j = 1:n
        C = [];
        for i = 1:m
            C = [C; A(i,j)*A(:,j)];
        end
    G = [G, C];
    end
    err_cov = pinv(G)* y ;
    std = sqrt(abs(err_cov));    %standard deviations of errors
end
end

