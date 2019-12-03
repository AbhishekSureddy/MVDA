function [ rmse ] = loocv( C_data, Z_data,flag )
%   Detailed explanation goes here
[nm,nw]= size(Z_data);                                  % nm = number of mixtures
                                                        % nw = number of wavelengths
if flag==1                                              % for OLS
    for i = 1:nm
        C = [C_data([1:i-1],:);C_data([i+1:nm],:)];     % C = concentration matrix
        Z = [Z_data([1:i-1],:);Z_data([i+1:nm],:)];     % Data matrix
        est_S = pinv(C)*Z;                              % From Z = C*S + E
        C_pred  = Z_data(i,:)*inv(est_S);               % predicting on the left out sample
        err = C_data(i,:)-C_pred;                       % error 
        rmse(i) = sqrt(sum(err.^2)/nw);                 % RMSE value of error
    end
else %for PCR
    for i = 1:nm
        C = [C_data([1:i-1],:);C_data([i+1:nm],:)];     % C = concentration matrix
        Z = [Z_data([1:i-1],:);Z_data([i+1:nm],:)];     % Data matrix
        [u,s,v] = svds(Z);                              % singular value decomposition
        for pcs = 1:5                                   % varying the number of pcs from 1 to 5
            T = Z*v(:,[1:pcs]);                         % scores matrix for the pcs considered
           % B = pinv(C)*T
            B = pinv(T)*C;                              % From C = T*B + E
            t_pred = Z_data(i,:)*v(:,[1:pcs]);
           % C_pred = t_pred*B'*inv(B*B');
            C_pred = t_pred*B;                          % predicting on left out sample
            err = C_data(i,:)-C_pred;                   % error
            rmse(i,pcs) = sqrt(sum(err.^2)/3);          % 3 in denominator because,
                                                        % no. of species in C matrix=3          
        end
    end
end
end

