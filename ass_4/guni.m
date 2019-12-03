% total number of species
m = 3; 
% structure of the mixing matrix as derived in part-(a)
Astruct = [1 1 0 1 1 1 0; 1 0 1 0 1 0 1;0 1 1 1 0 1 1]'; 
[u,s,v] = svd(measabs);          % Applying pCA
% Mixing matrix obtained in a rotated basis (not the exact basis of purespecies space) A0 = measabs*v(:,1:nspecies);    

% finding the non-singular rotation matrix M(all diagonal elements are 1)
% Procedure for finding M (explained in the report with an example)
M = [];
for i = 1:m
    Z = [];
    %finding the index of zero elements in each column of Astruct
    ind = find(Astruct(:,i) == 0);
        for j = 1:(length(ind))
            Z = [Z ; A0(ind(j),:)];
        end
            Y = - Z(:,i);
            Z(:,i) = [];
            M = [M , pinv(Z)*Y];
end
% reshaping the rotation matrix obtained into 3x3 matrix by putting 1's in the diagonal
d = diag([1 1 1]);
M1 = M(:);
r1 = d(:);
ind1 = find(r1 == 0);
for k = 1: length(ind1)
    r1(ind1(k),1) = M1(k,1);
end

M_rot = reshape(r1 , [3,3]); %final non_singular rotation matrix
disp('The rotation Matrix is: ' );
disp(M_rot);
%Actual mixing matrix in accordance with Astruct
A_PCAR = A0*M_rot;
disp('The MIXING MATRIX obtained with PCA-Rotation is: ');
disp(A_PCAR);
%denoised purespecies spectra taking A_PCAR mixing Matrix
est_pspec = pinv(A_PCAR)*measabs;  % OLS between A_PCAR and measabs
%correlation between estimated denoised purespecies spectra and actual purespecies spectra
corr_coef = corr(est_pspec' , pureabs');
disp('The Correlation between the estimated and true pure component spectra for all three species using PCA-Rotation is: ')
disp(corr_coef);
