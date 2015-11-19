%This code is for randomization

function zd = randomization(Vc,Vt,z,L,N,K)
%Covariance
sigma = sqrt(z.*(1-z));
Sigma = diag(sigma);
%generate randomized selection vector with mean z and covariance Sigma
Zt = mvnrnd(z,Sigma,L);
%projection onto the 0-1 subspace
[Zs, Ind]= sort(Zt,2,'descend');
Z = zeros(L,N);
%Choose the one with maximum SCNR
SCNR = zeros(L,1);
for l = 1:L
    Z(l,Ind(l,1:K)) = 1;
    zt = transpose(Z(l,:));
    SCNR(l) = real(det(Vt'*diag(zt)*Vt)/det(Vc'*diag(zt)*Vc));
end

[Ym,Im] = max(SCNR);
zd = transpose(Z(Im,:));
