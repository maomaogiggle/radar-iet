% ======================================================================%
% this program is MCM for multiple interferences case
% author: Xiangrong Wang%
% date: 21/06/2014 %
% ======================================================================%

function [z0,opt_scnr] = mcm(Vc,Vt,N,K)
z0 = ones(N,1);
z = ones(N,1);
index = 1:N;

for k = 1:N-K
    alpha = zeros(N+1-k,1);
   
    if (k==1)
        for i = 1:N
            z0(i) = 0;
            alpha(i) = det(Vt'*diag(z0)*Vt)/det(Vc'*diag(z0)*Vc);
            z0(i) = 1;
        end  
        [Y,I] = sort(real(alpha),'descend');
        z(I(1)) = 0;
        z0 =z; 
    else
        indexx = index(logical(z));
        for i = 1:length(indexx)
            z0(indexx(i)) = 0;
            alpha(i) = det(Vt'*diag(z0)*Vt)/det(Vc'*diag(z0)*Vc);
            z0(indexx(i)) = 1;
        end
        [Y,I] = sort(real(alpha),'descend');
        z(indexx(I(1))) = 0;
        z0 = z;
    end
end
opt_scnr = Y(1);
end


