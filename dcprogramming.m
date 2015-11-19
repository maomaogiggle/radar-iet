% ======================================================================%
% this program is used to calculate the lower bound of the S2C2
% author: Xiangrong Wang%
% date: 02/Dec/2014 %
% ======================================================================%

function [z0,scnr] = dcprogramming(Vc,Vt,N,K)

%define constant vector and matrix
%initialization
e = ones(N,1);
z0 = ones(N,1);
zeta = 0.01;
error = 1;
nummax = 10;
itenum = 1;

while (error > zeta && itenum < nummax)
    %begin the optimization procedure
    cvx_begin
    variable z(N,1);
    Dt = Vt'*diag(z)*Vt;
    alpha1 = log_det(Dt);
    G = Vc*((Vc'*diag(z0)*Vc)\Vc');
    grad = real(diag(G));
    alpha2 = log_det(Vc'*diag(z0)*Vc)+grad'*(z-z0);
    minimize (alpha2-alpha1);
    subject to 
        e'*z == K;
        0<= z <=1;
    cvx_end
    error = norm(z-z0);
    z0 = z;  
    itenum = itenum + 1;
end
   scnr = abs(det(Vt'*diag(z0)*Vt)/det(Vc'*diag(z0)*Vc));
end




