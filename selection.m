%validation of MCM method
%3/Dec/2014
%Xiangrong

clc;clear;

radar = radar_init;
%inter-element space
d = radar.d;
%PRF
PRF = radar.PRF;
%platform velocity
vp = radar.vp;
%antenna index
nvs = radar.pos';
%pulse index
nvt = radar.pulse;
%antenna number
N = radar.N;
%pulse number
M = radar.M;
%wavelength
lambda = radar.lambda;
%platform angle
beta = radar.beta;
%number of DoFs
MN = M*N;
%number of selected antennas
K = 10;

%generating the clutter
L = 10000;
CMR = clutter_gen(radar,L);
C_f = (CMR*CMR')/L;
%clutter power
Pc = trace(C_f)/MN;
%noise power
Pn = Pc*(10^(-radar.CNR/10));
%calculate covariance matrix
C = (CMR*CMR')/L+Pn*eye(MN);

%eigenvalue decomposition
[V,D] = eig(C);
Vc = V(:,14:end);

%target
fs_t = 0;
as_t = exp(1i*2*pi*fs_t*nvs);
fd_t = [-0.5:0.01:0.5];

% calculate the selection vector by enumerate
index = 1:MN;
index_sub = nchoosek(index,K);
num = nchoosek(MN,K);
SCNR = zeros(1,num);

scnr_t = zeros(length(fd_t),1);
Z_t = zeros(MN,length(fd_t));
scnr_m = zeros(length(fd_t),1);
Z_m = zeros(MN,length(fd_t));
scnr_p = zeros(length(fd_t),1);
Z_p = zeros(MN,length(fd_t));
scnr_d = zeros(length(fd_t),1);
Z_d = zeros(MN,length(fd_t));


for i=1:length(fd_t)    
    at_t = exp(1i*2*pi*nvt*fd_t(i));
    A_t = as_t*at_t;
    t = reshape(A_t,MN,1); 
    Vt =[t,Vc]; 
    
   for q = 1:num
        %determinant objective
        z = zeros(MN,1);
        z(index_sub(q,:)) = 1;
        SCNR(q) = det(Vt'*diag(z)*Vt)/det(Vc'*diag(z)*Vc);
    end
    
    [scnr_t(i),I] = max(abs(SCNR));
    index_opt = index_sub(I,:);
    z = zeros(MN,1);
    z(index_opt) = 1;
    Z_t(:,i) = z;
    
    %MCM
    [Z_m(:,i),scnr_m(i)] = mcm(Vc,Vt,MN,K);  
    %DC programming
    [Z_p(:,i),scnr_p(i)] = dcprogramming(Vc,Vt,MN,K); 
    %randomization
    zd = randomization(Vc,Vt,Z_p(:,i),100,MN,K);
    scnr_d(i) = abs(det(Vt'*diag(zd)*Vt)/det(Vc'*diag(zd)*Vc));
    Z_d(:,i) = zd;    
    %lower bound
%     [Z_b(:,i),scnr_b(i)] = bound(Vc,Vt,MN,K);

end

figure;
plot(fd_t,scnr_t);
hold on;
plot(fd_t,scnr_m,'k');
hold on;
plot(fd_t,scnr_p,'r');
hold on;
plot(fd_t,scnr_d,'c');

figure;
plot(fd_t,(scnr_t-scnr_m)./scnr_t,'k');
hold on;
plot(fd_t,(scnr_t-scnr_d)./scnr_t,'c');




