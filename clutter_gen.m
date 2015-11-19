function CMR = clutter_gen(radar,K)
lambda = radar.lambda;
Nphi = 720;
dphi = 2*pi/Nphi;
d = radar.d;
vp = radar.vp;
PRF = radar.PRF;
beta = radar.beta;
mv = radar.pulse;
nv = radar.pos';
M = radar.M;
N = radar.N;

%----------------------------------------------------
%generate the clutter for the required range
Amp = sqrt(radar.Pt/2); %No free space loss
CMR = zeros(M*N,K);
theta = radar.theta;
phiL = radar.phiL;

SM = zeros(M*N,Nphi);

for nphi = 1:Nphi
        phi = (nphi-1)*dphi;
        %Doppler related angle
        fr = cos(phi)*cos(theta);
        fd = (2*vp/lambda/PRF)*fr;
        %Spaital related angle
        fs = (d/lambda)*cos(theta)*cos(phi-beta);     
        s_t = exp(1i*2*pi*mv*fd);
        s_s = exp(1i*2*pi*nv*fs);
        S = s_s*s_t;
        %directivity
        D = 0.5*(1+cos(2*(phi-phiL)))*0.5*(1+cos(2*theta));
        G = abs(sinc(8*(phi-phiL)/pi))*abs(sinc(1*theta/pi));
        SM(:,nphi) = D*G*S(:)*dphi;
end

for k = 1:K
    Rf = Amp*(randn(Nphi,1)+1i*randn(Nphi,1));
    CMR(:,k) = SM*Rf;
end

