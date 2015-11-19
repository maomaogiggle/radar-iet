function radar = radar_init
%Initalise parameters
c = 3e8;                        %Speed of light
fc= 1.488e9;                      %The Carrier frequency
lambda = c/fc;                 %The carrier wavelength  
PRF = 1984;                    %The Pulse Repetition Frequency
vp = 100.3;                       %Platform velocity (m/s)
H = 3e3;                        %Platform altitude in m
Rs = 12e3;                      %Start range gate
Rf = 24.8e3;                    %End Range gate
dR = 100;                       %Range resolution
phiL = 90*pi/180;                %Azimuth of Look direction in radians
theta = 15*pi/180;              %Elevation of Look direction
beta = 0*pi/180;             %platform angle in degrees (0 - broadside, 90 - forward looking);
                                      %Flight path is horizontal, parallel to the ground in the x-direction
N = 4;                          %Number of antenna elements
d = lambda/2;                  %Sensor Spacing in m
CNR = 20;                      %Clutter to noise ratio in dB
M = 4;                            %Number of temporal pulses
pulse = 0:M-1;                 %selected pulses
pos   = 0:N-1;                 %selected antennas
Pt = 1;                              %Transmitted Power
                                        %Clutter model used is that given by Klemm in his book (page 56)
                                        %The doppler frequency of the individual clutter arrival is f_D = 2vpcos(phi)cos(theta)/lambda
%----------------------------------------------
%Set up the radar details
radar.d = d;
radar.PRF = PRF;
radar.vp = vp;
radar.H = H;
radar.Rs = Rs;
radar.Rf = Rf;
radar.dR = dR;
radar.pulse = pulse;
radar.pos = pos;
radar.lambda = lambda;
radar.beta = beta;
radar.N = N;
radar.M = M;
radar.phiL = phiL;
radar.theta = theta;
radar.Pt = Pt;
radar.CNR = CNR;
%--------------------------------------------
