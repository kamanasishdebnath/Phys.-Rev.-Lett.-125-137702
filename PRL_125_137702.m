% The following code can be used to generate the self stimulated echoes 
% reported in PHYSICAL REVIEW LETTERS 125, 137702 (2020)
% This is the main code. All figures can be generted by minor modifications

clear
clc


% Averaging of 100,000 ensembles with total spins= 10^10
num_of_ensemble=100000;
N= 1e10/num_of_ensemble;


% Resonance condition is considered hence cavity frequency = 0.0
wc= 0.0;
kappa= 2*pi*150*1000;   % linewidth of the cavity
gamma= 0;               % linewidth of the atoms
Gamma= 2*pi*0.5*1000;   % dephasing rate

% Normal distributed detunings for averaging
delta= normrnd(wc,1e7,[1,num_of_ensemble]);

% Coupling is homegeneous, g= 2pi * 8 HHz
g= zeros(1,num_of_ensemble) + 50.27;
number_of_equations= 1 + 2*(num_of_ensemble);

% Initial condition for the mean field equations
xinit= zeros(1,number_of_equations);
for i= 0:(num_of_ensemble-1)
    index= 2*i + 2;
    xinit(1,index)= 0;
    xinit(1,index+1)= -1;
end

% Specifying the simulation time
tstart= 0.0;             % start
tend=0.001;              % end
npoints= 10000;          % Number of points for plotting
tspan = linspace(tstart,tend,npoints);


% Intensity of the pulses
pulse1= 50000000000;
pulse2= pulse1;


% pi/2 pulse
T1= 0.0000002;
T2= 4.172e-7;
 
% Pi pulse separated by 45 \mu sec
separation= 45*1e-6;
T3= T2 + (separation);
T4= T3 + 2*(T2-T1);
options = odeset('RelTol',1e-4,'AbsTol',1e-6);


tic
detuning(i,:)= delta;
[t,y] = ode45(@(t,y) equations(t, y, g, N, wc, delta, kappa, Gamma, num_of_ensemble, pulse1, pulse2, T1, T2, T3, T4), tspan, xinit, options);
toc

photon= y(:,1).*conj(y(:,1));   % Photon number
real_field= real(y(:,1));       % Sigma X
imag_field= imag(y(:,1));       % Sigma Y
