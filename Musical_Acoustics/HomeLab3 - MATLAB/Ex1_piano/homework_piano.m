%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Modeling of musical instruments homework.                               %
% Numerical simulation of piano strings                                   %
% Physical model for a struck string using finite difference.             %
%                                                                         %
% Musical Acoustics course                                                %
% Mirco Pezzoli                                                           %
% 2023                                                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all
clc

%% Setup
% define the string parameters and the simulation variables defined
% according to the provided values and the numerical implementation.
% We want to implement the finite difference scheme of a piano string tuned
% as C2.

% Temporal sampling parameters
Fs=44100;
T=1/Fs;
TimeLength = 8;
N=TimeLength*Fs;

% Fundamental note
f1=65.4;

% Boundary 
Zb = 1e3;      
Zl = 1e20; 

% String parameters
b1 = 0.5;
b2 = 6.25e-9;       
k = 7.5e-6;         
L = 1.92;           
Ms = 35e-3;       
rho = Ms/L;    
c=2*L*f1;
Ts= c^2*rho;
%Ts = 4*L^2*rho*f1^2;
%c= sqrt(Ts/rho);

% Spatial sampling parameters
% Aliasing condition
Fnyq = Fs/2;
gamma = Fnyq/f1;

% Number of maximum spatial steps
X_courant=c*T;
M_courant=L/X_courant;
M_chaigne = ((-1+(1+16*k*gamma^2)^(1/2))/(8*k))^(1/2);
X_chaigne = L/M_chaigne;

% Integer values
M = floor(M_chaigne);
X = L/M_chaigne;
lambda = c*T/X;  % Courant Number must be <1

% Spatial sampling
x = 0:L/(M-1):L;   %X*M=L

% FD parameters
mu = k^2/(c^2*X^2); 
nu = 2*b2*T/(X^2);

% Hammer parameters
MH = 4.9e-3;		 
bH = 1e-4;           
a = 0.12;            
Vh0 = 2.5;           

% Hammer contact window definition
w = 0.2;             
m0 = floor((a*L)/X);    %central striking spatial sample
w_s = floor(w/X);       %contact spatial window width in samples
g_m_m0 = zeros(1, M);   %initialization
g_m_m0(1, m0-w_s/2:m0+w_s/2-1) = hann(w_s);   %filling the array by centering a hanning window in m0


%PDE Coefficients:
a1 = -lambda^2*mu / (1+b1*T);
a2 = (lambda^2+4*mu*lambda^2+nu) / (1+b1*T);
a3 = (2-2*lambda^2-6*lambda^2*mu-2*nu) / (1+b1*T);
a4 = (-1+b1*T+2*nu) / (1+b1*T);
a5 = -nu / (1+b1*T);
aF = (T^2/rho) / (1+b1*T);

% Bridge boundary coefficients
bR1 = (2-2*lambda^2*mu-2*lambda^2) / (1+b1*T+Zb*lambda);
bR2 = (4*lambda^2*mu+2*lambda^2) / (1+b1*T+Zb*lambda);
bR3 = (-2*lambda^2*mu) / (1+b1*T+Zb*lambda);
bR4 = (-1+b1*T+Zb*lambda) / (1+b1*T+Zb*lambda);
bRF = ((T^2)/rho) / (1+b1*T+Zb*lambda);

% Left hand (hinged string end) boundary coefficients
bL1 = (2-2*lambda^2*mu-2*lambda^2) / (1+b1*T+Zl*lambda);
bL2 = (4*lambda^2*mu+2*lambda^2) / (1+b1*T+Zl*lambda);
bL3 = (-2*lambda^2*mu) / (1+b1*T+Zl*lambda);
bL4 = (-1+b1*T+Zl*lambda) / (1+b1*T+Zl*lambda);
bLF = ((T^2)/rho) / (1+b1*T+Zl*lambda);

%convenient coefficients for the hammer
d1 = 2 / (1+bH*T/(2*MH));
d2 = (-1+bH*T/(2*MH)) / (1+bH*T/(2*MH));
dF = (-T^2/MH) / (1+bH*T/(2*MH));

% Hammer felt parameters
p = 2.3;
K = 4e8;    

%% Computation of the FD scheme
% Initialization
Fh=zeros(1, N);
F = zeros(M, N);
eta=zeros(1, N);
y=zeros(M, N);
avg = zeros(1, N);

%special cases
%n=0;
eta(1, 1)=0;
Fh(1, 1)=0;

%n=1
eta(1, 2)=Vh0*T;
Fh(1, 2)=K.*abs(eta(1, 1)-y(m0, 1)).^p;

%n=2
for mi=2:M-1
    y(mi, 3)=y(mi-1, 2)+y(mi+1, 2)-y(mi, 1)+(T^2*M*Fh(1, 2).*g_m_m0(1, mi))/Ms;
end
eta(1, 3)=2*eta(1, 2)-eta(1, 1)-(T^2).*Fh(1, 2)/MH;
Fh(1, 3)=K*abs(eta(1, 3)-y(m0, 3))^p;
eta(1, 4)=eta(1, 3);


% Computation loop
for ni=4:N
    if eta(1, ni) < y(m0, ni)
        Fh(1, ni)=0;
    else
        Fh(1, ni)=K.*abs(eta(1, ni)-y(m0, ni)).^p;
    end

    F(:,ni) = Fh(1,ni).*g_m_m0(1,:);

    y(1, ni+1)=bL1*y(1, ni)+bL2*y(2, ni)+bL3*y(3, ni)+bL4*y(1, ni-1)+bLF*F(1, ni);
    y(2, ni+1)=a1*(y(4, ni)-y(2, ni)+2*y(1, ni))+a2*(y(3, ni)+y(1, ni)+a3*y(2, ni)+a4*y(2, ni-1)+a5*(y(3, ni-1)+y(1, ni-1))+aF*F(2, ni));
    for mi=3:M-2
        y(mi, ni+1)=a1*(y(mi+2, ni)+y(mi-2, ni))+a2*(y(mi+1, ni)+y(mi-1, ni))+a3*y(mi, ni)+a4*y(mi, ni-1)+a5*(y(mi+1, ni-1)+y(mi-1, ni-1))+aF*F(mi, ni);
    end
    y(M-1, ni+1)=a1*(2*y(M, ni)-y(M-1, ni)+y(M-3, ni))+a2*((y(M, ni)+y(M-2, ni))+a3*y(M-1, ni)+a4*y(M-1, ni-1)+a5*(y(M, ni-1)+y(M-2, ni-1))+aF*F(M-1, ni));
    y(M, ni+1)=bR1*y(M, ni)+bR2*y(M-1, ni)+bR3*y(M-2, ni)+bR4*y(M, ni-1)+bRF*F(M, ni);

    eta(1, ni+1)=d1*eta(1, ni)+d2*eta(1, ni-1)+dF*Fh(1, ni);

    for avgi = M-m0-6:M-m0+6
        avg(1,ni) = avg(1,ni) + y(avgi, ni);
    end
    avg(1,ni) = avg(1,ni)/12;
end

%% Plot the displacement in time

figure(1)

for ni = 1:50:size(y,2)
    plot(x, y(:,ni),'LineWidth',2);
    ylim([-5e-6, 5e-6]);
    xlim([0,L]);
    xlabel('$x[m]$', Interpreter='Latex');
    ylabel('$y[m]$', Interpreter='Latex');
    title(ni*T,'[s]');
    pause(0.000001);
end


% dist=25000;
% plots=[dist*1, dist*2, dist*3, dist*4, dist*5, dist*6];
% figure(4)
% subplot(3, 2, 1);
% plot(x, y(:,plots(1)),'LineWidth',2);
% ylim([-2.5e-6, 2.5e-6]);
% xlim([0,L]);
% xlabel('$x[m]$', Interpreter='Latex');
% ylabel('$y[m]$', Interpreter='Latex');
% title(plots(1)*T,'[s]');
% 
% subplot(3, 2, 2);
% plot(x, y(:,plots(2)),'LineWidth',2);
% ylim([-2.5e-6, 2.5e-6]);
% xlim([0,L]);
% xlabel('$x[m]$', Interpreter='Latex');
% ylabel('$y[m]$', Interpreter='Latex');
% title(plots(2)*T,'[s]');
% 
% subplot(3, 2, 3);
% plot(x, y(:,plots(3)),'LineWidth',2);
% ylim([-2.5e-6, 2.5e-6]);
% xlim([0,L]);
% xlabel('$x[m]$', Interpreter='Latex');
% ylabel('$y[m]$', Interpreter='Latex');
% title(plots(3)*T,'[s]');
% 
% subplot(3, 2, 4);
% plot(x, y(:,plots(4)),'LineWidth',2);
% ylim([-2.5e-6, 2.5e-6]);
% xlim([0,L]);
% xlabel('$x[m]$', Interpreter='Latex');
% ylabel('$y[m]$', Interpreter='Latex');
% title(plots(4)*T,'[s]');
% 
% subplot(3, 2, 5);
% plot(x, y(:,plots(5)),'LineWidth',2);
% ylim([-2.5e-6, 2.5e-6]);
% xlim([0,L]);
% xlabel('$x[m]$', Interpreter='Latex');
% ylabel('$y[m]$', Interpreter='Latex');
% title(plots(5)*T,'[s]');
% 
% subplot(3, 2, 6);
% plot(x, y(:,plots(6)),'LineWidth',2);
% ylim([-2.5e-6, 2.5e-6]);
% xlim([0,L]);
% xlabel('$x[m]$', Interpreter='Latex');
% ylabel('$y[m]$', Interpreter='Latex');
% title(plots(6)*T,'[s]');

%% Plot the synthesized signal play it and save it on the disk
t=0:T:TimeLength-T;      %temporal axis

figure(2);
plot(t, avg);
xlabel('$t[s]$', Interpreter='Latex');
ylabel('$y[m]$', Interpreter='Latex');
title('Average displacement');
grid on;

% Play the sound
out=avg./max(avg);
sound(out, Fs);

% Save on disk
audiowrite('10937966_DallaMora_piano.wav', out, Fs);
audiowrite('10883559_Murciano_piano.wav', out, Fs);

%Frequency analysis
spectrum = fft(avg,Fs);
freq = 0:Fs-1;
figure(3)
subplot(2,1,1)
plot(freq,db(abs(spectrum)), LineWidth=2);
xlabel("Frequency [Hz]", Interpreter='Latex');
xlim([0 500])
ylabel("Magnitude [dB]", Interpreter='Latex');
grid on;
subplot(2,1,2)
plot(freq,angle(spectrum)*180/pi, LineWidth=2);
xlim([0 500])
ylabel("Phase [deg]", Interpreter='Latex');
xlabel("Frequency [Hz]", Interpreter='Latex');
grid on;








