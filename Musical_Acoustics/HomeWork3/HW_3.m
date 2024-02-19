%% DATA AND ELECTRICAL EQUIVALENCE

%top plate
k_p=1.41e5;
m_p=0.128*0.385;
A_p=0.0375*0.385;
R_p=32;
%soundhole
m_h=0.000804;
A_h=0.00785;
R_h=30;
%air cavity
v_v=0.0172;
r_v=eps;
rho_air=1.204;
c=343;


M_p=m_p/(A_p^2);
C_p=(A_p^2)/k_p;
M_h=m_h/(A_h^2);
C_v=v_v/(rho_air*c^2);

Fs=44100;
time=10;
out=sim('HW_3.slx', time);
f_max=500;
input=squeeze(out.impulse_voltage.Data);
output=squeeze(out.current.Data);
f=[0:Fs/length(input):Fs-1/length(input)]';
omega=2*pi*f;
Z_bridge=fft(input)./fft(output);

%% BRIDGE IMPEDANCE
figure(1);
subplot(2, 1, 1);
plot(out.impulse_voltage.Time, squeeze(out.impulse_voltage.Data), 'LineWidth', 2);
ylabel('Voltage [V]');
xlabel('Time [s]');
title('Input voltage');
hold on;
subplot(2, 1, 2);
plot(out.current.Time, squeeze(out.current.Data), 'LineWidth', 2);
ylabel('Current [A]');
xlabel('Time [s]');
title('Output current');


figure(2);
subplot(2, 1, 1);
plot(f, db(abs(Z_bridge)), 'LineWidth', 2 );
ylabel('Magnitude [dB]');
xlabel('Frequency [Hz]');
title('$|Z(\omega)|$', interpreter='latex', FontSize=15);
xlim([1, f_max]);
subplot (2, 1, 2);
plot(f, angle(Z_bridge)*180/pi, 'LineWidth', 2 );
title('$\angle Z(\omega)$', interpreter='latex', FontSize=15);
ylabel('Angle [deg]');
xlabel('Frequency [Hz]');
xlim([1, f_max]);
ylim([-180, 180]);

%% ALTERNATIVE ANALITICAL BRIDGE IMPEDANCE
% Fs=44100;
% f_max=500;
% f=1:f_max;
% omega=2*pi*f;
% 
% Z_p=1i*omega*M_p+1./(1i*omega*C_p)+R_p;
% Z_v=1./(1i*omega*C_v);
% Z_h=1i*omega*M_h+R_h;
% 
% Z_bridge=Z_p+(Z_h.*Z_v)./(Z_h+Z_v);
% 
% figure(3);
% sgtitle('Bridge Impedance, $Z(\omega)=\frac{F}{v}$', interpreter='latex');
% subplot(2, 1, 1);
% plot(f, db(abs(Z_bridge)), 'LineWidth', 2 );
% title('$|Z(\omega)|=\frac{F}{v}$', interpreter='latex');
% xlim([0, f_max]);
% subplot (2, 1, 2);
% plot(f, angle(Z_bridge), 'LineWidth', 2 );
% title('$\angle Z(\omega)=\frac{F}{v}$', interpreter='latex');
% xlim([0, f_max]);
%% transfer function from the plucking point to the bridge

Z=Z_bridge./max(abs(Z_bridge));
T=1/Fs;
zeta=exp(1i*omega*T);
notes=[82.41, 110, 146.83, 196, 246.94, 329.63];
N_s=floor(Fs./(2*notes));
g=5;
beta=1/g;
%p=-1;
%gain=.85;
%Hz=gain*(zeta+1)./(zeta-(p));

for id=1:length(N_s)
        N_left=floor(beta*N_s(id));
        N_right=N_s(id)-N_left;
        R_f=0.99;
        R_b=0.99;
        H_E2R1=zeta.^(-N_left)*(-R_f).*zeta.^(-N_left).*zeta.^(-N_right);
        H_loop=(-R_b)*(-R_f)*zeta.^(-2*N_right).*zeta.^(-2*N_left);
        H_E1R1=zeta.^(-N_right);
        H_EB(id, :)=0.5*(1+H_E2R1).*(H_E1R1./(1-H_loop)).*Z./zeta*(1-(-R_b));
end

figure(4);
subplot(3, 2, 1);
plot(f, db(abs(H_EB(1, :))), 'LineWidth', 2 );
title('$f_0=82.41 \ Hz$', interpreter='latex', FontSize=13);
ylabel('Magnitude [dB]');
xlim([0, f_max]);

subplot(3, 2, 2);
plot(f, db(abs(H_EB(2, :))), 'LineWidth', 2 );
title('$f_0=110 \ Hz$', interpreter='latex', FontSize=13);
ylabel('Magnitude [dB]');
xlim([0, f_max]);

subplot(3, 2, 3);
plot(f, db(abs(H_EB(3, :))), 'LineWidth', 2 );
title('$f_0=146.83 \ Hz$', interpreter='latex', FontSize=13);
ylabel('Magnitude [dB]');
xlim([0, f_max]);

subplot(3, 2, 4);
plot(f, db(abs(H_EB(4, :))), 'LineWidth', 2 );
title('$f_0=196 \ Hz$', interpreter='latex', FontSize=13);
ylabel('Magnitude [dB]');
xlim([0, f_max]);

subplot(3, 2, 5);
plot(f, db(abs(H_EB(5, :))), 'LineWidth', 2 );
title('$f_0=246.94 \ Hz$', interpreter='latex', FontSize=13);
xlabel('Frequency [Hz]');
ylabel('Magnitude [dB]');
xlim([0, f_max]);

subplot(3, 2, 6);
plot(f, db(abs(H_EB(6, :))), 'LineWidth', 2 );
title('$f_0=329.63 \ Hz$', interpreter='latex', FontSize=13);
xlabel('Frequency [Hz]');
ylabel('Magnitude [dB]');
xlim([0, f_max]);

%% TIME RESPONSE
t=0:1/Fs:time;
impulse=zeros(length(notes), length(t));
L=0.65;
x=1/5*L;
h=0.003;
d=1/5*L;
y_t=zeros(length(notes), length(t));
dt=t(2)-t(1);

for i=1:length(notes)
    c(i)=2*L.*notes(i);
    for n=1:1000
        y_t(i, :)=y_t(i, :)+(2*h*L^2*sin(d*n*pi/L)*sin(n*pi*x/L)*cos(2*pi*n*c(i).*t/(2*L)))/(n^2*pi^2*d*(L-d));
    end
    y_1(i, :)=diff(y_t(i, :))/dt;
    amp(i)=min(y_1(i, :));
    impulse(i, 1)=amp(i);
end



t_d=t(1:end-1);

% y_2=diff(y_1)/dt;
% t_d1=t(1:end-2);







for id=1:length(notes)
    %H_EB(id, (floor(500 * length(H_EB)/Fs) + 1):length(H_EB))=0;
    imp(id, :)=fft(impulse(id, :));
    F(id, :)=imp(id, :).*H_EB(id, :);
    iF(id, :)=ifft(F(id, :));
    
    
    figure(5);
    hold on;
    subplot(3, 2, id);
    plot(t, real(iF(id, :)));
    title("$f_0=$"+notes(id)+"$\ Hz$", interpreter='Latex');
    xlim([0, 3]);
    xlabel('$Time \ [s]$','Interpreter','latex','FontSize',12);
    ylabel('$F(t) \ [N]$','Interpreter','latex','FontSize',12);

    % figure(6);
    % subplot(2, 1, 1);
    % plot(t, y_2);
    % grid on;
    % subplot(2, 1, 2);
    % plot(f, imp);
    % grid on;
    
    audiowrite("string"+id+".wav", real(iF(id, :))/max(abs(real(iF(id, :)))), Fs);
    
end