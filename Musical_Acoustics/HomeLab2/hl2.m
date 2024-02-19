%% DIAMETER
rho=1.204;
c=343;
f_0=300;

f_0=300;
d_1=0.04;
l_1=0.01;
L=l_1+0.93*(d_1/2);
S_1=pi*(d_1/2)^2;

max_freq=10000;
max_mesh_element_size=c/(max_freq*10);

ka=2*pi*f_0*d_1/(2*c);

D=(3/2*c^2*S_1/(pi^3*L*f_0^2))^(1/3);

mu=1.95e-5;
f=10;
delta=sqrt(mu/(pi*rho*f));
ciao=0.01*sqrt(2)/delta;

%% IMPEDANCE
impedance1234=readtable('impedance1234.csv');
freq4=table2array(impedance1234(1:1000,2));
P4=table2array(impedance1234(1:1000,3));
U4=table2array(impedance1234(1:1000,4));
freq1=table2array(impedance1234(1001:2000,2));
P1=table2array(impedance1234(1001:2000,3));
U1=table2array(impedance1234(1001:2000,4));
freq3=table2array(impedance1234(2001:3000,2));
P3=table2array(impedance1234(2001:3000,3));
U3=table2array(impedance1234(2001:3000,4));
freq8=table2array(impedance1234(3001:4000,2));
P8=table2array(impedance1234(3001:4000,3));
U8=table2array(impedance1234(3001:4000,4));

figure(1);
subplot(2, 1, 1)
y=db(abs(P4./-U4));
idx=freq4(islocalmin(y));
plot(freq4, y, 'LineWidth', 2);
xline(idx(1), '-r', idx(1) + " Hz");
grid on;
title('Input impedance, $d_1=4 \ cm$', 'Interpreter','latex', 'FontSize', 20);
ylabel('$$|Z_{in}|$$ [dB]', 'Interpreter','latex', 'FontSize', 20);
xlabel('Frequency [Hz]', 'FontSize', 12);

subplot(2, 1, 2)
plot(freq4, angle(P4./-U4), 'LineWidth', 2);
xline(idx(1), '-r', idx(1) + " Hz");
grid on;
ylabel('$$\angle Z_{in}$$ [rad]', 'Interpreter','latex', 'FontSize', 20);
xlabel('Frequency [Hz]', 'FontSize', 12);

figure(2);
subplot(2, 1, 1)
y=db(abs(P1./-U1));
idx=freq1(islocalmin(y));
plot(freq1, y, 'LineWidth', 2)
xline(idx(1), '-r', idx(1) + " Hz");
grid on;
title('Input impedance, $d_1=1 \ cm$', 'Interpreter','latex', 'FontSize', 20);
ylabel('$$|Z_{in}|$$ [dB]', 'Interpreter','latex', 'FontSize', 20);
xlabel('Frequency [Hz]', 'FontSize', 12);

subplot(2, 1, 2)
plot(freq1, angle(P1./-U1), 'LineWidth', 2)
xline(idx(1), '-r', idx(1) + " Hz");
grid on;
ylabel('$$\angle Z_{in}$$ [rad]', 'Interpreter','latex', 'FontSize', 20);
xlabel('Frequency [Hz]', 'FontSize', 12);

figure(3);
subplot(2, 1, 1)
y=db(abs(P3./-U3));
idx=freq3(islocalmin(y));
plot(freq3, y, 'LineWidth', 2)
xline(idx(1), '-r', idx(1) + " Hz");
grid on;
title('Input impedance, $d_1=3 \ cm$', 'Interpreter','latex', 'FontSize', 20);
ylabel('$$|Z_{in}|$$ [dB]', 'Interpreter','latex', 'FontSize', 20);
xlabel('Frequency [Hz]', 'FontSize', 12);

subplot(2, 1, 2)
plot(freq3, angle(P3./-U3), 'LineWidth', 2)
xline(idx(1), '-r', idx(1) + " Hz");
grid on;
ylabel('$$\angle Z_{in}$$ [rad]', 'Interpreter','latex', 'FontSize', 20);
xlabel('Frequency [Hz]', 'FontSize', 12);

figure(4);
subplot(2, 1, 1);
y=db(abs(P8./-U8));
idx=freq8(islocalmin(y));
plot(freq8, y, 'LineWidth', 2)
xline(idx(1), '-r', idx(1) + " Hz");
grid on;
title('Input impedance, $d_1=8 \ cm$', 'Interpreter','latex', 'FontSize', 20);
ylabel('$$|Z_{in}|$$ [dB]', 'Interpreter','latex', 'FontSize', 20);
xlabel('Frequency [Hz]', 'FontSize', 12);

subplot(2, 1, 2)
plot(freq8, angle(P8./-U8), 'LineWidth', 2)
xline(idx(1), '-r', idx(1) + " Hz");
grid on;
ylabel('$$\angle Z_{in}$$ [rad]', 'Interpreter','latex', 'FontSize', 20);
xlabel('Frequency [Hz]', 'FontSize', 12);

%% MODAL SUPERPOSITION DIFFERENT REPRESENTATION

modal=readtable('modal_sweep.csv');
d = table2array(modal(:,1));
m = table2array(modal(:,2));
diameter=unique(d, 'stable');
mode = unique(m);
freq = unique(table2array(modal(:,3)));
for j=1:numel(diameter)
    figure(j)
    sgtitle(diameter(j) + " cm", 'FontWeight', 'bold', 'FontSize',20);
    for i=1:numel(mode)
        P(:,i) = table2array(modal(m==mode(i) & d==diameter(j), 4));
        U(:,i) = table2array(modal(m==mode(i) & d==diameter(j), 5));
        plot(freq, db(abs(P(:,1)./U(:,i))), ':', 'LineWidth', 1.5);
        ylabel('$$|Z_{in}|$$ [dB]','Interpreter','latex','FontSize',20)
        xlabel('Frequency [Hz]','Interpreter','latex')
        hold on
        grid on 
    end
    summ(j, :)=db(abs(sum(P(:,1)./U,2)));
    plot(freq, summ(j, :), '-', 'LineWidth', 1.5, 'color', 'b')
    legend("Mode 0","Mode 1", "Mode 2","Mode 3",'sum');
    
    figure(10);
    plot(freq, summ(j, :), '-', 'LineWidth', 1);
    hold on;
    legend("4","1", "3","8");
end
%% ELECTRIC ANALOGUE
rho=1.204;
c=343;
f_0=300;
Fs=2*max_freq;

d_1=0.01;
l_1=0.01;
S_1=pi*(d_1/2)^2;

d_2=0.04;
l_2=2*l_1;
S_2=pi*(d_2/2)^2;

L=l_1+0.93*(d_1)/2;
D2=0.1586;
%1.8 alpha per far coincidere le freq tra comsol e simulink
L2=l_2+1.7*(d_2/2); 

Ma1=rho*L/S_1;
Ma2=rho*L2/S_2;
V1=(D/2)^3*pi*4/3;
Ca1=V1/(rho*c^2);

V2=(D2/2)^3*pi*4/3;
Ca2=V2/(rho*c^2);

R_1=rho*c/S_1;
R_2=rho*c/S_2;

figure(1);
subplot(2, 1, 1);
plot(out.simout.Time, squeeze(out.simout.Data), 'LineWidth', 2);
title('Input voltage');
ylabel('Voltage [V]','Interpreter','latex','FontSize',15)
xlabel('Time [s]','Interpreter','latex')
hold on;
subplot(2, 1, 2);
plot(out.simout1.Time, squeeze(out.simout1.Data), 'LineWidth', 2);
title('Output current');
ylabel('Current [A]','Interpreter','latex','FontSize',15)
xlabel('Time [s]','Interpreter','latex')

in=squeeze(out.simout.Data);
output=squeeze(out.simout1.Data);
f=[0:Fs/length(in):Fs-1/length(in)]';
H=fft(in)./fft(output);

figure(2);
idx=f(islocalmin(db(abs(H))));
plot(f, db(abs(H)), 'LineWidth', 2 );
grid on;
ylabel('$$|Z_{in}|$$ [dB]', 'Interpreter','latex', 'FontSize', 20);
xlabel('Frequency [Hz]', 'FontSize', 12);
xline(idx(1), '-r', idx(1) + " Hz");
xline(idx(2), '-r', idx(2) + " Hz");
xlim([0, Fs/2]);
%xlim([0, 500]);

%% Double comsol

impedance=readtable('double.csv');
freq=table2array(impedance(:,1));
P=table2array(impedance(:,4));
U=table2array(impedance(:,5));

figure(1);
%subplot(2, 1, 1)
y=db(abs(P./-U));
idx=freq(islocalmin(y));
plot(freq, y, 'LineWidth', 2);
xline(idx(1), '-r', idx(1) + " Hz");
xline(idx(2), '-r', idx(2) + " Hz");
grid on;
%title('Input impedance, $d_1=4 \ cm$', 'Interpreter','latex', 'FontSize', 20);
ylabel('$$|Z_{in}|$$ [dB]', 'Interpreter','latex', 'FontSize', 20);
xlabel('Frequency [Hz]', 'FontSize', 12);
xlim([0, 500]);

% subplot(2, 1, 2)
% plot(freq, angle(P./-U), 'LineWidth', 2);
% xline(idx(1), '-r', idx(1) + " Hz");
% grid on;
% ylabel('$$\angle Z_{in}$$ [rad]', 'Interpreter','latex', 'FontSize', 20);
% xlabel('Frequency [Hz]', 'FontSize', 12);
% xlim([0, 1000]);
