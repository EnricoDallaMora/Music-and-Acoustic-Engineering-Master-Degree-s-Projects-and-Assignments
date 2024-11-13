t = 0:0.01:4;
% alpha=2;
% omega=40;
% phi=1;
% x=exp(-alpha.*t);
% x1=exp(-alpha.*t).*cos(omega.*t+phi);
% figure(1);
% plot(t,x,'r','LineWidth',1.2);
% hold on;
% plot(t,x1,'b','LineWidth',1.2);
% legend('$e^{-\alpha t}$','$e^{-\alpha t}(\tilde{A}_1 e^{j\omega_d t} \tilde{A}_2 e^{-j\omega_d t})$', interpreter="latex")
% ylim([-2 2]);
% grid on;
% xlabel('Time','Fontsize',12); 
% ylabel('Amplitude','Fontsize',12)
% %title('Free motion with adimensional damping ratio h','Fontsize',18)

m=0.1;
k=2.53*10^4;

w_0=sqrt(k/m);
f_0=w_0/(2*pi); 
T=2*pi/w_0;

t1=0.576;
tau=-t1/log(10^(-1/2));
alpha=1/tau;
Q=w_0*tau/2;
R=2*m/tau;
bandwidth=2*1/tau;

xi=1/(w_0*tau);
xi1=R/(2*sqrt(k*m));
xi2=alpha/(sqrt(alpha^2+w_0^2));

w = 0:1:1000;
%adm = zeros(length(w));
Z = zeros(length(w));
Y = zeros(length(w));
for ii = 1:length(w)
    %adm(ii) = (1i*w(ii))/ ( w(ii)^2 - 2*1i*alpha*w(ii) - w_0^2) ;

    Z(ii) = R + 1i*(w(ii)*m - k/(w(ii)));
    Y(ii) = 1 / Z(ii);

end
figure(1)
plot(w, abs(Y), 'b', LineWidth=1);
xlabel('$\omega$','Fontsize',12, interpreter="latex"); 
ylabel('$|Y(\omega)|$','Fontsize',12, interpreter="latex");
title('Module of the Admittance $Y(\omega)$','Fontsize',18, interpreter="latex");

w_d=sqrt(w_0^2-alpha^2)

%% EXTERNAL FORCE

t = 0:0.0001:0.5;
f = [60, 80, 100, 120, 140, 160]';
omega = f*2*pi;

for jj = 1:length(f)
    index(jj) = find( abs(w-omega(jj))==min(abs(w-omega(jj))));
end

x = zeros(length(f), length(t));

for ii = 1:length(t)
    for jj = 1:length(f)
        v(jj)=0.1/abs(Z(index(jj)));

        phi(jj)=atan( ...
            (       v(jj)*cos(angle(Z(index(jj))))    -alpha.*((-v(jj)/index(jj))*sin(angle(Z(index(jj)))))         ) ...
            / ...
            (              w_d.*((-v(jj)/index(jj))*sin(angle(Z(index(jj)))))                ) ...
            );









        A(jj)=(-0.1*sin(angle(Z(index(jj)))))/(cos(phi(jj))*abs(Z(index(jj)))*index(jj));

        x(jj,ii) = exp(-alpha.*t(ii))*A(jj)*cos(w_d.*t(ii)+phi(jj)) + 0.1 * sin(2*pi*f(jj)*t(ii) + angle(Z(index(jj)))) / (2*pi*f(jj)*abs(Z(index(jj))));
    end
end


figure(2)
sgtitle('Complete time response for different values of excitation force frequency','Fontsize',18);

subplot(3, 2, 1)
plot(t, x(1,:), LineWidth=1.5)
xlabel('time [s]','Fontsize',12); 
ylabel('x(t) [m]','Fontsize',12);
title('$f_1=60 Hz$','Fontsize',16, interpreter="latex");



subplot(3, 2, 2)
plot(t, x(2,:), LineWidth=1.5)
xlabel('time [s]','Fontsize',12); 
ylabel('x(t) [m]','Fontsize',12);
title('$f_2=80 Hz$','Fontsize',16, interpreter="latex");


subplot(3, 2, 3)
plot(t, x(3,:), LineWidth=1.5)
xlabel('time [s]','Fontsize',12); 
ylabel('x(t) [m]','Fontsize',12);
title('$f_3=100 Hz$','Fontsize',16, interpreter="latex");


subplot(3, 2, 4)
plot(t, x(4,:), LineWidth=1.5)
xlabel('time [s]','Fontsize',12); 
ylabel('x(t) [m]','Fontsize',12);
title('$f_4=120 Hz$','Fontsize',16, interpreter="latex");


subplot(3, 2, 5)
plot(t, x(5,:), LineWidth=1.5)
xlabel('time [s]','Fontsize',12); 
ylabel('x(t) [m]','Fontsize',12);
title('$f_5=140 Hz$','Fontsize',16, interpreter="latex");


subplot(3, 2, 6)
plot(t, x(6,:), LineWidth=1.5)
xlabel('time [s]','Fontsize',12); 
ylabel('x(t) [m]','Fontsize',12);
title('$f_6=160 Hz$','Fontsize',16, interpreter="latex");



figure(3)
xlabel('time [s]','Fontsize',12); 
ylabel('x(t) [m]','Fontsize',12);
hold on
plot(t, x(1,:), LineWidth=1.5)
plot(t, x(2,:), LineWidth=0.5)
plot(t, x(3,:), LineWidth=1.5)
plot(t, x(4,:), LineWidth=1.5)
plot(t, x(5,:), LineWidth=1.5)
plot(t, x(6,:), LineWidth=1.5)
legend('60 Hz', '80 Hz', '100 Hz', '120 Hz','140 Hz','160 Hz')
