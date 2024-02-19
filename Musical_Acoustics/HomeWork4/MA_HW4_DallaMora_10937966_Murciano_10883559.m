clear
close all
clc
L=0.45;                %m
alpha=deg2rad(0.75);   %rad
c=343;                 %m/s
delta_L_mouth=0.04;          %m
rho=1.225;             %kg/m^3
%% Primo punto

f_0=329.63;            %Hz
w_0=2*pi*f_0;          %rad/s
k_0=w_0/c;             %wavenumber

r_in = linspace(0, 0.06, 1000); %independent variable
r_out = r_in - L*tan(alpha);

S_in = r_in.^2*pi;
S_out = r_out.^2*pi;

M = (delta_L_mouth*rho)./S_in;
Z_mouth = 1i*w_0*M;

L1 = L + (0.85*r_out); %Lp
x_out = r_out./tan(alpha);
theta_1 = atan(k_0*x_out)./k_0;
Z_cone = ((1i*rho*c)./(S_in)).*((sin(k_0.*L1).*sin(k_0.*theta_1))./(sin(k_0.*(L1+theta_1))));

Z_tot=Z_mouth+Z_cone;
Z_tot_db = db(Z_tot);

figure(1)
plot(r_in,Z_tot_db,LineWidth=1.2)
xlabel("$R_{in}\ [m]$", Interpreter='latex', FontSize=20); 
ylabel("$|Z_{tot}|\ [dB]$", Interpreter='latex', FontSize=20)

hold on;
r_in = r_in(Z_tot_db == min(Z_tot_db));
xline(r_in, 'k--', LineWidth=1.4, Color='r')
text(r_in*1.01, min(Z_tot_db)/1.1, ["$R_{in}=$"+num2str(r_in)+" m"], Interpreter="latex", FontSize=14)
r_out = r_in - L*tan(alpha);
%% Secondo punto

d_h = linspace(0, 0.45, 1000);

f1 = 349.23; 
w1 = 2*pi*f1;
k_1 = w1/c;

x_out = r_out/tan(alpha);
S_o = r_out^2*pi;
S_i = r_in^2*pi;

deltaL=0.85*r_out;
L1 = L + deltaL; 
delta_red = d_h + (deltaL.^2./(d_h+2*deltaL));
L11 = L1-delta_red; 

M = delta_L_mouth*rho./S_i;
x_out_1 = delta_red+x_out-deltaL;

theta_1 = atan(k_1*x_out_1)./k_1;

Z_in = 1i*w1*M + ((1i*rho*c)./(S_i)).*((sin(k_1.*L11).*sin(k_1.*theta_1))./(sin(k_1.*(L11+theta_1))));
Zindb = db(Z_in);
    
figure(2)
plot(d_h, Zindb, LineWidth=1.4)
hold on
xlabel("$d_h\ [m]$", Interpreter='latex', FontSize=20); 
ylabel("$|Z_{in}|\ [dB]$", Interpreter='latex', FontSize=20)

d_h = d_h(Zindb == min(Zindb));
xline(d_h, 'k--', LineWidth=1.4, Color='r')
text(d_h*1.1, min(Zindb)/1.1, ["$d_h=$"+num2str(d_h)+" m"], Interpreter="latex", FontSize=14)
hole1_coord = L-d_h;


%% Terzo punto

delta_red_h1 = d_h + (deltaL^2/(d_h+2*deltaL)); %acoustic reduction delta con l'ultimo buco aperto
L11_h1 = L1 - delta_red_h1;                
f2 = 392;
w2 = 2*pi*f2;
k2 = w2/c;

d_h2 = linspace(0, L11_h1, 1000);
delta2 = d_h2 - (d_h2*deltaL)./(d_h2+deltaL);
L11_h2 = L11_h1 - delta2;
theta2 = atan(k2*(x_out+delta_red_h1+delta2-deltaL))./k2;
Z_in2 = 1i*w2*M + ((1i*rho*c)./(S_i)).*((sin(k2.*L11_h2).*sin(k2.*theta2))./(sin(k2.*(L11_h2+theta2))));
Zindb2 = db(Z_in2);
figure(3)
plot(d_h2, Zindb2, LineWidth=1.4)
hold on
xlabel("$d_{h2}\ [m]$", Interpreter='latex', FontSize=20); 
ylabel("$|Z_{in}|\ [dB]$", Interpreter='latex', FontSize=20)
xlim([0 L11_h1])

d_h2= d_h2(Zindb2 == min(Zindb2));
xline(d_h2, 'k--', LineWidth=1.4,Color='r')
text(d_h2*1.1, min(Zindb2)*1.3, ["$d_{h2}=$"+num2str(d_h2)+" m"], Interpreter="latex", FontSize=14)
hole2_coord = L11_h1-d_h2;


%% Quarto punto 

b = r_out;


S_in = pi*r_in^2;
S_out = pi*r_out^2;

f=20:4000;
omega=2*pi*f;
Z_mouth = 1i*omega*M;

l_1_in = -r_in/tan(alpha);
l_1_out = l_1_in + L - d_h2;
l_2_in = l_1_out;
l_2_out = l_2_in + d_h2 - d_h;
l_3_in = l_2_out;
l_3_out = l_3_in + d_h;

r_1_in = r_in;
r_1_out = -l_1_out*tan(alpha);
r_2_in = r_1_out;
r_2_out = -l_2_out*tan(alpha);
r_3_in = r_2_out;
r_3_out = r_out;

S_1_in = pi*r_1_in^2;
S_2_in = pi*r_2_in^2;
S_3_in = pi*r_3_in^2;

S_1_out = pi*r_1_out^2;
S_2_out = pi*r_2_out^2;
S_3_out = pi*r_3_out^2;

L1_in = (rho*l_1_in)/S_1_in; ZL_1_in = 1i*omega*L1_in;
L2_in = (rho*l_2_in)/S_2_in; ZL_2_in = 1i*omega*L2_in;
L3_in = (rho*l_3_in)/S_3_in; ZL_3_in = 1i*omega*L3_in;
L1_out = (rho*l_1_out)/S_1_out; ZL_1_out = 1i*omega*L1_out;
L2_out = (rho*l_2_out)/S_2_out; ZL_2_out = 1i*omega*L2_out;
L3_out = (rho*l_3_out)/S_3_out; ZL_3_out = 1i*omega*L3_out;

a_1 = r_1_out;
a_2 = r_2_out;
k = omega/c;
t_1 = deltaL + 0.125*b*(b/a_1)*(1 + 0.172*(b/a_1)^2);
t_2 = deltaL + 0.125*b*(b/a_2)*(1 + 0.172*(b/a_2)^2);
t_e_1 = ((1./k).*tan(k.*t_1) + b*(1.4 - 0.58*(b/a_1)^2))./(1 - 0.61*k*b.*tan(k*t_1));
t_e_2 = ((1./k).*tan(k.*t_2) + b*(1.4 - 0.58*(b/a_2)^2))./(1 - 0.61*k*b.*tan(k*t_2));
t_a_1 = (0.47*b*(b/a_1)^4)/(tanh(1.84*t_1./b) + 0.62*(b/a_1)^2 + 0.64*(b/a_1));
t_a_2 = (0.47*b*(b/a_2)^4)/(tanh(1.84*t_2./b) + 0.62*(b/a_2)^2 + 0.64*(b/a_2));

theta_11 = atan((k)*(-l_3_out+deltaL))./(k);

Z_a1_mezzi = (1/2)*((rho*c)/(pi*b^2))*(-1i*(k).*t_a_1);
Z_a2_mezzi = (1/2)*((rho*c)/(pi*b^2))*(-1i*(k).*t_a_2);
Z_s_cl_1 = ((rho*c)/(pi*b^2))*(-1i*cot((k).*t_1));
Z_s_cl_2 = ((rho*c)/(pi*b^2))*(-1i*cot((k).*t_2));
Z_s_op_1 = ((rho*c)/(pi*b^2))*(-1i*(k).*t_e_1);
Z_s_op_2 = ((rho*c)/(pi*b^2))*(-1i*(k).*t_e_2);
Z_L = 1i*((rho*c)/S_3_out).*(sin(k*deltaL).*sin(k.*theta_11)./sin(k.*(deltaL+theta_11)));
%Z_L = 1i*((rho*c)/S_3_out)*tan((omega./c).*deltaL);

z_0_1 = 1i*((rho*c)/S_1_out)*tan(k.*(-l_1_in+l_1_out + 0.6*r_1_out));
z_0_2 = 1i*((rho*c)/S_2_out)*tan(k.*(-l_2_in+l_2_out + 0.6*r_2_out));
z_0_3 = 1i*((rho*c)/S_3_out)*tan(k.*(-l_3_in+l_3_out + 0.85*r_3_out));

% z_0_1 = 1i*((rho*c)/S_1_out)*tan(k.*(-l_1_in+l_1_out));
% z_0_2 = 1i*((rho*c)/S_2_out)*tan(k.*(-l_2_in+l_2_out));
% z_0_3 = 1i*((rho*c)/S_3_out)*tan(k.*(-l_3_in+l_3_out));


% open
%Z_1_o = Z_L.*ZL_3_out./(Z_L + ZL_3_out);
Z_1_o = 0;
%Z_2_o = (Z_1_o + z_0_3).*ZL_3_in./(Z_1_o + z_0_3 + ZL_3_in);
Z_2_o = 0;
%Z_3_o = (Z_2_o + Z_a2_mezzi).*Z_s_op_2./(Z_2_o + Z_a2_mezzi + Z_s_op_2);
Z_3_o = 0;
%Z_4_o = (Z_3_o + Z_a2_mezzi).*ZL_2_out./(Z_3_o + Z_a2_mezzi + ZL_2_out);
Z_4_o = 0;
%Z_5_o = (Z_4_o + z_0_2).*ZL_2_in./(Z_4_o + z_0_2 + ZL_2_in);
Z_5_o = 0;
Z_6_o = (Z_5_o + Z_a1_mezzi).*Z_s_op_1./(Z_5_o + Z_a1_mezzi + Z_s_op_1);
Z_7_o = (Z_6_o + Z_a1_mezzi).*ZL_1_out./(Z_6_o + Z_a1_mezzi + ZL_1_out);
Z_fin_open = (Z_7_o + z_0_1).*ZL_1_in./(Z_7_o + z_0_1 + ZL_1_in)+ Z_mouth;

% closed
Z_1_c = Z_L.*ZL_3_out./(Z_L + ZL_3_out);
Z_2_c = (Z_1_c + z_0_3).*ZL_3_in./(Z_1_c + z_0_3 + ZL_3_in);
Z_3_c = (Z_2_c + Z_a2_mezzi).*Z_s_cl_2./(Z_2_c + Z_a2_mezzi + Z_s_cl_2);
Z_4_c = (Z_3_c + Z_a2_mezzi).*ZL_2_out./(Z_3_c + Z_a2_mezzi + ZL_2_out);
Z_5_c = (Z_4_c + z_0_2).*ZL_2_in./(Z_4_c + z_0_2 + ZL_2_in);
Z_6_c = (Z_5_c + Z_a1_mezzi).*Z_s_cl_1./(Z_5_c + Z_a1_mezzi + Z_s_cl_1);
Z_7_c = (Z_6_c + Z_a1_mezzi).*ZL_1_out./(Z_6_c + Z_a1_mezzi + ZL_1_out);
Z_fin_closed = (Z_7_c + z_0_1).*ZL_1_in./(Z_7_c + z_0_1 + ZL_1_in)+ Z_mouth;

% open-closed 
Z_1_oc = 0;
Z_2_oc = 0;
Z_3_oc = (Z_2_oc + Z_a2_mezzi).*Z_s_op_2./(Z_2_oc + Z_a2_mezzi + Z_s_op_2);
Z_4_oc = (Z_3_oc + Z_a2_mezzi).*ZL_2_out./(Z_3_oc + Z_a2_mezzi + ZL_2_out);
Z_5_oc = (Z_4_oc + z_0_2).*ZL_2_in./(Z_4_oc + z_0_2 + ZL_2_in);
Z_6_oc = (Z_5_oc + Z_a1_mezzi).*Z_s_cl_1./(Z_5_oc + Z_a1_mezzi + Z_s_cl_1);
Z_7_oc = (Z_6_oc + Z_a1_mezzi).*ZL_1_out./(Z_6_oc + Z_a1_mezzi + ZL_1_out);
Z_fin_open_closed = (Z_7_oc + z_0_1).*ZL_1_in./(Z_7_oc + z_0_1 + ZL_1_in) + Z_mouth;

figure(5)
subplot(3, 1, 1);
plot(f, db(abs(Z_fin_closed)), LineWidth=1.4);
title("Impedance with two closed holes");
xlabel("$Frequency \ [Hz]$", Interpreter='latex', FontSize=12); 
ylabel("$|Z_{in}| \ [dB]$", Interpreter='latex', FontSize=12);

subplot(3, 1, 2);
plot(f, db(abs(Z_fin_open_closed)), LineWidth=1.4);
title("Impedance with only the last hole open");
xlabel("$Frequency \ [Hz]$", Interpreter='latex', FontSize=12); 
ylabel("$|Z_{in}| \ [dB]$", Interpreter='latex', FontSize=12);

subplot(3, 1, 3);
plot(f, db(abs(Z_fin_open)), LineWidth=1.4);
title("Impedance with two open holes");
xlabel("$Frequency \ [Hz]$", Interpreter='latex', FontSize=12); 
ylabel("$|Z_{in}| \ [dB]$", Interpreter='latex', FontSize=12);



