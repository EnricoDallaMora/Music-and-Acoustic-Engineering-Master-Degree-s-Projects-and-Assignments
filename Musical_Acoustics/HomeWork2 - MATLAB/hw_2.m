%DATA
a=0.15; %side length [m]
th=0.001; %thickness [m]
E=69*10^9; % Young Modulus [Pa]
rho=2700; %density [Kg/m^3]
ni=0.334; %poisson's ratio

c_ql=sqrt(E/(rho*(1-ni^2))); %quasi-longitudinal waves propagation speed [m/s]
c_l=sqrt(E*(1-ni)/(rho*(1+ni)*(1-2*ni))); %longitudinal waves propagation speed [m/s]

%define omega axis
f = 0:1:10000;
c_b=sqrt(2*pi.*f.*th*c_ql/sqrt(12)); %bending waves propagation speed [m/s]
figure(1);
plot(f, c_b,'r','LineWidth',1.2);
xlabel('Frequency [Hz]','Fontsize',12); 
ylabel('Speed [m/s]','Fontsize',12);
title('Propagation speed of bending waves','Fontsize',18);


%modal freqs    
ratios=[1; 2.04; 2.04; 3.01; 3.66; 3.67];
f_00=1.654*c_ql*th/a^2; %freq of mode (0,0) [Hz]
f_string=f_00;
f_ii=f_00.*ratios; %frequencies for the ij modes [Hz]


%Stikaz spruce
E_l=12; %Young modulus logitudinally
E_r=0.9; %Young modulus radially
length_ratio=(E_l/E_r)^(1/4); 
b=length_ratio*a; %length of the b side [m]



%STRING DATA
rho_string=5000; %density of the string [kg/m^3]
L=0.45; %length of the string [m]
radius=0.0011; %radius of the string [m]

mu_string=rho_string*pi*radius^2; %linear density of the string [kg/m]
%f=c/2L;
%c=(T/mu)^1/2 

T=(f_string*2*L)^2*mu_string;


%COUPLING CONDITION
Q=25;
m_string=rho_string*L*pi*radius^2; %string mass [kg]
m_plate=rho*a^2*th;

if m_string/(1*m_plate)>(pi^2)/(4*Q^2)
    coupling="strong";
else
    coupling="weak";
end

ratio11=(f_string*1-f_ii(1))/f_ii(1);   %0.02
ratio22=(f_string*2-f_ii(2))/f_ii(2);   %+0.032/-0.013
ratio23=(f_string*2-f_ii(3))/f_ii(3);   %+0.032/-0.013
ratio34=(f_string*3-f_ii(4))/f_ii(4);   %0.02
ratio45=(f_string*4-f_ii(5))/f_ii(5);   %fuori dal grafico
ratio46=(f_string*4-f_ii(6))/f_ii(6);   %fuori dal grafico

f_system11_1=f_ii(1)*0.98;
f_system11_2=f_ii(1)*1.02;

f_system22_1=f_ii(2)*0.968;
f_system22_2=f_ii(2)*1.013;

f_system23_1=f_ii(3)*0.968;
f_system23_2=f_ii(3)*1.013;

f_system34_1=f_ii(4)*0.98;
f_system34_2=f_ii(4)*1.02;

f_system45_1=f_string*4;
f_system45_2=f_ii(5);

f_system46_1=f_string*4;
f_system46_2=f_ii(6);
