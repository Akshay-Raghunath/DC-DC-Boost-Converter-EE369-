clc;
clear all;


vm=230*sqrt(2);
f=50;
wc=2*pi*20;

%boost converter parameters
Vin = 24;
Vo = 48;
D = 1-(Vin/Vo);
R = 28;
Po = (Vo^2)/R;
fsw = 20e03;
L = 500e-6;
rL = 0.4;
C = 680e-06;
Iin = Po/Vin;

%boost converter modelling 

A1 = [-rL/L 0;0 -1/(C*R)];
B1 = [1/L;0];
C1 = [0 1];

A2 = [-rL/L -1/L;1/C -1/(C*R)];
B2 = [1/L;0];
C2 = [0 1];

A = A1*D + A2*(1-D);
B = B1*D + B2*(1-D);
%C = C1*D + C2*(1-D);

X = -(A^(-1))*B*Vin; %Calculates the steady-state operating point of the converter (inductor current and output voltage)
X_ver = [Iin;Vo];

s = tf('s');
sI = [s 0 ;0 s];
x_by_d = ((sI-A)^-1)*((A1-A2)*X + (B1-B2)*Vin); %Calculates the transfer function of the inductor current with respect to duty cycle (i_by_d)
i_by_d = x_by_d(1);
%v0_by_d = (C*((sI-A)^-1)*((A1-A2)*X + (B1-B2)*Vin)) + (C1-C2)*X;

%LPF
wc_lpf=2*pi*30;
lpf_tf=wc_lpf/(s+wc_lpf);
tc_lpfd=c2d(lpf_tf,10e-6,'tustin');


%k-factor method

wc=2*pi*1000;
DPM=55;

s=tf('s');
plant=i_by_d *lpf_tf;

[mag,phi]=bode(plant,wc)

boost=DPM-90-phi;
k=tand((boost/2)+45);
wz=wc/k;
wp=wc*k;
A=(k*wc)/mag;

tc1=A/s;
tc2=(s+wz)/(s+wp);
tc=tc1*tc2;
tcd1=c2d(tc1,10e-6,'tustin');
tcd2=c2d(tc2,10e-6,'tustin');

margin(plant*tc)
