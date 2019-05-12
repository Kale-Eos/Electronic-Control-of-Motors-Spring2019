%% Lab 2: Example 3-9
%
clear all, close all, clc;

Vrms = 120;         % volts
f = 60;             % hertz
r = 500;            % ohms
c = 100*10^-6;      % Farads

% Part A: Find the voltage output, Vout
Vm = Vrms*sqrt(2);              % volts
W = 377;                        % omega
Wrc = W*r*c;                    % radians
theta = -atan(Wrc)+pi;          % radians
Vtheta = Vm*sin(theta);         % volts

Fs = 1000;                      % sample frequency
dt = 1/Fs;                      % seconds per sample
t = 0:dt:2*pi+theta;            % time duration
t1 = 0:dt:2*pi+theta;           % time duration (duplicate for fig(2))

Vout1 = @(t) Vm*sin(t);
Vout2 = @(t) Vtheta*exp(-(t-theta)/Wrc);
Vout1_1 = Vm*sin(t1);
Vout2_1 = Vtheta*exp(-(t1-theta)/Wrc);

% with C = 100uF
figure(1)
subplot(2,2,1)
syms t
Vout = piecewise(0<=t<=theta, Vout1, theta<=t<=(2*pi+theta), Vout2);
fplot(Vout)
grid minor
title('Measurement of Voltage Output where C = 100uF')
legend('Voltage Output where C = 100uF')
xlabel('Frequency'), xlim([0,8])
ylabel('Amplitude'), ylim([-200,200])

% with C = 100uF
subplot(2,2,2)
hold on
plot(t1,Vout1_1)
plot(t1,Vout2_1)
hold off
grid minor
title('Complete measurement of Voltage Output where C = 100uF')
legend('Vout1','Vout2')
xlabel('Frequency'), ylabel('Amplitude')

c1 = 200*10^-6;                 % Farads
Wrc1 = W*r*c1;                  % radians
theta1 = -atan(Wrc1)+pi;        % radians
Vtheta1 = Vm*sin(theta1);       % volts
t3 = 0:dt:2*pi+theta1;          % time duration (duplicate for fig(4))
Vout3 = @(t) Vm*sin(t);
Vout4 = @(t) Vtheta1*exp(-(t-theta1)/Wrc1);
Vout3_1 = Vm*sin(t3);
Vout4_1 = Vtheta1*exp(-(t3-theta)/Wrc1);

% with C = 200uF
subplot(2,2,3)
syms t
Vout0 = piecewise(0<=t<=theta1, Vout3, theta1<=t<=(2*pi+theta1), Vout4);
fplot(Vout0)
grid minor
title('Measurement of Voltage Output where C = 200uF')
legend('Voltage Output where C = 200uF')
xlabel('Frequency'), xlim([0,8])
ylabel('Amplitude'), ylim([-200,200])

% with C = 200uF
subplot(2,2,4)
hold on
plot(t3,Vout3_1)
plot(t3,Vout4_1)
hold off
grid minor
title('Complete measurement of Voltage Output where C = 200uF')
legend('Vout3','Vout4')
xlabel('Frequency'), ylabel('Amplitude')

% Part 2: Find peak-to-peak ripple, Delta Vout
alpha_fun = @(t) sin(t)-sin(theta)*exp((-2*pi+t-theta)/Wrc);
alpha = fzero(alpha_fun,0.1);

deltaVout = Vm*(1-sin(alpha));      % standard answer using equation 3-49
deltaVout1 = Vm*((2*pi)/(Wrc));     % approximated with equation 3-51

alpha_funC = @(t) sin(t)-sin(theta1)*exp((-2*pi+t-theta1)/Wrc1);
alphaC = fzero(alpha_funC,0.1);

% with C = 200uF
deltaVout2 = Vm*(1-sin(alphaC));     % standard answer using equation 3-49
deltaVout3 = Vm*((2*pi)/(Wrc1));     % approximated with equation 3-51

% Part 3: Find the capacitor current, Ic
Ic1 = @(t) (-Vtheta/r)*exp(-(t-theta)/Wrc);
Ic2 = @(t) W*c*Vm*cos(t);
Ic1_1 = (-Vtheta/r)*exp(-(t1-theta)/Wrc);
Ic2_1 = W*c*Vm*cos(t1);

% with C = 100uF
figure(2)
subplot(2,2,1)
syms t
IcOut = piecewise(0<=t<=theta, Ic1, theta<=t<=(2*pi+theta), Ic2);
fplot(IcOut)
grid minor
title('Measurement of Capacitor Current where C = 100uF')
legend('Capacitor Current where C = 100uF')
xlabel('Frequency'), xlim([0,8])
ylabel('Amplitude'), ylim([-16,16])

% with C = 100uF
subplot(2,2,2)
hold on
plot(t1,Ic1_1)
plot(t1,Ic2_1)
hold off
grid minor
title('Complete measurement of Capacitor Current where C = 100uF')
legend('Ic1','Ic2')
xlabel('Frequency')
ylabel('Amplitude'), ylim([-16,16])

Ic3 = @(t) (-Vtheta1/r)*exp(-(t-theta1)/Wrc1);
Ic4 = @(t) W*c1*Vm*cos(t);
Ic3_1 = (-Vtheta1/r)*exp(-(t3-theta1)/Wrc1);
Ic4_1 = W*c1*Vm*cos(t3);

% with C = 200uF
subplot(2,2,3)
syms t
IcOut0 = piecewise(0<=t<=theta1, Ic3, theta<=t<=(2*pi+theta1), Ic4);
fplot(IcOut0)
grid minor
title('Measurement of Capacitor Current where C = 200uF')
legend('Capacitor Current where C = 200uF')
xlabel('Frequency'), xlim([0,8])
ylabel('Amplitude'), ylim([-16,16])

% with C = 200uF
subplot(2,2,4)
hold on
plot(t3,Ic3_1)
plot(t3,Ic4_1)
hold off
grid minor
title('Complete measurement of Capacitor Current where C = 200uF')
legend('Ic3','Ic4')
xlabel('Frequency')
ylabel('Amplitude'), ylim([-16,16])

% Part 4: Find the peak diode current, Id
% Uses equation 3-48
Id = Vm*(W*c*cos(alpha)+(sin(alpha)/r));        % peak diode current in amps

Id1 = Vm*(W*c1*cos(alphaC)+(sin(alphaC)/r));    % peak diode current in amps

% Part 5: Find the value of c such that Delta Vout is 1% of Vm
% deltaVout = 0.01*Vm
% Uses equation 3-51
C = Vm/(f*r*(0.01*Vm));    % Farads