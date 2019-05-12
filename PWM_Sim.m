%% EE450 - Pulse Width Modulation (PWM) Simulation
% MATLAB Code developed by Paco Ellaga
% SID: 009602331

clc, clear all, close all;

% Variables
t = 0:0.001:1;  % time interval
fc = 1200;      % carrier frequency
fm = 1000;      % modulated frequency
a = 10;         % amplitude of carrier
b = 5;          % amplitude of modulated -- less than value of a
C = a.*sawtooth(2*pi*fc*t);    % carrier signal
M = b.*sin(2*pi*fm*t);         % modulated signal
L = length(C);

% Simulation loop
for i = 1:L
    if (M(i)>=C(i))
        PWM(i) = 1;     % if M value is greater than C return pulse
    else
        PWM(i) = 0;     % else return null signal
    end
end

figure(1)
% Representation of the Message Signal
subplot(4,1,1)
plot(t,M,'black')
xlabel('Time'), ylabel('Amplitude')
title('Message Signal'), legend('Message Signal')
grid minor

% Representation of the Carrier Signal
subplot(4,1,2)
plot(t,C)
xlabel('Sample'), ylabel('Amplitude')
title('Carrier Signal'), legend('Carrier Signal')
grid minor

% Representation of the PWM Signal
subplot(4,1,3);
plot(t,PWM,'red')
xlabel('Sample'), ylabel('Amplitude')
title('PWM Signal'), legend('PWM Signal')
axis([0 1 -1 2]), grid minor

% Power Spectral Density
Welch = pwelch(PWM);
subplot(4,1,4)
plot(Welch)
xlabel('Sample'), ylabel('Amplitude')
title('PWM Power Spectral Density'), legend('Power Spectral Density')
grid minor

figure(2)
subplot(3,1,1)
spectrogram(t,PWM)
subplot(3,1,2)
spectrogram(t,M)
subplot(3,1,3)
spectrogram(t,C)