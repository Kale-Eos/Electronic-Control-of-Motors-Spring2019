%% EE450 - Pulse Width Modulation (PWM) Simulation
% MATLAB Code developed by Paco Ellaga
% SID: 009602331

clc, clear all, close all;

% Variables
t = 0:0.001:1;      % time interval
fTri = 120;         % Triangular frequency
fSin = 100;         % Sinusoidal frequency
a = 10;             % amplitude of Triangular
b = 5;              % amplitude of Sinusoidal -- less than value of a

Tri = a.*sawtooth(2*pi*fTri*t);    % Triangular signal
Sin = b.*sin(2*pi*fSin*t);         % Sinusoidal signal
L = length(Tri);

% Simulation loop
for i = 1:L
    if (Sin(i) >= Tri(i))
        PWM(i) = 1;     % if M value is greater than C return pulse
    else
        PWM(i) = 0;     % else return null signal
    end
end

figure(1)
% Representation of the Message Signal
% subplot(4,1,1)
subplot(3,1,1)
plot(t,Sin,'black')
xlabel('Time'), ylabel('Amplitude')
title('Message Signal'), legend('Message Signal')
grid minor

% Representation of the Carrier Signal
% subplot(4,1,2)
subplot(3,1,2)
hold on
plot(t,Tri)
% xlabel('Sample'), ylabel('Amplitude')
% title('Carrier Signal'), legend('Carrier Signal')
% grid minor

% Representation of the PWM Signal
% subplot(4,1,3);
plot(t,PWM,'red')
hold off
xlabel('Sample'), ylabel('Amplitude')
xlim([0 0.5]), ylim([min(Tri+2) max(Tri-2)])
title('PWM Signal'), legend('PWM Signal')
% axis([0 1 -1 2])
grid minor

% Power Spectral Density
Welch = pwelch(PWM);
% subplot(4,1,4)
subplot(3,1,3)
plot(Welch)
xlabel('Sample'), ylabel('Amplitude')
title('PWM Power Spectral Density'), legend('Power Spectral Density')
grid minor

figure(2)
subplot(3,1,1)
spectrogram(t,PWM)
subplot(3,1,2)
spectrogram(t,Sin)
subplot(3,1,3)
spectrogram(t,Tri)