%% EE450 - Pulse Width Modulation (PWM) Simulation
% MATLAB Code developed by Paco Ellaga
% SID: 009602331

clc, clear all, close all;

% Variables
t = 0:0.001:1;      % time interval and spacing
fTri = 120;         % Triangular frequency
fSin = 120;         % Sinusoidal frequency
a = 10;             % amplitude of Triangular
b = 10;             % amplitude of Sinusoidal -- less than value of a
thetaTri = -pi/2;   % Triangular phase delay
thetaSin = -pi/2;   % Triangular phase delay

Tri = a.*sawtooth(2*pi*fTri*t+thetaTri);    % Triangular signal
Sin = b.*sin(2*pi*fSin*t+thetaSin);         % Sinusoidal signal
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
% Representation of the Sinusoidal Signal
subplot(3,1,1)
hold on
plot(t,Sin,'black','LineWidth',2)
plot(t,Tri)
hold off
xlabel('Time'), xlim([0 0.1])
ylabel('Amplitude'), ylim([min(Sin) max(Sin)])
title('Message Signal'), legend('Message Signal','Triangular Signal')
grid minor

% Representation of Triangular Signal
subplot(3,1,2)
hold on
plot(t,Tri)
plot(t,PWM,'red','LineWidth',2)
hold off
xlabel('Sample'), xlim([0 0.1])
ylabel('Amplitude'), ylim([min(Tri) max(Tri)])
title('PWM Signal Output'), legend('Triangular Signal','PWM Signal')
grid minor

% Power Spectral Density
Welch = pwelch(PWM);
subplot(3,1,3)
plot(Welch)
xlabel('Sample'), ylabel('Amplitude')
title('PWM Power Spectral Density'), legend('Power Spectral Density')
grid minor

figure(2)
% Spectrogram representations
subplot(3,1,1)
spectrogram(t,PWM)
title('PWM SPectrogram')
subplot(3,1,2)
spectrogram(t,Sin)
title('Sinusoid Signal Spectrogram')
subplot(3,1,3)
spectrogram(t,Tri)
title('Triangular Signal Spectrogram')