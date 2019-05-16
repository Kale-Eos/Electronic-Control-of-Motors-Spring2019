%% EE450 - Pulse Width Modulation (PWM) Simulation
% MATLAB Code developed by Paco Ellaga
% SID: 009602331

clc, clear all, close all;

% Base Variables (temp values, change if needed)
fs = 1000;          % sampling frequency
dt = 1/fs;          % differentional of time
t = 0:dt:1-dt;      % time interval and spacing

Vdc = 100;          % Source Voltage
RLoad = 10;         % Resistor Load
LLoad = 20e-3;      % Inductor Load

% Signal Variables (temp values, change if needed)
fTri = 1260;        % Triangular frequency (amplitude)
fSin = 60;          % Sinusoidal frequency (amplitude)
aTri = 10;          % amplitude of Triangular
aSin = 8;           % amplitude of Sinusoidal -- less than value of aTri
thetaTri = 0;       % Triangular phase delay
thetaSin = 0;       % Triangular phase delay

% Signal ratios 
mFreq = fTri/fSin;  % frequency modulation ratio
mAmp = aSin/aTri;   % amplitude modulation ratio

% Signal Outputs
Vout1 = mAmp*Vdc;   % Output voltage of the fundamental frequency of 1st harmonic

I_1 = Vout1/sqrt((RLoad^2)+(1*2*pi*fSin*LLoad)^2);                  % Current amplitude of 1st harmonic
I_mFreq = Vout1/sqrt((RLoad^2)+(mFreq*2*pi*fSin*LLoad)^2);          % Current amplitude of mFreq value
I_mFreqP2 = Vout1/sqrt((RLoad^2)+((mFreq+2)*2*pi*fSin*LLoad)^2);    % Current amplitude of mFreq value + 2
I_mFreqM2 = Vout1/sqrt((RLoad^2)+((mFreq-2)*2*pi*fSin*LLoad)^2);    % Current amplitude of mFreq value - 2

I_1RMS = I_1/sqrt(2);                       % Current RMS of 1st harmonic
Pabs_1 = ((I_1RMS)^2)*RLoad;                % Power absorbed of 1st harmonic
I_mFreqRMS = I_mFreq/sqrt(2);               % Current RMS of mFreq harmonic
Pabs_mFreq = ((I_mFreqRMS)^2)*RLoad;        % Power absorbed of mFreq harmonic
I_mFreqP2RMS = I_mFreqP2/sqrt(2);           % Current RMS of (mFreq+2) harmonic
Pabs_mFreqP2 = ((I_mFreqP2RMS)^2)*RLoad;    % Power absorbed of (mFreq+2) harmonic
I_mFreqM2RMS = I_mFreqM2/sqrt(2);           % Current RMS of (mFreq+2) harmonic
Pabs_mFreqM2 = ((I_mFreqM2RMS)^2)*RLoad;    % Power absorbed of (mFreq-2) harmonic

Pabs = Pabs_1 + Pabs_mFreq + Pabs_mFreqP2 + Pabs_mFreqM2;   % Total Power absorbed

THD_1 = 100*((sqrt((I_mFreqRMS^2)+(I_mFreqP2RMS^2)+(I_mFreqM2RMS^2)))/I_1RMS);  % Total harmonic distortion of the load current

% Actual Signals
Tri = aTri.*sawtooth(2*pi*fTri*t+thetaTri);     % Triangular signal
Sin = aSin.*sin(2*pi*fSin*t+thetaSin);          % Sinusoidal signal
L = length(Tri);

% need to eval following:
% THD - total harmonic distortion
% Power - reflective of load

% Simulation loop
for i = 1:L
    if (Sin(i) > Tri(i))
        PWM(i) = -1;     % if M value is greater than C return pulse
    else
        PWM(i) = 1;     % else return positive PWM signal
    end
end

figure('Name','PWM Signal Generation and Output Representation'),figure(1)

% Representation of the Sinusoidal Signal
subplot(3,1,1)
hold on, plot(t,Sin,'black','LineWidth',2)
plot(t,Tri), hold off
xlabel('Time (Snapshot)')
xlim([0 0.1])   % may be adjusted based on various input scenarios
ylabel('Amplitude'), % ylim([min(Sin) max(Sin)])
title('Message Signal'), legend('Message Signal','Triangular Signal')
grid minor, box on

% Representation of Triangular Signal
subplot(3,1,2)
hold on, plot(t,Tri)
plot(t,PWM,'red','LineWidth',2), hold off
xlabel('Sample (Snapshot)')
xlim([0 0.1])   % may be adjusted based on various input scenarios
ylabel('Amplitude'), ylim([min(Tri) max(Tri)])
title('PWM Signal Output'), legend('Triangular Signal','PWM Signal')
grid minor, box on

% Power Spectral Density
Welch = pwelch(PWM);    % allocation of Welch's PSD estimate
subplot(3,1,3)
plot(Welch)             % Spectrum of PWelch returned values
xlabel('Sample'), ylabel('Amplitude')
title('PWM Power Spectral Density'), legend('Power Spectral Density')
grid minor, box on

PSD_Min = min(Welch);   % Power Spectral Density Absolute Min 
PSD_Max = max(Welch);   % Power Spectral Density Absolute Max

% Spectrogram representations
figure('Name','Spectragram Evaluations'),figure(2)
subplot(3,1,1)
spectrogram(t,Sin)                      % time versus Sinusoid
title('Sinusoid Signal Spectrogram')
subplot(3,1,2)
spectrogram(t,Tri)                      % time versus Triangular
title('Triangular Signal Spectrogram')
subplot(3,1,3)
spectrogram(t,PWM)                      % time versus PWM
title('PWM SPectrogram')