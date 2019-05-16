%% EE450 - Pulse Width Modulation (PWM) Simulation
% MATLAB Code developed by Paco Ellaga
% SID: 009602331

clc, clear all, close all;

% Base Variables (temp values, change if needed)
fs = 1000;          % sampling frequency
dt = 1/fs;          % differential of time
t = 0:dt:1-dt;      % time interval and spacing

Vdc = 100;          % Source Voltage
RLoad = 10;         % Resistor Load
LLoad = 20e-3;      % Inductor Load

fSin = 60;          % Sinusoidal frequency (1st amplitude harmonic)
aSin = 8;           % amplitude of Sinusoidal (less than or equal to value of aTri)
aTri = 10;          % amplitude of Triangular
thetaTri = 0;       % Triangular phase delay in terms of pi
thetaSin = 0;       % Sinusoid phase delay in terms of pi

% Signal modulation ratios (temp values, change if needed) 
mFreq = 21;             % fTri/fSin;  % frequency modulation ratio and target harmonic
mFreqP2 = mFreq + 2;    % harmonic + 2
mFreqM2 = mFreq - 2;    % harmonic - 2

mAmp = aSin/aTri;       % amplitude modulation ratio
% Following 2 variables can be plugged by matrix
mAmp_mFreq = 0.82;      % Note: Just plugged these from table 8-3
mAmp_mFreq2 = 0.22;     % I can't do this anymore, will come back later

fTri = mFreq*fSin;      % Triangular modulation frequency (amplitude harmonic)
fTriP2 = mFreqP2*fSin;  % harmonic + 2
fTriM2 = mFreqM2*fSin;  % harmonic - 2

% Signal Outputs
Vout_1 = mAmp*Vdc;              % Output voltage of the fundamental frequency of 1st harmonic
Vout_mFreq = mAmp_mFreq*Vdc;    % From target harmonic
Vout_mFreq2 = mAmp_mFreq2*Vdc;  % From target harmonic +- 2

I_1 = Vout_1/sqrt((RLoad^2)+(1*2*pi*fSin*LLoad)^2);                     % Current amplitude from 1st harmonic
I_mFreq = Vout_mFreq/sqrt((RLoad^2)+(mFreq*2*pi*fSin*LLoad)^2);         % From mFreq value
I_mFreqP2 = Vout_mFreq2/sqrt((RLoad^2)+((mFreqP2)*2*pi*fSin*LLoad)^2);  % From mFreq value + 2
I_mFreqM2 = Vout_mFreq2/sqrt((RLoad^2)+((mFreqM2)*2*pi*fSin*LLoad)^2);  % From mFreq value - 2

I_1RMS = I_1/sqrt(2);                   % Current RMS from 1st harmonic
I_mFreqRMS = I_mFreq/sqrt(2);           % From mFreq harmonic
I_mFreqP2RMS = I_mFreqP2/sqrt(2);       % From mFreq harmonic + 2
I_mFreqM2RMS = I_mFreqM2/sqrt(2);       % From mFreq harmonic -2

P_1 = ((I_1RMS)^2)*RLoad;               % Power absorbed from 1st harmonic
P_mFreq = ((I_mFreqRMS)^2)*RLoad;       % From mFreq harmonic
P_mFreqP2 = ((I_mFreqP2RMS)^2)*RLoad;   % From mFreq harmonic + 2
P_mFreqM2 = ((I_mFreqM2RMS)^2)*RLoad;   % From mFreq harmonic - 2

Pabs = P_1 + P_mFreq + P_mFreqP2 + P_mFreqM2;   % Total Power absorbed

THD_1 = 100*((sqrt((I_mFreqRMS^2)+(I_mFreqP2RMS^2)+(I_mFreqM2RMS^2)))/I_1RMS);  % Total harmonic distortion of the load current

% Actual Graphed Signals
Tri = aTri.*sawtooth(2*pi*fTri*t+thetaTri);     % Triangular signal
Sin = aSin.*sin(2*pi*fSin*t+thetaSin);          % Sinusoidal signal
L = length(Tri);

% Simulation loop
for i = 1:L
    if (Sin(i) > Tri(i))
        PWM(i) = -1;    % if M value is greater than C return pulse
    else
        PWM(i) = 1;     % else return positive PWM signal
    end
end

figure('Name','PWM Signal Generation and Output Representation'),figure(1)

% Representation of the Sinusoidal Signal
subplot(2,1,1)
hold on, plot(t,Sin,'black','LineWidth',2)
plot(t,Tri), hold off
xlabel('Time (Snapshot)')
xlim([0 0.2])   % may be adjusted based on various input scenarios
ylabel('Amplitude'), % ylim([min(Sin) max(Sin)])
title('Message Signal'), legend('Message Signal','Triangular Signal')
grid minor, box on

% Representation of Triangular Signal
subplot(2,1,2)
hold on, plot(t,Tri)
plot(t,PWM,'red','LineWidth',2), hold off
xlabel('Sample (Snapshot)')
xlim([0 0.2])   % may be adjusted based on various input scenarios
ylabel('Amplitude'), ylim([min(Tri) max(Tri)])
title('PWM Signal Output'), legend('Triangular Signal','PWM Signal')
grid minor, box on

figure('Name','Spectragram Evaluations'),figure(2)
% Spectrogram representations
subplot(3,1,1)
spectrogram(t,Sin)                      % time versus Sinusoid
title('Sinusoid Signal Spectrogram')
subplot(3,1,2)
spectrogram(t,Tri)                      % time versus Triangular
title('Triangular Signal Spectrogram')
subplot(3,1,3)
spectrogram(t,PWM)                      % time versus PWM
title('PWM Spectrogram')

figure('Name','PSD Versus PWM'),figure(3)
% Power Spectral Density (PWelsh)
Welch = pwelch(PWM);    % allocation of Welch's PSD estimate
subplot(2,1,1)
plot(Welch)             % Spectrum of PWelch returned values
xlabel('Sample'), ylabel('Amplitude')
title('Power Spectral Density of PWM'), legend('Power Spectral Density')
grid minor, box on

PSD_Min = min(Welch);   % Power Spectral Density Absolute Min 
PSD_Max = max(Welch);   % Power Spectral Density Absolute Max

subplot(2,1,2)
spectrogram(t,PWM)                      % time versus PWM
title('PWM Spectrogram (Reshown for comparison)')