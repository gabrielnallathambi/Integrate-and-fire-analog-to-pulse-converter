clear;
%% Input Signal
Fs = 10000; % sampling frequency
signalTime = 0:1/Fs:1; % input time
signalAmplitude =1*sin(2*pi*signalTime*12);
%% IF parameters
pThreshold=0.005;
nThreshold=-pThreshold;
refractoryPeriod =0;
decayRate=0;
timeResolution = 1e-8;
%% Convert analog signals to pulses
outputPulses = apcTimeApprox(signalTime,signalAmplitude,pThreshold,nThreshold,refractoryPeriod,decayRate,timeResolution);
%% Plotting
h=zeros(2,1);
h(1)=subplot (2,1,1);
plot(signalTime,signalAmplitude)
title('Input signal')
h(2)=subplot (2,1,2);
stem(outputPulses(:,1),outputPulses(:,2))
title('Output pulse train');
linkaxes(h,'x')

