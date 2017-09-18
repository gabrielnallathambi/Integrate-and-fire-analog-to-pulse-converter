function outputPulses = apcTimeApprox(signalTime,signalAmplitude,pThreshold,nThreshold,refractoryPeriod,decayRate,timeResolution)
%APCTIMEAPPROX Apc converter.
% apcTimeApprox.m : Implements biphasic Integrate and Fire Analog to Pulse
% converter by approximating the equation: Integrate from t1 to t2
% [x(t)exp(-a(t2-t)).dt], where 'a' is rate of decay.
%
% Inputs:
% Signal description:
%     signalTime - time of the input signal.
%     signalAmplitude - amplitude of the input signal
% Analog to pulse converter parameters:
%     pThreshold - positive threshold 
%     nThreshold - negative threshold 
%     refractoryPeriod - refractory period (accepts zero or non-zero scalars)
%     decayRate - rate of decay (accepts zero or non-zero scalars)
%     timeResolution - time resolution of integration
%    
% Outputs:
% outputPulses - provides time of occurence of pulses (column 1) and its corresponding threshold polarity (column 2)
%
%   % Example:
%   %   For data sampled at 10 KHz, convert the data into
%   %   pulses with the following parameters.
%
%   Fs = 10000; % sampling frequency
%   signalTime = 0:1/Fs:0.1;            % input time
%   signalAmplitude=2*sin(2*pi*t*12);   % input signal
%   pThreshold=0.001;                   % positive threshold
%   nThreshold=-pThreshold;             % negative threshold
%   refractoryPeriod =0;                % refractory period
%   decayRate=40;                       % decay rate
%   timeResolution = 1e-8;              % time resolution
%   outputPulses = apcTimeApprox(signalTime,signalAmplitude,pThreshold,nThreshold,refractoryPeriod,decayRate,timeResolution);
% Last Updated : 1-18-2017
% References:
%     [1] Alexander singh alvarado, TIME ENCODED COMPRESSION AND CLASSIFICATION USING THE INTEGRATE AND FIRE SAMPLER Jan 2012, University of Florida.

% Acknowledments: Special thanks to Alex for walking me through the
% implementation of his IFC code, which helped me to write this modified
% code.

disp(sprintf(['\n' , '%s', '%d', '%s'], 'Analog to Pulse conversion - Status: Initializing'));
%% 1a. Initialize input and output parameters
x=signalTime; % x axis of input signal (time)
y=signalAmplitude; % y axis of input signal (amplitude)
outputPulses=zeros(1e6,2); % initialize to a big number
integratedArea = 0; % area under the curve
%% 1b. Set indexes and flags
signalTimeIndex=1; % index number for input signal time and amplitude
pulseIndex=0; % index # for the output pulses
zeroCrossing = false; % flag denoting change in polarity of amplitude
outOfBounds = false; % flag denoting "signalTimeIndex" exceeding length of input signal

%% 2. Integration
while(signalTimeIndex < length(signalTime));
    pulseDetected = false; % flag denoting the detection of output pulse
    %% 2a. Determine area co-ordinates {(x1,y1),(x2,y2)} --> boundary points
    x1 = x(signalTimeIndex);
    y1 = y(signalTimeIndex);
    x2 = x(signalTimeIndex+1);
    y2 = y(signalTimeIndex+1);
    %% Check for Zero Crossing and update area co-ordinates (x2,y2)
    if(y1*y2 < 0) % change in polarity
        m= (y2-y1)/(x2-x1); % slope
        c= y2 - m*x2; % intercept
        x2 = -c/m; % time of zero-crossing (0=mx+b)
        y2 = 0;
        zeroCrossing = true;
    end
    %% 2b. Calculate area under the curve (coarse approximation)
    previousIntegratedArea=integratedArea;
    integratedArea = apcIntegration( x1,x2,y1,y2,decayRate,previousIntegratedArea);
    %% 3. Positive Comparator
    if (integratedArea >= pThreshold)
        %% 3a. Calculate area under the curve (fine approximation)
        pulseTime = apcFineIntegration(x1,y1,x2,y2,decayRate,previousIntegratedArea,pThreshold,timeResolution); %fine time at which the threshold is reached
        %% 3b. Capture the pulse time and its polarity
        if (pulseIndex==size(outputPulses,1)) % if size of output pulses exceeds the initial value, then extend its size
            outputPulses=[outputPulses;zeros(1e6,2)];
        end
        pulseIndex=pulseIndex+1;
        outputPulses(pulseIndex,1) = pulseTime;
        outputPulses(pulseIndex,2) = pThreshold;
        %% 3c. Blank for refractory period and reset area
        x1 = pulseTime+refractoryPeriod;
        integratedArea=0;
        %% 3d. Check for the end of the signal (out of bounds)
        if(x1 > signalTime(end))
            outOfBounds = true;
            break;
        end
        %% 3e. Update x-axis and y-axis with new co-ordinates
        remainingTimeVector = find(signalTime>x1);
        if(isempty(remainingTimeVector))
            break;
        end
        % find y-axis point (y1) by interploation
        xint = [];
        yint = [];
        xint = signalTime([max(remainingTimeVector(1)-3,1):min(remainingTimeVector(1)+3,length(signalTime))]); % find time points for interpolation (max:choose nearby point or the first point, min: choose nearby point or last point)
        yint = signalAmplitude([max(remainingTimeVector(1)-3,1):min(remainingTimeVector(1)+3,length(signalTime))]); % find amplitude points for interpolation
        y1 = interp1(xint,yint,x1); % interpolate to find amplitude value for new time x1
        signalTimeIndex = remainingTimeVector(1)-1; 
        x(signalTimeIndex)=x1; % update x-axis with new time after blanking with refractory period
        y(signalTimeIndex)=y1; % update y-axis with interpolated value
        %% 3f. Update indices and flags
        signalTimeIndex = signalTimeIndex-1; % decrement the index as it will be incremented outside the loop
        pulseDetected = true;
    end
    %% 4. Negative Comparator
    if (integratedArea <= nThreshold)
        %% 4a. Calculate area under the curve (fine approximation)
        pulseTime = apcFineIntegration(x1,y1,x2,y2,decayRate,previousIntegratedArea,nThreshold,timeResolution); %fine time at which the threshold is reached
        %% 4b. Capture the pulse time and its polarity
        if (pulseIndex==size(outputPulses,1)) % if size of output pulses exceeds the initial value, then extend its size
            outputPulses=[outputPulses;zeros(1e6,2)];
        end
        pulseIndex=pulseIndex+1;
        outputPulses(pulseIndex,1) = pulseTime;
        outputPulses(pulseIndex,2) = nThreshold;
        %% 4c. Blank for refractory period and reset area
        x1 = pulseTime+refractoryPeriod;
        integratedArea=0;
        %% 4d. Check for the end of the signal (out of bounds)
        if(x1 > signalTime(end))
            outOfBounds = true;
            break;
        end
        %% 4e. Update x-axis and y-axis with new co-ordinates 
        remainingTimeVector = find(signalTime>x1);
        if(isempty(remainingTimeVector))
            break;
        end
        % find y-axis point (y1) by interploation
        xint = [];
        yint = [];
        xint = signalTime([max(remainingTimeVector(1)-5,1):min(remainingTimeVector(1)+5,length(signalTime))]); % find time points for interpolation (max:choose nearby point or the first point, min: choose nearby point or last point)
        yint = signalAmplitude([max(remainingTimeVector(1)-5,1):min(remainingTimeVector(1)+5,length(signalTime))]); % find amplitude points for interpolation
        y1 = interp1(xint,yint,x1); % interpolate to find amplitude value for new time x1
        signalTimeIndex = remainingTimeVector(1)-1;
        x(signalTimeIndex)=x1; % update x-axis with new time after blanking with refractory period
        y(signalTimeIndex)=y1; % update y-axis with interpolated value
        %% 4f. Update indices and flags
        signalTimeIndex = signalTimeIndex-1; % decrement the index as it will be incremented outside the loop
        pulseDetected = true;
    end
    %% 5. Miscellaneous
    %% 5a. Update x-axis and y-axis with new co-ordinates after zero crossing
    if(zeroCrossing && ~pulseDetected)
        x(signalTimeIndex)=x2;
        y(signalTimeIndex)=0;
        signalTimeIndex = signalTimeIndex-1; % decrement the index as it will be incremented outside the loop
        zeroCrossing = false;
    end
    signalTimeIndex=signalTimeIndex+1;
    %% 5b. Stop if "out of bounds" is detected within the loop
    if(outOfBounds)
        break;
    end
end
outputPulses(pulseIndex+1:end,:)=[]; % remove the unused initialized values
disp(sprintf(['%s', '%d', '%s'], 'Conversion complete: ',pulseIndex,' pulses generated'));