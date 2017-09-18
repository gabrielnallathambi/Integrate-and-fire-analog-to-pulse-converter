function pulseTime = apcFineIntegration(x1,y1,x2,y2,decayRate,previousArea,threshold,timeResolution)
%APCFINEINTEGRATION Fine Integrator of APC converter and computes time occurence
% of output pulse with resolution of timeResolution.
% apcFineIntegration.m : Divides time axis into cummulatively increasing chunks 
% specified by time resolution, uses "apcIntegration.m" to calculate the
% area area of the chunks, and finds the exact time at which the comparator 
% threshold is crossed  
%
% Inputs:
% x1, x2 - first and second x co-ordinates (scalar)
% y1, y2 - first and second y co-ordinates (scalar)
% decayRate - rate of decay of Analog to pulse converter
% previousArea - previous area under curve 
% threshold - threshold of APC
% timeResolution - time resolution of integration
%
% Outputs:
% pulseTime -  time occurence of output pulse
%
%   % Example:
%   %   refer function call in apcTimeApprox.m 
% Last Updated : 1-18-2017
% References:
%     [1] Alexander singh alvarado, TIME ENCODED COMPRESSION AND 
% CLASSIFICATION USING THE INTEGRATE AND FIRE SAMPLER Jan 2012, University 
% of Florida 


N = floor((x2-x1)/timeResolution); % total # of fine intervals for integration
if ((N ==0)) % check if interval is too small
    pulseTime=x2;
else
    intervalVector = 1:N;
    intervalWidth = (x2-x1)/N; % width of each interval
    m = (y2-y1)/(x2-x1); % slope and intercept for interploating y2
    b= y2 - m*x2;
    x2FineTimeVectors = x1 + intervalWidth*intervalVector; % cummulatively increasing intervals
    y2FineAmplitudeVectors = (x2FineTimeVectors)*m+b;
    % calculate areas of the fine intervals
    integratedAreas = apcIntegration(x1*ones(N,1),x2FineTimeVectors(1:end)',y1*ones(N,1),y2FineAmplitudeVectors(1:end)',decayRate,previousArea); 
    errors = threshold-integratedAreas; 
    [~,ind]=min(abs(errors)); % find the interval at which threshold is crossed
    pulseTime = x2FineTimeVectors(ind(1));
end
