function integratedArea = apcIntegration( t1,t2,y1,y2,decayRate,previousArea)
%APCINTEGRATION Integrator of APC converter.
% apcIntegration.m : Implements the integrator equation of Analog to Pulse 
% converter 
%
% Inputs:
% t1, t2 - first and second x co-ordinates (maybe scalar or vector)
% y1, y2 - first and second y co-ordinates (maybe scalar or vector)
% decayRate - rate of decay of Analog to pulse converter
% previousArea - previous area under curve 
%    
% Outputs:
% integratedArea - total area under curve
%
%   % Example:
%   %   refer function call in apcTimeApprox.m (scalar inputs) and
%   apcFineIntegration.m (vector inputs)
%
% Last Updated : 1-18-2017
% References:
%     [1] Alexander singh alvarado, TIME ENCODED COMPRESSION AND 
% CLASSIFICATION USING THE INTEGRATE AND FIRE SAMPLER Jan 2012, University 
% of Florida.

%% Check if t1 and t2 are too close during sub-routine from "fineIntegration.m"
if(t2(1)-t1(1) < 1e-20)
    integratedArea =zeros(length(t1),1);
else
    %% calculate equation of straight line between t1 and t2
    x1 = 0; % shifted things to the origin
    x2 = t2-t1;
    m  = (y2-y1)./(x2-x1);
    b   = y2 - m.*x2; % equation of straight line between t1 and t2 is y=m(t-t1)+b
    %% Area calculation
    if decayRate==0
        % Integrate from t1 to t2 [m(t-t1)+b].dt,
        integratedArea =((m/2).*(t2.^2-t1.^2))-(m.*t1.*(t2-t1))+(b.*(t2-t1))+previousArea;
    else
        % Integrate from t1 to t2 [{m(t-t1)+b}exp(-a(t2-t)).dt]
        integratedArea = (1/decayRate).*exp(-1*decayRate*(t2-t1)).*((m/decayRate)-b)+(1/decayRate).*((m.*(t2-t1))+b-(m/decayRate))+ previousArea.*exp(-1*decayRate*(t2-t1));
    end
    
end