function [P_Pc,M,A_At,T_Tc] = subsonic(k,Pb_Pc,Ae_At,Ac_At,N)
% SUBSONIC Computes isentropic flow properties at stations of a nozzle for
%   when the flow is completely subsonic.
%
% Inputs:
%   k - Ratio of specific heats
%   Pb_Pc - Back to chamber (total) pressure ratio
%   Ae_At - Nozzle area expansion ratio
%   Ac_At - Chamber area contraction ratio
%   N - Number of stations
%
% Outputs:
%   P_Pc - Station to chamber pressure ratios
%   M - Station Mach numbers
%   A_At - Station area ratios
%   T_Tc - Station to chamber temperature ratios
%
% Notes:
%   Author - Matthew Mader
%   Contact - matthewjmader@gmail.com
%

% initializations
A_At = linspace(Ac_At,1,N+1);
A_At = [A_At(1:end-1) linspace(1,Ae_At,N)];
M = zeros(size(A_At));
P_Pc = zeros(size(A_At));
T_Tc = zeros(size(A_At));

% compute the sonic flow area (not throat area for this case)
Me = PMR(k,Pb_Pc);
Ae_As = MAR(k,Me);

% convert area ratios
A_As1 = A_At*Ae_As/Ae_At;

% ratio cannot be lower than 1 (rounding errors)
A_As1 = max(A_As1,1);

for idx = 1:length(A_At)
    % Evaluate Mach and pressure ratio at each sub-sonic area station
    M(idx) = AMR(k,A_As1(idx),'sub');
    T_Tc(idx) = MTR(k,M(idx));
    P_Pc(idx) = MPR(k,M(idx));
end % for

end