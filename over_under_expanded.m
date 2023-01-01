function [P_Pc,M,A_At,T_Tc] = over_under_expanded(k,Ae_At,Ac_At,N)
% OVER_UNDER_EXPANDED Computes isentropic flow properties at stations of a
%   nozzle when no shocks exist in the nozzle.
%
% Inputs:
%   k - ratio of specific heats
%   Pb_Pc - Back to chamber pressure ratio
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
A_At_conv = linspace(Ac_At,1,N+1);
A_At_div = linspace(1,Ae_At,N);
A_At = [A_At_conv(1:end-1), A_At_div];
M = zeros(size(A_At));
P_Pc = zeros(size(A_At));
T_Tc = zeros(size(A_At));

% converging section (subsonic choked)
for idx = 1:N
    [M(idx),T_Tc(idx),P_Pc(idx)] = flowisentropic(k,A_At(idx),'sub');
end

% diverging section (supersonic choked)
for idx = N + (1:N)
    [M(idx),T_Tc(idx),P_Pc(idx)] = flowisentropic(k,A_At(idx),'sup');
end

end