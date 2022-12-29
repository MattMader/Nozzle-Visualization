function [P_Pc,M,A_At,T_Tc] = subsonic(k,Pb_Pc,Ae_At,Ac_At,N)
% SUBSONIC_UNCHOKED Computes isentropic flow properties at stations of a
%   conical nozzle for when the flow is completely subsonic.
%
% Inputs:
%

% initializations
A_At = linspace(Ac_At,1,N+1);
A_At = [A_At(1:end-1) linspace(1,Ae_At,N)];
M = zeros(size(A_At));
P_Pc = zeros(size(A_At));
T_Tc = zeros(size(A_At));

% compute the sonic flow area (not throat area for this case)
[~,~,~,~,Ae_As] = flowisentropic(k,Pb_Pc,'pres');

% convert area ratios
A_As1 = A_At*Ae_As/Ae_At;

% ratio cannot be lower than 1 (rounding errors)
A_As1 = max(A_As1,1);

for idx = 1:length(A_At)
    % Evaluate Mach and pressure ratio at each sub-sonic area station
    [M(idx),T_Tc(idx),P_Pc(idx)] = flowisentropic(k,A_As1(idx),'sub');
end % for

end