function [P_Pc,M,A_At,T_Tc] = normal_shock(k,Pb_Pc,Ae_At,Ac_At,N)
% NORMAL_SHOCK solves for the normal shock location within the nozzle and
% returns the exit Mach and pressure ratio.
%
% Inputs:
%

% initializations
A_At_conv = linspace(Ac_At,1,N+1);
A_At_div = linspace(1,Ae_At,N);
A_At = [A_At_conv(1:end-1), A_At_div];
M = zeros(size(A_At));
P_Pc = zeros(size(A_At));
T_Tc = zeros(size(A_At));

% location of shock
Ans_At = fzero(@(guess)rootfun(k,Ae_At,Pb_Pc,guess),[1,Ae_At]);

% Mach number before the normal shock
M1 = flowisentropic(k,Ans_At,'sup');

% normal Mach and stagnation pressure drop
[~,~,~,~,M2,P02_Pc] = flownormalshock(k, M1);

% normal shock area to new sonic area ratio
[~,~,~,~,Ans_As2] = flowisentropic(k, M2);

% throat to new sonic area ratio
At_As2 = Ans_As2/Ans_At;

% converging section (subsonic choked)
for idx = 1:N
    [M(idx),T_Tc(idx),P_Pc(idx)] = flowisentropic(k,A_At(idx),'sub');
end

% diverging section (with normal shock)
for idx = N + (1:N)
    if A_At(idx) < Ans_At
        [M(idx),T_Tc(idx),P_Pc(idx)] = flowisentropic(k,A_At(idx),'sup');
    else
        [M(idx),T_Tc(idx),P_P02] = flowisentropic(k,At_As2*A_At(idx),'sub');
        P_Pc(idx) = P_P02*P02_Pc;
    end
end


end

function f = rootfun(k,Ae_At,Pb_Pc,Ans_At)
% ROOTFUN For root solver to locate the location of the normal shock within
% the nozzle. Iteratively solves to equalize the exit pressure to the back
% pressure.

% Mach number before the normal shock
M1 = flowisentropic(k,Ans_At,'sup');

% normal shock
[~,~,~,~,M2,P02_Pc] = flownormalshock(k, M1);

% normal shock area to new sonic area ratio
[~,~,~,~,Ans_As2] = flowisentropic(k, M2);

% exit area to new sonic area ratio
Ae_As2 = Ans_As2/Ans_At*Ae_At;

% exit pressure ratio with stagnation pressure loss
[~,~,Pe_P02] = flowisentropic(k,Ae_As2,'sub');

% exit pressure ratio to chamber pressure
f = Pb_Pc - Pe_P02*P02_Pc;


end % function