function [f,Me] = NS_rootfun(k,Ae_At,Pb_Pc,Ans_At)
% NS_ROOTFUN For root solver to locate the location of the normal shock
%   within the nozzle and return exit Mach number. Iteratively solve to
%   equilize the exit pressure equaling the back pressure.
%
% Inputs:
%   k - Ratio of specific heats
%   Ae_At - Nozzle area expansion ratio
%   Pb_Pc - Back to total pressure ratio
%   Ans_At - Guess for normal shock station [1, Ae_At]
%
% Outputs:
%   f - Objective for root finding algorithm
%   Me - Exit Mach number
%

% Mach number before the normal shock
M1 = AMR(k,Ans_At,'sup');

% normal shock relation
[P02_Pc,M2] = NSR(k,M1);

% normal shock area to new sonic area ratio
Ans_As2 = MAR(k,M2);

% exit area to new sonic area ratio
Ae_As2 = Ans_As2/Ans_At*Ae_At;

% exit Mach
Me = AMR(k,Ae_As2,'sub');

% exit pressure ratio with stagnation pressure loss
Pe_P02 = MPR(k,Me);

% exit pressure ratio to chamber pressure
f = Pb_Pc - Pe_P02*P02_Pc;

end % function