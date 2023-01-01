function [P_Pc,M,A_At,T_Tc] = normal_shock(k,Pb_Pc,Ae_At,Ac_At,N)
% NORMAL_SHOCK Computes isentropic flow properties at stations of a nozzle
%   when a normal shock exists in the diverging section.
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

% location of shock
Ans_At = fzero(@(guess)NS_rootfun(k,Ae_At,Pb_Pc,guess),[1,Ae_At]);

% Mach number before the normal shock
M1 = AMR(k,Ans_At,'sup');

% normal Mach and stagnation pressure drop
[P02_Pc,M2] = NSR(k,M1);

% normal shock area to new sonic area ratio
Ans_As2 = MAR(k,M2);

% throat to new sonic area ratio
At_As2 = Ans_As2/Ans_At;

% converging section (subsonic choked)
for idx = 1:N
    M(idx) = AMR(k,A_At(idx),'sub');
    T_Tc(idx) = MTR(k,M(idx));
    P_Pc(idx) = MPR(k,M(idx));
end

% diverging section (with normal shock)
for idx = N + (1:N)
    if A_At(idx) < Ans_At
        M(idx) = AMR(k,A_At(idx),'sup');
        T_Tc(idx) = MTR(k,M(idx));
        P_Pc(idx) = MPR(k,M(idx));
    else
        M(idx) = AMR(k,At_As2*A_At(idx),'sub');
        T_Tc(idx) = MTR(k,M(idx));
        P_P02 = MPR(k,M(idx));
        P_Pc(idx) = P_P02*P02_Pc;
    end
end


end