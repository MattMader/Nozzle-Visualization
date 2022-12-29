function [m_dot,P_e,u_e] = nozzle(k,epsilon,Tc,Pc,Pa)

%% Constants
Pb_Pc = Pa/Pc;

%% Delineations

% subsonic
[~,~,Pb_Pc_sub] = flowisentropic(k,epsilon,'sub');

% supersonic
[Me_sup,~,Pb_Pc_sup] = flowisentropic(k,epsilon,'sup');

% shock at nozzle exit
[~,~,P2_P1,~,~] = flownormalshock(k,Me_sup);
Pb_Pc_nse = P2_P1*Pb_Pc_sup;

end