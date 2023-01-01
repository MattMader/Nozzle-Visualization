function [m_dot, Pe, ue] = nozzle(Pc,Pb,k,epsilon,Tc,R,A)
% NOZZLE - Computes ideal nozzle properties given various parameters
%
% Inputs:
%   Pc - Chamber (stagnation) pressure
%   Pb - Back (ambient) pressure
%   k - Ratio of specific heats
%   epsilon - Nozzle area expansion ratio
%   Tc - Chamber (stagnation) temperature
%   R - Specific gas constant
%   A - Throat area
%
% Outputs:
%   m_dot - Mass flow rate
%   Pe - Exit pressure
%   ue - Exit velocity
%
% Notes:
%   Author: Matthew Mader
%   Contact: matthewjmader@gmail.com
%

% catch equilibrium and negative pressures
if Pb/Pc >= 1
    m_dot = 0;
    Pe = Pb;
    ue = 0;
    return
end % if

% subsonic diffuser threshhold
Me_sub = AMR(k,epsilon,'sub');
Pe_Pc_sub = MPR(k,Me_sub);

% perfectly expanded supersonic flow
Me_sup = AMR(k,epsilon,'sup');
Pe_Pc_sup = MPR(k,Me_sup);

% normal shock at exit
[P02_Pc,Me_nse] = NSR(k,Me_sup);
Pe_P02 = MPR(k,Me_nse);
Pe_Pc_nse = Pe_P02*P02_Pc;

% flow is subsonic and unchoked (except edge case)
if Pb/Pc >= Pe_Pc_sub

    % exit pressure equals back pressure
    Pe = Pb;

    % pressure Mach relation to get exit Mach
    Me = PMR(k,Pe/Pc);

    % Mach temperature relation to get exit temperature
    Te = Tc*MTR(k,Me);

    % using speed of sound for calorically perfect gas
    ue = Me*sqrt(k*R*Te);

    % obtain exit to sonic area ratio from area Mach relation
    Ae_As = MAR(k,Me);

    % compute throat to sonic area ratio and bound minimum to 1
    At_As = max(Ae_As/epsilon,1); % (Ae/As)*(At/Ae) => (At/As)
    
    % Mach at throat from area Mach relations
    Mt = AMR(k,At_As,'sub');

    % mass flow with specific throat Mach
    m_dot = mass_flow(k,R,A,Tc,Pc,Mt);

% flow is choked but normal shock between throat and exit
elseif Pb/Pc >= Pe_Pc_nse

    % exit pressure normalizes
    Pe = Pb;

    % compute choked mass flow
    m_dot = mass_flow(k,R,A,Tc,Pc);

    % use root finding to locate (area ratio) normal shock in nozzle
    Ans_At = fzero(@(Ans_At)NS_rootfun(k,epsilon,Pb/Pc,Ans_At),[1 epsilon]);
    [~,Me] = NS_rootfun(k,epsilon,Pb,Ans_At);

    % Mach temperature relation to get exit temperature
    Te = Tc*MTR(k,Me);

    % using speed of sound for calorically perfect gas
    ue = Me*sqrt(k*R*Te);

% flow is choked with no shocks within the nozzle
else

    % use previously computed supersonic exit pressure ratio
    Pe = Pc*Pe_Pc_sup;

    % choked mass flow
    m_dot = mass_flow(k,R,A,Tc,Pc);

    % exit temperature from Mach temperature relation
    Te = Tc*MTR(k,Me_sup);

    % using speed of sound for calorically perfect gas
    ue = Me_sup*sqrt(k*R*Te);

end % if/elseif/elseif/else

end % function