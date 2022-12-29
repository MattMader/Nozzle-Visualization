function [P_Pc,M,A_At,T_Tc] = over_under_expanded(k,Ae_At,Ac_At,N)

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