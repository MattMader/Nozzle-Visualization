%% Setup

% clear
clear
close all
clc

% add needed tools
addpath ./tools/Isentropic/

%% Inputs

% ratio of specific heats
k = 1.2; % [-]

% nozzle expansion area
Ae_At = 3; % [-]

% nozzle contraction area
Ac_At = 4; % [-]

% number of area stations
N = 100;

% diverging half angle
alpha_d = deg2rad(15); % [rad]

% converging half angle
alpha_c = deg2rad(45); % [rad]

%% Nozzle Geometry

% normalized length of converging section
Lc_Dt = (sqrt(Ac_At) - 1)/tan(alpha_c); % [-]

% normalized length of diverging section
Ld_Dt = (sqrt(Ae_At) - 1)/tan(alpha_d); % [-]

% convert station points to normalized distance
x1 = linspace(-Lc_Dt,0,N);
x2 = linspace(0,Ld_Dt,N);
x = [x1 x2];

% conical nozzle momentum correction factor
lambda = 0.5*(1+cos(alpha_d));

%% Delineations

% subsonic
M_sub = AMR(k,Ae_At,'sub');
Pb_Pc_sub = MPR(k,M_sub);

% supersonic
Me_sup = AMR(k,Ae_At,'sup');
Pb_Pc_sup = MPR(k,Me_sup);

% shock at nozzle exit
[P02_Pc,Me_nse] = NSR(k,Me_sup);
Pb_P02 = MPR(k,Me_nse);
Pb_Pc_nse = Pb_P02*P02_Pc;

%% Subsonic

Pb_Pc = linspace(1,Pb_Pc_sub,50);

f = figure(1);

for idx = 1:length(Pb_Pc)

    [P_Pc,M,~,T_Tc] = subsonic(k,Pb_Pc(idx),Ae_At,Ac_At,N);

    subplot 311
    hold on
    plot(x,P_Pc,'b-')

    subplot 312
    hold on
    plot(x,M,'b-')

    subplot 313
    hold on
    plot(x,T_Tc,'b-')

end

%% Normal Shock in Nozzle

Pb_Pc = linspace(Pb_Pc_sub,Pb_Pc_nse+1e-3,25);

for idx = 1:length(Pb_Pc)
    
    [P_Pc,M,~,T_Tc] = normal_shock(k,Pb_Pc(idx),Ae_At,Ac_At,N);

    subplot 311
    hold on
    plot(x,P_Pc,'g-')

    subplot 312
    hold on
    plot(x,M,'g-')

    subplot 313
    hold on
    plot(x,T_Tc,'g-')
end

%% Shocks outside nozzle

[P_Pc,M,A_At,T_Tc] = over_under_expanded(k,Ae_At,Ac_At,N);

subplot 311
hold on
plot(x,P_Pc,'r-')
grid on
ylabel("Pressure Ratio [P/P_c]")
qw{1} = plot(nan, 'b-');
qw{2} = plot(nan, 'g-');
qw{3} = plot(nan, 'r-');
legend([qw{:}], [sprintf("Unchoked (1>P_b/P_0>%0.2f)",Pb_Pc_sub) sprintf("Shock In Nozzle (%0.2f>P_b/P_0>%0.2f)",Pb_Pc_sub,Pb_Pc_nse) sprintf("Shock Outside Nozzle (%0.2f>P_b/P_0)",Pb_Pc_nse)], 'location', 'southwest')
title(sprintf("Nozzle Performance (%0.2f deg Conv, %0.2f deg Div)With Varying Backpressure Ratio",rad2deg(alpha_c),rad2deg(alpha_d)))

subplot 312
hold on
plot(x,M,'r-')
hold off
grid on
ylabel("Mach")

subplot 313
hold on
plot(x,T_Tc,'r-')
hold off
grid on
ylabel("Temperature Ratio [T/T_c]")
xlabel("Normalized Distance From Throat [L/D_t]")