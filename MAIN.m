%% Setup

% clear
clear
close all
clc

%% Inputs

% ratio of specific heats
k = 1.2; % [-]

% nozzle expansion area
Ae_At = 25; % [-]

% nozzle contraction area
Ac_At = 4; % [-]

% number of area stations
N = 100;

% diverging half angle
alpha_d = deg2rad(15); % [rad]

% converging half angle
alpha_c = deg2rad(30); % [rad]

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
[~,~,Pb_Pc_sub] = flowisentropic(k,Ae_At,'sub');

% supersonic
[Me_sup,~,Pb_Pc_sup] = flowisentropic(k,Ae_At,'sup');

% shock at nozzle exit
[~,~,P2_P1,~,~] = flownormalshock(k,Me_sup);
Pb_Pc_nse = P2_P1*Pb_Pc_sup;

%% Subsonic

Pb_Pc = linspace(1,Pb_Pc_sub,50);

for idx = 1:length(Pb_Pc)

    [P_Pc,M,~,T_Tc] = subsonic(k,Pb_Pc(idx),Ae_At,Ac_At,N);

    figure(2)
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

    figure(2)
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

figure(2)
subplot 311
hold on
plot(x,P_Pc,'r-')
grid on
ylabel("Pressure Ratio [P/P_c]")
title("Nozzle Performance With Varying Backpressure Ratio")
qw{1} = plot(nan, 'b-');
qw{2} = plot(nan, 'g-');
qw{3} = plot(nan, 'r-');
legend([qw{:}], ["Unchoked" "Shock In Nozzle" "Shock Outside Nozzle"], 'location', 'southwest')

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