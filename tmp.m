%%

clear
clc

%%

k = 1.2;
epsilon = 3;
Tc = 2800;
R = 287;
A = 0.01;

alpha = deg2rad(15);

%%

N = 1e4;

Pb = 101200;

Pc = linspace(0,convpres(600,'psi','pa'),N);

%%

m_dot = zeros(size(Pc));
Pe = zeros(size(Pc));
ue = zeros(size(Pc));
F = zeros(size(Pc));

lambda = 0.5*(1 + cos(alpha));

tic
for idx = 1:length(Pc)
    [m_dot(idx), Pe(idx), ue(idx)] = nozzle(Pc(idx),Pb,k,epsilon,Tc,R,A);
    F(idx) = lambda*0.95*m_dot(idx)*ue(idx) + (Pe(idx)-Pb)*A*epsilon;
end
dt = toc(tic);
dt/N

%%

% figure(1)
% plot(convpres(Pc,'pa','atm'),F)
% grid on