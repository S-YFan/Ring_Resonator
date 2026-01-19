clc; clear; close all;

%20250102, This code demonstration no time modulated single resonator case, by varing the coupling coefficient we can understand the different dip about the resonator transmittance. 

omega_0 = 1550;
omega = linspace(omega_0-10,omega_0+10,201);
tau_w = 12;
tau_inner1 = 12;

transfer = 1-(2./tau_w)./(1i.*(omega-omega_0)+ 1./tau_w+ 1 ./tau_inner1);

figure;
plot(omega,abs(transfer).^2);
xlabel('omega');
ylabel('intensity');
grid on;
title('Transmission');
