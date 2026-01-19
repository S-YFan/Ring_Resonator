clc; clear; close all;

%plot 2Ring1Bus With Coupled EachOther No Modulated

lambda_0 = 1550e-9; % central wavelength (m)
c = 3e8;
omega_0 = 2*pi*c / lambda_0; % change to angular frequency (rad/s)
lambda = linspace(lambda_0-0.3e-9,lambda_0+0.3e-9,201);
omega = 2*pi*c ./ lambda;

Q_c = 2*1e4;
Q_w = 4*1e4;
Q_i = 2*1e5;

tau_inner1 = 2 * Q_i / omega_0;
tau_w = 2 * Q_w / omega_0;
mu = omega_0/(2*Q_c);
phi = -pi/2;

%version one
A_1 = 1i.*(omega-omega_0)+1/tau_w+1/tau_inner1;
B = (exp(1i.*phi).*A_1 - (2./tau_w).*exp(1i.*phi) - 1i.*mu) ./ (A_1 - 1i.*mu.*exp(1i.*phi));
C = A_1+1i.*mu*B;

transfer_me = exp(1i.*phi) .* (1 - (2./tau_w)./C) - (2./tau_w).*B./C; %verson one transfer function

%Derived from the matrix
C_matrix = [-sqrt(2/tau_w)*exp(1i*phi), -sqrt(2/tau_w)];
K_matrix = [sqrt(2/tau_w); sqrt(2/tau_w)*exp(1i*phi)];
D_matrix = exp(1i*phi);
M_matrix = zeros(2, 2, length(omega));
M_matrix(1,1,:) = 1i*(omega_0 - omega) - 1/tau_w - 1/tau_inner1;
M_matrix(1,2,:) = -1i*mu * ones(length(omega), 1);
M_matrix(2,1,:) = (-1i*mu - (2/tau_w)*exp(1i*phi)) * ones(length(omega), 1);
M_matrix(2,2,:) = 1i*(omega_0 - omega) - 1/tau_w - 1/tau_inner1;

transfer_matrix = zeros(size(omega));
for i=1:length(omega)
transfer_matrix(i) = C_matrix * (M_matrix(:,:,i) \ (-K_matrix)) + D_matrix;
end

%This transfer function is from Liang, et al. paper
delta = (omega - omega_0) ./ omega_0;
gamma0 = 1 ./ (2 * Q_w .* (1i * delta + 1/(2 * Q_i) + 1/(2 * Q_w)));
num = 1 - 2*gamma0 - 2*gamma0 + (4 + 1i*2*(Q_w/Q_c)*exp(-1i*phi) + (Q_w/Q_c)^2) .* gamma0.^2;
den = 1 - 1i*(Q_w/Q_c) .* (2*exp(1i*phi) + 1i*(Q_w/Q_c)) .* gamma0.^2;

transfer_paper_T5 = exp(1i*phi) .* (num ./ den);

figure;
plot(lambda,abs(transfer_matrix).^2, 'b');
hold on;
plot(lambda,abs(transfer_me).^2,'r--');
plot(lambda,abs(transfer_paper_T5).^2,'g*');
hold off;
%Three line perfectly overlape each other.
xlabel('lambda');
ylabel('intensity');
grid on;
xlim([1549.7e-9,1550.3e-9]);
lgd = legend('transfer_matrix','transfer_me','transfer_paper_T5');
lgd.FontSize = 20; 
lgd.ItemTokenSize = [40, 25]; 
title('Transmission');
