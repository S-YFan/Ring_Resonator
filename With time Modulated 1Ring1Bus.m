clc; clear; close all;

small_omega = 2*pi*299792458/(1542.07*10^-9);
Freq_FSR =2*pi*15.04*10^6;
%Paper experiment setup, from Dutt, et al. Nat. Comm. (2019)

J = 0.116 * Freq_FSR;
gamma=pi*300*1000; 
%EXPERIMENT VALUE IS γ=π⋅300kHz 
%THEORITICAL THIS VALUE IS tau=J/10.
%Three different case, J=0,0.116,0.233

% NN means simulate mode's number
NN = 5;

matrix_C = eye(2);
matrix_D = zeros(2,NN);
matrix_D(1,(NN+1)/2)= 1i * sqrt(gamma);
matrix_D(2,(NN+1)/2)= 1i * sqrt(gamma);


m_CapOmega = diag((-2:2) * Freq_FSR);
%We interest in FSR varing value, so we ignore the background freq. We
%conduct the background in further calculation.
m_CapitalGamma = -gamma*eye(NN);

time_row = linspace(-pi,pi,1000);
omega_row = linspace(-2*Freq_FSR,2*Freq_FSR,1000);

%This include a whole period.
%We need to know that the different TIME will have different Coupling
%coeffiecient, so we write the kappa matrix in the loop while I GUESS that
%will not affect the value because the final value is the intensity.

a_field =zeros(5,length(time_row),length(omega_row));

a_field_2_to_1 =zeros(5,length(time_row),length(omega_row));


Transmission_Tminus2_0 = zeros(length(time_row),length(omega_row));
Transmission_Tminus1_0 = zeros(length(time_row),length(omega_row));
Transmission_T0_0 = zeros(length(time_row),length(omega_row));
Transmission_Tplus1_0 = zeros(length(time_row),length(omega_row));
Transmission_Tplus2_0 = zeros(length(time_row),length(omega_row));

Transmission_Case2_0 = zeros(length(time_row),length(omega_row));
Transmission_Case2_2_to_1 = zeros(length(time_row),length(omega_row));



for index_time=1:length(time_row)
   epsilon_k = J * (cos(time_row(index_time)));
   m_CapitalGamma = -gamma*eye(NN) + 1i * epsilon_k;
   m_kappaTranspose = zeros(NN,2);
   m_kappaTranspose((NN+1)/2,1) = 1i * sqrt(gamma) * exp(-1i * small_omega * time_row(index_time));
   m_kappaTranspose((NN+1)/2,2) = 1i * sqrt(gamma) * exp(-1i * small_omega * time_row(index_time));
   s_in_one_to_two = [exp(-1i * small_omega * time_row(index_time)); 0];
   s_in_two_to_one = [0; exp(-1i * small_omega * time_row(index_time))];


   for CouplingIndex=1:NN-1
      %J1,J-1
      m_CapitalGamma(CouplingIndex,CouplingIndex+1)=1i*epsilon_k;
      m_CapitalGamma(CouplingIndex+1,CouplingIndex)=1i*epsilon_k;
   end
   for CouplingIndex=1:NN-2
      %J2,J-2
      m_CapitalGamma(CouplingIndex,CouplingIndex+2)=1i*epsilon_k;
      m_CapitalGamma(CouplingIndex+2,CouplingIndex)=1i*epsilon_k;
   end
   for CouplingIndex=1:NN-3
      %J3,J-3
      m_CapitalGamma(CouplingIndex,CouplingIndex+3)=1i*epsilon_k;
      m_CapitalGamma(CouplingIndex+3,CouplingIndex)=1i*epsilon_k;
   end
   for CouplingIndex=1:NN-4
      %J4,J-4
      m_CapitalGamma(CouplingIndex,CouplingIndex+4)=1i*epsilon_k;
      m_CapitalGamma(CouplingIndex+4,CouplingIndex)=1i*epsilon_k;
   end
   for index_omega=1:length(omega_row)
      m_SmallOmega = eye(5,5)*omega_row(index_omega);
      a_field(:, index_time, index_omega) = (1i*(m_SmallOmega-m_CapOmega)-m_CapitalGamma) \ (m_kappaTranspose * s_in_one_to_two);      

      Transmission_Case2_0(index_time, index_omega) = abs(sum(1i* sqrt(2*gamma) * a_field(:, index_time, index_omega))/ exp(-1i * small_omega * time_row(index_time))).*2;
      
      Transmission_Tminus2_0(index_time, index_omega) = abs(1i*sqrt(2*gamma)*a_field(1,index_time,index_omega) / exp(-1i * small_omega * time_row(index_time))).*2;
      Transmission_Tminus1_0(index_time, index_omega) = abs(1i*sqrt(2*gamma)*a_field(2,index_time,index_omega) / exp(-1i * small_omega * time_row(index_time))).*2;
      Transmission_T0_0(index_time, index_omega) = abs(1i*sqrt(2*gamma)*a_field(3,index_time,index_omega) / exp(-1i * small_omega * time_row(index_time))).*2;
      Transmission_Tplus1_0(index_time, index_omega) = abs(1i*sqrt(2*gamma)*a_field(4,index_time,index_omega) / exp(-1i * small_omega * time_row(index_time))).*2;
      Transmission_Tplus2_0(index_time, index_omega) = abs(1i*sqrt(2*gamma)*a_field(4,index_time,index_omega) / exp(-1i * small_omega * time_row(index_time))).*2;

      a_field_2_to_1(:, index_time, index_omega) = (1i*(m_SmallOmega-m_CapOmega)-m_CapitalGamma) \ (m_kappaTranspose * s_in_two_to_one);
      Transmission_Case2_2_to_1(index_time, index_omega) = abs(sum(1i* sqrt(2*gamma) * a_field_2_to_1(:, index_time, index_omega))/ exp(-1i * small_omega * time_row(index_time))).*2;


   end
end

% case 1, means the different mode's transmission.


figure;
imagesc(time_row,omega_row,Transmission_Tminus2_0');
set(gca,'YDir','normal');
xlabel('time');
ylabel('frequency (Hz)');
title('-2 mode J=', J/Freq_FSR);
grid on;
colormap('plasma');
colorbar;

figure;
plot(mean(Transmission_Tminus2_0, 1), omega_row);
xlabel('intensity');
ylabel('omega');
grid on;
title('time average T -20 mode J=', J/Freq_FSR);


figure;
imagesc(time_row,omega_row,Transmission_Tminus1_0');
set(gca,'YDir','normal');
xlabel('time');
ylabel('frequency (Hz)');
title('-1 mode J=', J/Freq_FSR);
grid on;
colormap('plasma');
colorbar;

figure;
plot(mean(Transmission_Tminus1_0, 1), omega_row);
xlabel('intensity');
ylabel('omega');
grid on;
title('time average T -10 mode J=', J/Freq_FSR);


figure;
imagesc(time_row,omega_row,Transmission_T0_0');
set(gca,'YDir','normal');
xlabel('time');
ylabel('frequency (Hz)');
title('0 mode J=', J/Freq_FSR);
grid on;
colormap('plasma');
colorbar;

figure;
plot(mean(Transmission_T0_0, 1), omega_row);
xlabel('intensity');
ylabel('omega');
grid on;
title('time average T 00 mode J=', J/Freq_FSR);



figure;
imagesc(time_row,omega_row,Transmission_Tplus1_0');
set(gca,'YDir','normal');
xlabel('time');
ylabel('frequency (Hz)');
title('1 mode J=', J/Freq_FSR);
grid on;
colormap('plasma');
colorbar;

figure;
plot(mean(Transmission_Tplus1_0, 1), omega_row);
xlabel('intensity');
ylabel('omega');
grid on;
title('time average T +10 mode J=', J/Freq_FSR);





figure;
imagesc(time_row,omega_row,Transmission_Tplus2_0');
set(gca,'YDir','normal');
xlabel('time');
ylabel('frequency (Hz)');
title('2 mode J=', J/Freq_FSR);
grid on;
colormap('plasma');
colorbar;

figure;
plot(mean(Transmission_Tplus2_0, 1), omega_row);
xlabel('intensity');
ylabel('omega');
grid on;
title('time average T +20 mode J=', J/Freq_FSR);

% case 2, means the physically received transmission.

figure;
imagesc(time_row,omega_row,Transmission_Case2_0');
set(gca,'YDir','normal');
xlabel('time');
ylabel('frequency (Hz)');
title('case 2 sum fisrt and abs later J=', J/Freq_FSR);
grid on;
colormap('plasma');
colorbar;

figure;
plot(mean(Transmission_Case2_0, 1), omega_row, 'b');
hold on;
plot(mean(Transmission_Case2_2_to_1, 1), omega_row, 'r--');
hold off;
xlabel('intensity');
ylabel('omega');
legend('T12','T21');
grid on;
title('case 2 time average T sum  J=', J/Freq_FSR);

% case 3, means add these different mode's transfer function first, then using |tt*| to get its transmission later.


figure;
imagesc(time_row,omega_row,(Transmission_Tplus2_0+Transmission_Tplus1_0+Transmission_T0_0+Transmission_Tminus1_0+Transmission_Tminus2_0)');
set(gca,'YDir','normal');
xlabel('time');
ylabel('frequency (Hz)');
title('sum T-20+T-10+T00T+T+10+T+20 mode J=', J/Freq_FSR);
grid on;
colormap('plasma');
colorbar;


figure;
plot(mean(Transmission_Tplus2_0+Transmission_Tplus1_0+Transmission_T0_0+Transmission_Tminus1_0+Transmission_Tminus2_0, 1), omega_row);
xlabel('intensity');
ylabel('omega');
grid on;
title('case 3 time average T sum  J=', J/Freq_FSR);
