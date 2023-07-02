% Kramers-Kronig Relations
% Refractive index = n + i*k;

clear;
clc;

% Importing data from txt file
A = importdata('extinctionSpectrum.txt');
wavelength = A.data(:,1);     % Wavelength data imported
K = A.data(:,2);              % Extinction coefficient data imported

omega_P = wavelength;   
omega = omega_P; 


B = ones(size(wavelength)); % Creating Array of ONEs
N = length(wavelength);      
d_omegaP = omega_P(1:end-1)-omega_P(2:end); %Derivative of the Omega_prime

for j = 1:N
    for i = 1:N-1
        if i ~= j % Condition to implement integration
            B(j) = B(j)+[((2/pi)*(omega_P(i).*K(i)))./((omega_P(i).^2 - omega(j).^2))].*d_omegaP(i);  
        end
    end
end

% Output figure in full screen mode
fh = figure();
fh.WindowState = 'maximized';

% Plot
plot(wavelength,[K B])
legend('Imaginary','Real')
title('Plot of the spectra for the real and imaginary part of the refractive index')
xlabel('Wavelength')
ylabel("Refractive Index") 

