
% Assignment 5

clc; clear; 

%Define time
dt = 0.005;
t = -10:dt:10;
N = length(t);

%Define frequency domain
fs = 1/dt;
df = fs/N;
f = -(N/2-1/2)*df:df:(N/2-1/2)*df;

%% 1. In time domain evaluate and plot the following functions:
%  a) An information-bearing signal.
A_m = 1; f_m = 1; 
m =A_m* cos(2*pi*f_m*t);
txt1 = 'Plot of an Information-Bearing Signal';
plot(t,m), title(txt1);

%% b) A carrier wave.
A_c = 1; f_c = 10; 
c = A_c* cos(2*pi*f_c*t);
txt2 = 'Plot of a Carrier Wave';
plot(t,c), title(txt2);

%% c)  The modulated signal.
A_c = 1; f_c = 10; k_a = 0.5;
s = A_c*(1+k_a.*m).*cos(2*pi*f_c*t);
txt3 = 'Plot of The modulated signal';
plot(t,s), title(txt3);

%% 2. Calculate and plot the frequency domain representation of all the above functions, i.e.,
%  a) M(f)
M = fft(m); M = fftshift(M); M = 1/sqrt((N-1)/2)*M;
figure, plot(f,abs(M)), title('M(f)')

%% 
%  b) C(f)
C = fft(c); C = fftshift(C); C = 1/sqrt((N-1)/2)*C;
figure, plot(f,abs(C)), title('C(f)')

%% 
%  c) S(f)
S = fft(s); S = fftshift(S); S = 1/sqrt((N-1)/2)*S;
figure, plot(f,abs(S)), title('S(f)')