
% Assignment 3

clc; clear; 

%Define time
dt = 0.01;
t = -10:dt:10;
N = length(t);

%Define frequency domain
fs = 1/dt;
df = fs/N;
f = -(N/2-1/2)*df:df:(N/2-1/2)*df;

%% a) Impulse response
h = exp(-t).*(t>0);
figure, plot(t,h)
title('Impulse response h(t)'), xlabel('t'), ylabel('h(t)');

%% b) Input 
%  Note:There is other way to implement Unit Step Function
x = 1*((0<=t) & (t<=2));
figure, plot(t,x)
title('Input x(t)'), xlabel('t'), ylabel('x(t)');

%% c) Output
y = conv(x,h,'same');
figure,plot(t,y)
title('Output y(t)'), xlabel('t'), ylabel('y(t)');

%% d) Transfer function
H = fft(h); Hf = fftshift(H); Hf = 1/sqrt((N-1)/2)*Hf;
figure, plot(f,abs(Hf))
title('Transfer function H(f)'), xlabel('f'), ylabel('|H(f)|')

%% e) Input Fourier Transform
X = fft(x); Xf = fftshift(X); Xf = 1/sqrt((N-1)/2)*Xf;
figure, plot(f,abs(Xf))
title('Input Fourier Transfer X(f)'), xlabel('f'), ylabel('|X(f)|')

%% f) Output in Fourier Domain
Y = X.*H; Yf = fftshift(Y); Yf = 1/sqrt((N-1)/2)*Yf;
figure, plot(f,abs(Yf))
title('Output in Fourier Domain Y(f)'), xlabel('f'), ylabel('|Y(f)|')

%% g) Inverse Fourier Transform
yf = ifft(Y); yf = ifftshift(yf);
figure, plot(t,yf)
title('Output from Inverse Fourier Transform'), xlabel('t'), ylabel('y(t)')
