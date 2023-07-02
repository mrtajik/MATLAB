
% Assignment 2

% a) First, evaluate and plot these functions in the time domain (use −10<= t<=10)
clear;
clc;
%Define time
dt = 0.01;
t = -10:dt:10;
N = length(t);
%Define frequency domain
fs = 1/dt;
df = fs/N;
f = -(N/2-1/2)*df:df:(N/2-1/2)*df;

% i)
T = 2;
syms t
g = piecewise(abs((1/T)*t) > 0.5,0,abs((1/T)*t) <= 0.5,1);
figure, fplot(g)
title('Time function of rect(t/T)')

%% 
% ii)
dt = 0.01;
t = -10:dt:10;
N = length(t);

W = 1;
g1 = sinc(2*W*t);
plot(t,g1)
title('Time function of sinc(2Wt)')

%% 
% iii)
a = 1; 
g2 = exp(-a*t).*heaviside(t);
plot(t,g2);
txt = 'Graph of exp(-at)*u(t)';
title(txt)

%% 
% iv)
a = 1; 
g3 = exp(-a*abs(t));
figure, plot(t,g3);
txt = 'Graph of exp(-a|t|)';
title(txt)

%% 
% v)
a = 1; 
g4 = exp(-pi*(t.^2));
figure, plot(t,g4);
txt = 'Graph of exp(-\pi*t^2)';
title(txt)

%% 
% vi)
for i = 1:length(t)
    g(i) = triangle(t(i));
end
figure,plot(t,g)
title('Graph of Triangular wave')
 
%%
% vii)
for i = 1:length(t)
    g(i) = dirac_delta(t(i));
end
figure, plot(t,g)
title('Graph of Dirac Delta Function')

%% 
% viii)
g = ones(1,length(t));
figure, plot(t,g)
title('Plot of the function 1')

%% 
% ix)
t0 = 2;
for i = 1:length(t)
    g(i) = dirac_delta(t(i)-t0);
end
figure, plot(t,g)
title('Graph of Dirac Delta Function(t-to)')

%% 
% x)
fc = 2; 
g = exp(1i*2*fc*pi*t);
figure, plot(t,g);
txt = 'Graph of exp(j*2\pi*fc*t)';
title(txt)

%% 
% xi)
fc = 2; 
g = cos(2*fc*pi*t);
figure, plot(t,g);
txt = 'Graph of cos(2\pi*fc*t)';
title(txt)

%% 
% xii)
fc = 2; 
g = sin(2*fc*pi*t);
figure, plot(t,g);
txt = 'Graph of sin(2\pi*fc*t)';
title(txt)

%% 
% xiii)
g = sign(t);
figure, plot(t,g)
title('Function of sgn(t)')

%% 
% xiv)
g = (pi*t).^(-1);
figure, plot(t,g)
title('Function of 1/(\pi*t)')

%% 
% xv)
g = heaviside(t);
figure, plot(t,g)
title('Function of u(t)')

%% 
% xvi)

g = 10*ones(length(t));
figure, plot(t,g)
title('Summation of Delta function')

%% b)use the “fft” command to calculate their discrete Fourier transforms and plot
% i)
for i = 1:length(t)
    g(i) = rect(t(i));
end
G = fft(g); G = fftshift(G); G = 1/sqrt((N-1)/2)*G;
plot(f,abs(G))
title('Fourier Transform of rect(t/T)')

%% 
% ii)
W = 1;
g1 = sinc(2*W*t);
G = fft(g1);
plot(t,abs(fftshift(G)))
title('Fourier Transform of sinc(2Wt)')

%% 
% iii)
a = 1; 
g2 = exp(-a*t).*heaviside(t);
G = fft(g2);
plot(t,g2);
txt = 'Fourier Transform of exp(-at)*u(t)';
title(txt)

%% 
% iv)
a = 1; 
g = exp(-a*abs(t));
G = fft(g); G = fftshift(G); G = 1/sqrt((N-1)/2)*G;
figure, plot(t,abs(G));
txt = 'Fourier Transform of exp(-a|t|)';
title(txt)

%% 
% v)
a = 1; 
g = exp(-pi*(t.^2));
G = fft(g); G = fftshift(G); G = 1/sqrt((N-1)/2)*G;
figure, plot(t,abs(G));
txt = 'Fourier Transform of exp(-\pi*t^2)';
title(txt)

%% 
% vi)
for i = 1:length(t)
    g(i) = triangle(t(i));
end
G = fft(g); G = fftshift(G); G = 1/sqrt((N-1)/2)*G;
figure, plot(t,abs(G));
title('Fourier Transform of Triangular wave')
title(txt)

%%
% vii)
for i = 1:length(t)
    g(i) = dirac_delta(t(i));
end
G = fft(g); G = fftshift(G); G = 1/sqrt((N-1)/2)*G;
figure, plot(t,abs(G));
title('Fourier Transform of Dirac Delta Function')

%% 
% viii)
g = ones(1,length(t));
G = fft(g); G = fftshift(G); G = 1/sqrt((N-1)/2)*G;
figure, plot(t,abs(G));
title('Fourier Transform of the function 1')

%% 
% ix)
t0 = 2;
for i = 1:length(t)
    g(i) = dirac_delta(t(i)-t0);
end
G = fft(g); G = fftshift(G); G = 1/sqrt((N-1)/2)*G;
figure, plot(t,abs(G));
title('Fourier Transform of the Dirac Delta Function(t-to)')

%% 
% x)
g = exp(1i*2*fc*pi*t);
G = fft(g); G = fftshift(G); G = 1/sqrt((N-1)/2)*G;
figure, plot(t,abs(G)); txt = 'Fourier Transform of exp(j*2\pi*fc*t)';
title(txt);

%% 
% xi)
g = cos(2*fc*pi*t);
G = fft(g); G = fftshift(G); G = 1/sqrt((N-1)/2)*G;
figure, plot(t,abs(G));
txt = 'Fourier Transform of cos(2\pi*fc*t)';
title(txt)

%% 
% xii)
g = sin(2*fc*pi*t);
G = fft(g); G = fftshift(G); G = 1/sqrt((N-1)/2)*G;
figure, plot(t,abs(G));
txt = 'Fourier Transform of sin(2\pi*fc*t)';
title(txt)

%% 
% xiii)
g = sign(t);
G = fft(g); G = fftshift(G); G = 1/sqrt((N-1)/2)*G;
figure, plot(t,abs(G));
title('Fourier Transform of of sgn(t)')

%% 
% xiv)
g = (pi*t).^(-1);
G = fft(g); G = fftshift(G); G = 1/sqrt((N-1)/2)*G;
figure, plot(t,abs(G));
title('Fourier Transform of 1/(\pi*t)')

%% 
% xv)
g = heaviside(t);
G = fft(g); G = fftshift(G); G = 1/sqrt((N-1)/2)*G;
figure, plot(t,abs(G));
title('Fourier Transform of u(t)')

%% 
% xvi)
g = 10*ones(length(t));
G = fft(g); G = fftshift(G); G = 1/sqrt((N-1)/2)*G;
figure, plot(t,abs(G));
title('Fourier Transform of the Summation of Delta function')

%%
% Defining functions
function gg = triangle(t)
if abs(t) < 1
    gg = 1- abs(t);
else
    gg = 0;
end
end
%% 

function g = dirac_delta(t)
if t == 0
 g = 10;
else
 g = 0;
end
end

%% 

function g = rect(t)
if abs(t/2) <= 1/2
g = 1;
else
g = 0;
end
end

