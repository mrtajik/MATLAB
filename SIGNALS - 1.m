
% Assignment 1

% a)
clear;
clc;
t = -5:0.01:5;
g = sin(2*pi*t);
plot(t,g);
figure()

% b)
clear;
clc;
t = -5:0.01:5;
a = 2;
b = 1;
g = exp(-(((t-a)/b).^2));
plot(t,g)
figure()

% c)
clear;
clc;
t = -5:0.01:5;
a = -1;
b = 2;
g = exp(-(((t-a)/b).^2));
plot(t,g)
figure()

% d)
clc;
clear;
t = -5:0.01:5;
syms t
g = piecewise(t < 0,0,t >= 0,exp(-t));
fplot(g)
figure()

% e)
clc;
clear;
t = -5:0.01:5;
g = sin(2*pi*t).*exp(-(t-2).^2);
plot(t,g)
figure()

% f)
clc;
clear;
t = -5:0.01:5;
g = sinc(t);
plot(t,g)
figure()

% g)
clc;
clear;
t = -5:0.01:5;
syms t
g = piecewise(abs(t) > 0.5,0,abs(t) <= 0.5,1);
fplot(g)
figure()

% Another way to plot giving function
% g)
clc;
clear;
a = 1;
t = -5:0.001:5;
rect=@(t,a) ones(1,numel(t)).*(abs(t)<a/2) % a is the width of the pulse
g = rect(t,1);   
plot(t,g)
figure()
