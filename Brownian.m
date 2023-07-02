% Mubinjon Satymov
% Physics 661
% Project 1 - Brownian Motion

clear;
clc;
% 1D Uniform
N = [10;100;1000;2000;3000];

global x;
global t;
xAv = zeros(size(N));
x2Av = zeros(size(N));
for k = 1 : 5: 100
    for i = 1:length(N)
        [xAv(i),x2Av(i)] = Br1DUniform(k,5,N(i));
    end
end
plot(N,x2Av)
title("\Deltax^2 vs N")
xlabel("N")
ylabel("\Deltax^2")
figure,plot(N,xAv)
title("\Deltax vs N")
xlabel("N")
ylabel("\Deltax")
figure,plot(t,x)
title("1D Uniformly distributed Brownian Motion")
ylabel("dx")
xlabel("t")

%% 
% 1D Normal
for k = 1 : 1 :10
    for i=1:length(N)
        [xAv(i),x2Av(i)] = Br1DNormal(k,10,N(i));
    end
end

plot(N,x2Av)
title("\Deltax^2 vs N")
xlabel("N")
ylabel("\Deltax^2")
figure,plot(N,xAv)
title("\Deltax vs N")
xlabel("N")
ylabel("\Deltax")
figure,plot(t,x)
title("1D Normalized Brownian Motion")
ylabel("dx")
xlabel("t")

%% 
 
% 3D Uniform
clear

k = [1;1;10];
N = [100;500;1000];

global x;
global y;
global z;

xAv = zeros(size(N));
x2Av = zeros(size(N));
yAv = zeros(size(N));
y2Av = zeros(size(N));
zAv = zeros(size(N));
z2Av = zeros(size(N));

for i=1:length(N)
    [xAv(i),x2Av(i)] = Br3DUniform(k(i),10,N(i));
    [yAv(i),y2Av(i)] = Br3DUniform(k(i),10,N(i));
    [zAv(i),z2Av(i)] = Br3DUniform(k(i),10,N(i));
end

plot(N,x2Av)
title("\Deltax^2 vs N")
xlabel("N")
ylabel("\Deltax^2")
figure, plot (N,xAv)
title("\Deltax vs N")
xlabel("N")
ylabel("\Deltax")
figure,plot3(x,y,z)
xlabel("X")
ylabel("Y")
zlabel("Z")


%% 
% 3D Normal
clear

N = [100;500;1000];
global x;
global y;
global z;
xAv = zeros(size(N));
x2Av = zeros(size(N));
yAv = zeros(size(N));
y2Av = zeros(size(N));
zAv = zeros(size(N));
z2Av = zeros(size(N));

for k = 1 :5 : 50
    for i=1:length(N)
        [xAv(i),x2Av(i)] = Br3DNormal(k,10,N(i));
        [yAv(i),y2Av(i)] = Br3DNormal(k,10,N(i));
        [zAv(i),z2Av(i)] = Br3DNormal(k,10,N(i));
    end
    
end

plot(N,x2Av)
title("\Deltax^2 vs N")
xlabel("N")
ylabel("\Deltax^2")
figure, plot (N,xAv)
title("\Deltax vs N")
xlabel("N")
ylabel("\Deltax")
figure,plot3(x,y,z)
xlabel("X")
ylabel("Y")
zlabel("Z")
