% Coupled-Oscillators: Three Mass System

clear; clc;
%Assuming mass = 2kg and Sping Constant k = 10 N/m and Z = 1000N/m;
m = 2; k = 10; Z = 100; k1 = Z; k4 = k1; k2 = k; k3 = k2; m1 = m; m2 = m1; m3 = m2;

dt = 0.001;
T = 10;
t = 0:dt:T;
N = length(t);

x = zeros(N,3);
v = zeros(size(x));

x(1,:) = [1 1 1]; %initial condition
v(1,:) = [0 0 0]; %initial condition

for i = 1:N-1
    a1 = ((-k1)*x(i,1) - k2*(x(i,1) - (x(i,2)))) /m;
    a2 = ((-k2)*(x(i,2)-x(i,1)) - k3*(x(i,2) - x(i,3))) /m;
    a3 = ((-k3)*(x(i,3)- x(i,2)) - k4*x(i,3)) /m;

    v(i+1,1) = v(i,1) + a1*dt;
    v1 = v(i+1,1);
    
    v(i+1,2) = v(i,2) + a2*dt;
    v2 = v(i+1,2);
    
    v(i+1,3) = v(i,3) + a3*dt;
    v3 = v(i+1,3);

    x(i+1,1) = x(i,1) + v(i+1,1)*dt;
    x1 = x(i+1,1);
  
    x(i+1,2) = x(i,2) + v(i+1,2)*dt;
    x2 = x(i+1,2);
   
    x(i+1,3) = x(i,3) + v(i+1,3)*dt;
    x3 = x(i+1,3);

    KE1(i+1) = (1/2)*m1.*(v1).^2;
    KE2(i+1) = (1/2)*m2.*(v2).^2;
    KE3(i+1) = (1/2)*m3.*(v3).^2;

    U1(i+1) = 1/2*(k1*(x1.^2)+k2*0.5*(x2-x1).^2);
    U2(i+1) = 1/2*((k2*0.5*(x2-x1).^2)+(k3*0.5*(x3-x2).^2));
    U3(i+1) = 1/2*((k3*0.5*(x3-x2).^2)+(k4*(x3).^2));
end

KE = KE1 + KE2 + KE3;
U = U1 + U2 + U3;
E = KE  + U;

subplot(3,1,1)
p1 = plot(t,x);
hold on
title('Displacement vs Time')
xlabel("t"); ylabel("x");


p2 = plot(t(1,1),x(1,1),'o');
hold on
p2(1).MarkerFaceColor = p1(1).Color;

p3 = plot(t(1,1),x(1,2),'o');
p3(1).MarkerFaceColor = p1(2).Color;

p4 = plot(t(1,1),x(1,3),'o');
p4(1).MarkerFaceColor = p1(3).Color;
legend('x1','x2','x3','m1','m2','m3')

subplot(3,1,2);
p5 = plot(x(1,1),0,'o',x(1,2),0,'*',x(1,3),0,'o');
p5(1).MarkerFaceColor = p5(1).Color;
p5(2).MarkerFaceColor = p5(2).Color;
p5(3).MarkerFaceColor = p5(3).Color;
xlim([-15 15])

subplot(3,1,3);
p6 = plot(t,KE);
hold on
p7 = plot(t,U);
P8 = plot(t,E);
hold off

for i = 1:N
    pause(1/N)
    p2.XData = t(1,i);
    p2.YData = x(i,1);

    p3.XData = t(1,i);
    p3.YData = x(i,2);

    p4.XData = t(1,i);
    p4.YData = x(i,3);

    p5(1).XData = x(i,1) - 5;
    p5(2).XData = x(i,2);
    p5(3).XData = x(i,3) + 5;
end
