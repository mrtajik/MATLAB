% Mubinjon Satymov
% Physics 661
% Final Project 
% Part 3 - 1D Wave Equation



%% 1D DBC
clear; clc;
L = 1; T = 1;
h = 0.001; r = 1; k = r*h;
x = 0:h:L; t = 0:k:T;
BC = "DBC";

u = zeros(length(t),length(x));

% Setting BC
mu = 0.25; sigma = 0.01;
u(1,:) = f(x,mu,sigma);

% Set Initial Velocity of Wave
g2 = f(x,mu+k,sigma);
g1 = f(x,mu-k,sigma);
g = (g2-g1)/(2*k);

% Set Dirichlet BC (Fixed-Ends)
u0 = 0; uL = 0;
u(:,1) = u0; u(:,end) = uL;

for i = 2:length(x)-1
    u(2,i) = 1/2*(r^2*u(1,i+1) + 2*(1-r^2)*u(1,i) + r^2*u(1,i-1)) + k*g(i);
end

for j = 2:length(t)-1
    for i = 2:length(x)-1
       u(j+1,i) = r^2*u(j,i+1) + 2*(1-r^2)*u(j,i) + r^2*u(j,i-1) - u(j-1,i);
    end 
end

% Plotting
s1 = surf(x,t,u);
shading flat
view(2)
colorbar
title({"Boundary Condition: " + upper(BC)});
xlabel("Position")
ylabel("Time")

%% 1D NBC
clear; clc;
c = 1;
l = 4.5;
dt = 0.06;
dx = 0.05;
ldx = 4.5/dx;
x = 1:ldx;
tm = 100;
N = 0.4/dt;
n = 1:8*N;
n0 = 4*N;
En = zeros(tm,ldx);
En0 = zeros(tm,1);
En0(1:8*N) = exp(((-(n-n0).^2)/(2*N^2)));
Ent(1:tm) = 0;
En(1,1) = En0(1);
for p = 2:tm-1
   En(p,1) = En0(p);
   for q = 2:ldx-1
       
       if (q == (ldx/9)*4)||(q == (ldx/9)*5)
           c = sqrt(0.4);
       end
       if ((q > (ldx/9)*4) && (q < (ldx/9)*5))
           c = 0.5;
       end
       if ((q < (ldx/9)*4) || (q > (ldx/9)*5))
           C = 1;
       end
       En(p+1,q) = ((c*dx)/dt)^2*(En(p,q+1)-2*En(p,q)+En(p,q-1)) + 2*En(p,q) - En(p-1,q); 
   end
   plot(En(p,:))
   title(['Boundary Condition: NBC, Time = ',num2str(p+1)]);
   ylim([-1.5 1.5])
   pause(0.03);
end

%% 1D ABC
clear; clc;
L = 1; T = 1;
h = 0.001; r = 1; k = r*h;
x = 0:h:L; t = 0:k:T;
BC = "ABC";

u = zeros(length(t),length(x));

% Setting BC
mu = 0.25; sigma = 0.01;
u(1,:) = f(x,mu,sigma);

% Set Initial Velocity of Wave
g2 = f(x,mu+k,sigma);
g1 = f(x,mu-k,sigma);
g = (g2-g1)/(2*k);

% Set Dirichlet BC (Fixed-Ends)
u0 = 0; uL = 0;
u(:,1) = u0; u(:,end) = uL;

for i = 2:length(x)-1
    u(2,i) = 1/2*(r^2*u(1,i+1) + 2*(1-r^2)*u(1,i) + r^2*u(1,i-1)) + k*g(i);
end

for j = 2:length(t)-1
    % Can you replace inner loop with matrix calculations?
    for i = 2:length(x)-1
       u(j+1,i) = r^2*u(j,i+1) + 2*(1-r^2)*u(j,i) + r^2*u(j,i-1) - u(j-1,i);
    end 
  
    % Absorbing Boundaries
    u(j+1,1) = u(j,1)+(r-1/r+1)*(u(j+1,1)-u(j));
    u(j+1,end) = u(j,end-1)+(r-1/r+1)*(u(j+1,end-1)-u(j,end));
   
end

 p1 = plot(x,u(1,:));
        t1 = title({"Boundary Condition: " + upper(BC);...
            "\itct = " + compose("%5.3f",t(1))});
        ylim(2.2*max(abs(u(1,:)))*[-1 1])

        for p = 2:length(t)
            if mod(p,5) == 0
                pause(0.005)
                p1.YData = u(p,:);
                t1.String = {"Boundary Condition: " + upper(BC);...
                    "\itct = " + compose("%5.3f",t(p))};
                drawnow
            end
        end

        p1.YData = u(p,:);
        t1.String = {"Boundary Condition: " + upper(BC);...
            "\itct = " + compose("%5.3f",t(p))};

% Helping Function
function u0 = f(x,mu,sigma)
u0 = normpdf(x,mu,sigma) + 0.5*normpdf(x,1-mu,sigma);
end