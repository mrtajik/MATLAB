% 2D Heat Equation Simulation

clear; clc;
L = 1; W = 1; T = 0.05;

h = 0.01; r = 1/4; dt = r*h^2;

% Setting Boundary Conditions
x = -L/2:h:L/2;
y = -W/2:h:W/2;
t = 0:dt:T;

[X,Y] = meshgrid(x,y);

% Preallocate varaible for solution
u = zeros(length(y),length(x),length(t));

% Rectangular
ind = abs(X) < 0.25 & abs(Y) < 0.25;

u(:,:,1) = ind + u(:,:,1);

BC = "NBC"; 

s1 = surf(x,y,u(:,:,1));
shading flat
view(3)
colorbar
t1 = title({"Boundary Condition: " + upper(BC);...
    "\alpha \ t = 0"});
zlim([0 1])
clim([0 1])

%% Loop over time steps
for p = 1:length(t)-1
    for i = 2:length(y)-1
        for j = 2:length(x)-1
            u(i,j,p+1) = (1-4*r)*u(i,j,p) + r*u(i-1,j,p) + r*u(i+1,j,p) + r*u(i,j-1,p) + r*u(i,j+1,p);
        end
    end
   
    % Insert NBC
    if upper(BC) == "NBC"
        % sides
        u(:,1,1:end) = u(:,2,:);
        u(:,length(x),1:end) = u(:,length(x)-1,:);
        u(1,:,1:end) = u(2,:,:);
        u(length(y),:,1:end) = u(length(y)-1,:,:);

        % corners
        u(1,1,1:end) = u(1,1,:);
        u(1,length(x),1:end) = u(1,length(x),:);
        u(length(y),1,1:end) = u(length(y),1,:);
        u(length(y),length(x),1:end) = u(length(y),length(x),:);

    end

    if mod(p,0.001*T/dt) == 0
        s1.ZData = u(:,:,p);
        t1.String = {"Boundary Condition: " + upper(BC);...
            "\alpha\itt = " + compose("%5.3f",t(p))};
        drawnow
    end
end

% Plots last calculation
s1.ZData = u(:,:,end);
t1.String = {"Boundary Condition: " + upper(BC);...
    "\alpha\itt = " + compose("%5.3f",t(end))};
