% Mubinjon Satymov
% Physics 661
% Final Project 
% Part 3 - 2D Laplace Equation

%% 2D Laplace Equation
clear; clc;
L = 1; W = 1;
h = 1e-2; r = 1; k = r*h;

x = -L/2:h:L/2; y = -W/2:k:W/2;

numPoints =(length(x)-1)*(length(y)-1);
maxIter = 10*numPoints;

u = zeros(length(y),length(x),2);

% Setting Boundaries cosine
u(1,:,1) =sin(pi*x/L);
u(end,:,1) = sin(pi*x/L);

% Initialize BC for next iteration
u(:,:,2) = u(:,:,1);

% Create surf plot
s1=surf(x,y,u(:,:,1),'EdgeColor','none');       
shading flat
title('2-D Laplace''s equation');
xlabel('Spatial coordinate (x) \rightarrow')
ylabel('{\leftarrow} Spatial coordinate (y)')
zlabel('Solution profile (P) \rightarrow')

% Loop for iterative solution
for p = 1:maxIter
    for i = 2:length(x)-1
         for j = 2:length(y)-1
              u(i,j,2) = (u(i+1,j,1) + u(i-1,j,1) + u(i,j+1,1) +u(i,j-1,1))/4;     
         end
     end
    
    % Updating the plot every Nth interation
    if mod(p,floor(sqrt(numPoints))) == 0
        s1.ZData = u(:,:,2);
        t1.String = {"Interation: " + p;...
            "Error: " + compose("%5.2e",err)};
        drawnow
    end

    %Stopping Condition???
    deltaU = abs(u(i,j,2) - u(i,j,1));
    err = deltaU;
   
    if err  < 1e-10
        break;
    else
        u(:,:,1) = u(:,:,2);
    end
end

