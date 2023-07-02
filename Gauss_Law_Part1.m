% Gauss s Law
% A point charge is placed at the origin of an x,y,z coordinate system.

clc;
clear; 

dx = 0.05;
dy = 0.05;
AreaSquare = dx*dy;

x = 0:dx:0.5;
y = 0:dy:0.5;

Z = zeros(length(x),length(y))';

[X,Y] = meshgrid(x,y);

H = X(1:end-1,:) - dx/2;
I = Y(:,1:end-1) - dy/2;
J = H(:,2:end);
K = I(2:end,:);
L = zeros(length(J),width(J));

plot3(J,K,L, '.')
hold on 
surf(X,Y,Z)
plot3(0,0,0.5,'x') % Position of the Charge

for j = 1:width(J)
    for i = 1:length(J)
        plot3 ([J(i,j) J(i,j)],[K(i,j) K(i,j)],[0.1 L(i,j)],'b-') % Normal vectors at each partition
        plot3 ([0 J(i,j)],[0 K(i,j)],[0.5 L(i,j)],'r-')           % The electric field ad each partition
    end
end
hold off

R = sqrt(((0.-J).^2) + ((0.-K).^2)   + ((0.5-L).^2)); 
kq = (1.602 * (10^(-19))) / (4 * pi * (8.8541878128 *(10^(-12))));
E = kq ./ (R .^ 2);
Theta = acos(0.5 ./ R); 
Flux = E .* AreaSquare .* cos(Theta);
Totalflux  = sum(Flux,'all') * 24;
Gauss_law  = (1.602 * (10^(-19)))/(8.8541878128*(10^(-12)));

 