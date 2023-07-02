% Gauss’s Law
% Flux through the cube


clc;
clear;
% Note: if I decrease dx and dy step it takes a lot of time for PC to render
% the code
dx = 0.05;
dy = 0.05;
Area_Square = dx*dy;

x = 0:dx:1;
y = 0:dy:1;

xcharge = .2;
ycharge = .1;
zcharge = .5;

Z = zeros(length(x),length(y));
V = ones(length(x),length(y));

[X,Y] = meshgrid(x,y);

H = X(1:end-1,:) - dx/2;
I = Y(:,1:end-1) - dy/2;
J = H(:,2:end);
K = I(2:end,:);
L = zeros(length(J),length(K));
M = ones(length(J),length(K));
kq = (1.602 * (10^(-19))) / (4 * pi * (8.8541878128 *(10^(-12))));

normalvec = 0.1;
normalvec2 = 0.9;

% Output figure in full screen mode
fh = figure();
fh.WindowState = 'maximized';

surf(X,Y,Z) %surface 1
hold on
plot3(J,K,L, '.') %points 1
surf(X,Y,V,'FaceAlpha', 0.1) %surface 2
plot3(J,K,M, '.') %points 2
surf(Y,V,X,'FaceAlpha',0.1) %surface 3
plot3(K,M,J ,'.') %points 3
surf(V,X,Y,'FaceAlpha',0.1) %surface 4
plot3(M,J,K,'.') %points 4
surf(Z,X,Y) %surface 5
plot3(L,J,K,'.') %points 5
surf(Y,Z,X) %surface 6
plot3(K,L,J,'.') %points 6
plot3(xcharge,ycharge,zcharge,'k*')

for j = 1:width(J)
    for i = 1:length(J)

        plot3 ([J(i,j) J(i,j)],[K(i,j) K(i,j)],[normalvec L(i,j)],'b-')
        plot3 ([xcharge J(i,j)],[ycharge K(i,j)],[zcharge L(i,j)],'r-')
        
        plot3 ([K(i,j) K(i,j)],[normalvec L(i,j)],[J(i,j) J(i,j)],'b-')
        plot3 ([xcharge L(i,j)],[ycharge J(i,j)],[zcharge K(i,j)],'r-')
        
        plot3 ([normalvec L(i,j)],[J(i,j) J(i,j)],[K(i,j) K(i,j)],'b-')
        plot3([xcharge K(i,j)],[ycharge L(i,j)],[zcharge J(i,j)],'r-')
        
        plot3 ([J(i,j) J(i,j)],[K(i,j) K(i,j)],[normalvec2 M(i,j)],'b-')
        plot3 ([xcharge J(i,j)],[ycharge K(i,j)],[zcharge M(i,j)],'r-')

        plot3 ([normalvec2 M(i,j)],[J(i,j) J(i,j)],[K(i,j) K(i,j)],'b-')
        plot3 ([xcharge M(i,j)],[ycharge J(i,j)],[zcharge K(i,j)],'r-')

        plot3 ([K(i,j) K(i,j)],[normalvec2 M(i,j)],[J(i,j) J(i,j)],'b-')
        plot3 ([xcharge K(i,j)],[ycharge M(i,j)],[zcharge J(i,j)],'r-')  
    end
end

xlabel('X')
ylabel ('Y')
zlabel ('Z')
hold off

Q = sqrt(((0.-J).^2) + ((0.-K).^2)   + ((0-L).^2));

R1 = sqrt(((xcharge - J).^2) + ((ycharge - K).^2)   + ((zcharge - L).^2));
R2 = sqrt(((xcharge - J).^2) + ((ycharge-K).^2)   + ((zcharge-M).^2));
R3 = sqrt(((xcharge - K).^2) + ((ycharge-M).^2)   + ((zcharge-J).^2));
R4 = sqrt(((xcharge - M).^2) + ((ycharge-J).^2)   + ((zcharge-K).^2));
R5 = sqrt(((xcharge - L).^2) + ((ycharge-J).^2)   + ((zcharge-K).^2));
R6 = sqrt(((xcharge - K).^2) + ((ycharge-L).^2)   + ((zcharge-J).^2));

E1 = kq ./ (R1 .^ 2);
E2 = kq ./ (R2 .^ 2);
E3 = kq ./ (R3 .^ 2);
E4 = kq ./ (R4 .^ 2);
E5 = kq ./ (R5 .^ 2);
E6 = kq ./ (R6 .^ 2);

theta1 = acos(zcharge ./ R1);
theta2 = acos((1 - zcharge) ./ R2);
theta3 = acos((1 - ycharge) ./ R3);
theta4 = acos((1 - xcharge) ./ R4);
theta5 = acos(xcharge ./ R5);
theta6 = acos(ycharge ./ R6);

%theta = acos(Q ./ R); % dont really need rn 

flux1 = E1 .* Area_Square .* cos(theta1);
flux2 = E2 .* Area_Square .* cos(theta2);
flux3 = E3 .* Area_Square .* cos(theta3);
flux4 = E4 .* Area_Square .* cos(theta4);
flux5 = E5 .* Area_Square .* cos(theta5);
flux6 = E6 .* Area_Square .* cos(theta6);

Totalflux  = (sum(flux1,'all') + sum(flux2,'all')) + sum(flux3,'all') + sum(flux4,'all') + sum(flux5,'all') + sum(flux6,'all');

Gauss_law  = (1.602 * (10^(-19)))/(8.8541878128*(10^(-12)));


