clear % clear all stored variables
close all % close all open figures
clc % clear command window

%list variables for velocity field%
[x,y] = meshgrid(-3:.5:3,-3:.5:3); % x and y values for velocity field
u = 2*y; % u velocity function       
v = 1+(2*x); % v velocity function
xmarker = -0.5; %stagnation 'x' point
ymarker = 0; % stagnation 'y' point

%velocity plot%
figure % open a plot window
hold on; box on % hold current figure plot and place box around plot
quiver (x,y,u,v) % plot velocity vector field
plot(xmarker,ymarker,'r*') % stagnation point
axis([-3 3 -3 3]) % fix axis to desired range
title({'Aerodynamics Assignment 1';'Question 1 Velocity and Streamline plots'}) % plot title

%{
%stream function
x2 = [-3:.1:3]; y2 = [-3:.1:3];
psi = y2.^2 - x2.^2 - x2 -0.25; % stream function with constant determined from stagnation point
[C,h] = contour(x,y,psi,[1:1:6]); % plot stream lines from 1 to 6
[B,j] = contour(x,y,psi,[0:1:1],'-.'); % plot zero contour line dashed
clabel(C,h);clabel(B,j); % label contour values
%}

% streamline plots at finer resolution (0.14) to smooth out curves
x2 =[-3:.001:3];

% positive value curve set
psi0 = sqrt(x2.^2 +x2 + 0.25);
psi1 = sqrt(x2.^2 +x2 + 1.25);
psi2 = sqrt(x2.^2 +x2 + 2.25);
psi3 = sqrt(x2.^2 +x2 + 3.25);
psi4 = sqrt(x2.^2 +x2 + 4.25);
psi5 = sqrt(x2.^2 +x2 + 5.25);
psi6 = sqrt(x2.^2 +x2 + 6.25);
plot(x2,psi0,x2,psi1,x2,psi2,x2,psi3,x2,psi4,x2,psi5,x2,psi6)

% negative value curve set
psi00 = -sqrt(x2.^2 +x2 + 0.25);
psi11 = -sqrt(x2.^2 +x2 + 1.25);
psi22 = -sqrt(x2.^2 +x2 + 2.25);
psi33 = -sqrt(x2.^2 +x2 + 3.25);
psi44 = -sqrt(x2.^2 +x2 + 4.25);
psi55 = -sqrt(x2.^2 +x2 + 5.25);
psi66 = -sqrt(x2.^2 +x2 + 6.25);
plot(x2,psi00,x2,psi11,x2,psi22,x2,psi33,x2,psi44,x2,psi55,x2,psi66)

psi_half = sqrt(x2.^2 +x2 + 0.26);
psi0_half = -sqrt(x2.^2 +x2 + 0.26);
plot(x2,psi0_half,x2,psi_half);

%%
figure
[x,y,z] = meshgrid(-1.01:0.02:1.01,-1.01:0.02:1.01,-1.01:0.02:1.01);

U = 1; % x-wise free stream velocity
a = 0.5; % sphere radius
mu = -2*pi*a^3*U; % doublet strength required to create spherical streamlines

ux = U+mu.*(2.*x.^2-y.^2-z.^2)./(4*pi.*(x.^2+y.^2+z.^2).^(5/2));
uy = 3*mu.*x.*y./(4*pi.*(x.^2+y.^2+z.^2).^(5/2));
uz = 3*mu.*x.*z./(4*pi.*(x.^2+y.^2+z.^2).^(5/2));

[startx, starty, startz] = meshgrid(-1,-0.1:0.02:0.1,-0.1:0.02:0.1); 

H = streamline(x,y,z,ux,uy,uz,startx,starty,startz);
