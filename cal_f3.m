% ************************************************************************
% MATLAB Code: Force Calculation on a Charged Surface
% Author: Aslan Mohammadpour
% Affiliation: Imperial College London
% Email: Aslan.mohammadpour@gmail.com/
% A.mohammadpourshoorbakhlou@imperial.ac.uk
% Date: 10th August 2024
%
% Description:
% This MATLAB script computes the force vector acting on a charged surface
% by integrating the charge density distribution over the surface. The
% script calculates the normal vector using the cross product of the 
% partial derivatives of the position vector with respect to spherical 
% coordinates theta and phi. More detail of the theory is given in this
% paper: 
% https://www.sciencedirect.com/science/article/pii/S0304388615300334
%
% License:
% This file is provided "as is" without any express or implied warranties.
% Permission is granted to use, distribute, and modify this code, provided
% that proper citation is given to the original author. If you use this 
% code or its derivatives in your research, publications, or projects, 
% please cite the author:
%
% Aslan Mohammadpour
% Research associate in Computational Mechanics/Physics
% Imperial College London
%
% ************************************************************************

% Constants
epsilon_0 = 8.854187817e-12; % Vacuum permittivity (F/m)

% Define the parameters
phi_0 = 0.00001/3/(4*pi*epsilon_0*0.3);
r0 = 0.35; % Radius of the spherical conductor
a0 = phi_0*(4*pi*epsilon_0*r0); % Example value for a0

a1 = 0.25*(a0*r0); % Example value for a1

Psi = 1; % Dimensionless potential
k1 = a1 / (a0 * r0); % Dimensionless coefficient

syms Phi Theta

% Calculate the surface shape function xi(theta)
Xi = (1 + sqrt(1 + 4 * Psi * k1 * cos(Theta))) / (2 * Psi);



% Compute the partial derivatives of Xi with respect to theta
dXi_dTheta = - k1 * sin(Theta) / sqrt(1 + 4 * Psi * k1 * cos(Theta));

% Partial derivatives of r with respect to theta
dr_dTheta_x = dXi_dTheta * sin(Theta) * cos(Phi) + Xi * cos(Theta) * cos(Phi);
dr_dTheta_y = dXi_dTheta * sin(Theta) * sin(Phi) + Xi * cos(Theta) * sin(Phi);
dr_dTheta_z = dXi_dTheta * cos(Theta) - Xi * sin(Theta);

% Partial derivatives of r with respect to phi
dr_dPhi_x = -Xi * sin(Theta) * sin(Phi);
dr_dPhi_y = Xi * sin(Theta) * cos(Phi);
dr_dPhi_z = 0;

% Compute the normal vector using cross product
nx = dr_dTheta_y * dr_dPhi_z - dr_dTheta_z * dr_dPhi_y;
ny = dr_dTheta_z * dr_dPhi_x - dr_dTheta_x * dr_dPhi_z;
nz = dr_dTheta_x * dr_dPhi_y - dr_dTheta_y * dr_dPhi_x;


% Normalize the normal vector
norm_n = sqrt(nx^2 + ny^2 + nz^2);
nx = nx / norm_n;
ny = ny / norm_n;
nz = nz / norm_n;

% Compute the surface charge density sigma(theta)
sigma0 = a0 / (4 * pi * r0^2); % Surface charge density for a sphere
sigma_tilde = sqrt((1 / Xi^2 + 2 * k1 * cos(Theta) / Xi^3)^2 + (k1 * sin(Theta) / Xi^3)^2);
sigma = sigma0 *  sigma_tilde;



fx   = matlabFunction((sigma^2 / epsilon_0) * nx* (Xi*r0)^2* sin(Theta), 'Vars', [Theta, Phi]);
fy   = matlabFunction((sigma^2 / epsilon_0) * ny* (Xi*r0)^2* sin(Theta), 'Vars', [Theta, Phi]);
fz   = matlabFunction((sigma^2 / epsilon_0) * nz* (Xi*r0)^2* sin(Theta), 'Vars', [Theta, Phi]);

% Integral limits
theta_min = 0;
theta_max = pi;
phi_min = 0;
phi_max = 2 * pi;
% Compute the force components
Fx = integral2(fx, theta_min, theta_max, phi_min, phi_max);
Fy = integral2(fy, theta_min, theta_max, phi_min, phi_max);
Fz = integral2(fz, theta_min, theta_max, phi_min, phi_max);



% Display the results
fprintf('The force vector is: [Fx, Fy, Fz] = [%f, %f, %f] N\n', Fx, Fy, Fz);


%%
% Visualization of the surface and charge distribution
Theta_vals = linspace(theta_min, theta_max, 200);
Phi_vals = linspace(phi_min, phi_max, 100);
[Theta_mesh, Phi_mesh] = meshgrid(Theta_vals, Phi_vals);

% Convert the symbolic expressions to numeric values for plotting
Xi_vals = double(subs(Xi, Theta, Theta_mesh));
X_vals = Xi_vals .* sin(Theta_mesh) .* cos(Phi_mesh);
Y_vals = Xi_vals .* sin(Theta_mesh) .* sin(Phi_mesh);
Z_vals = Xi_vals .* cos(Theta_mesh);

% Compute the charge distribution for visualization
sigma_vals = double(subs(sigma, Theta, Theta_mesh));
%%

% Filter the points where y >= 0
positive_y_indices = Y_vals >= 0;

X_cut = X_vals(positive_y_indices);
Y_cut = Y_vals(positive_y_indices);
Z_cut = Z_vals(positive_y_indices);
sigma_cut = sigma_vals(positive_y_indices);

% Create a new figure for the cut plot
figure;
scatter3(X_cut, Y_cut, Z_cut, 10, sigma_cut, 'filled');
colorbar;
title('3D Surface Plot of Charged Surface with Contour of Charge Distribution (y >= 0)');
xlabel('X');
ylabel('Y');
zlabel('Z');
view(3);
axis equal;
shading interp;

% Enhance the plot with lighting
camlight;
lighting phong;


%%
% Plot the 3D surface with contour of charge distribution
figure;
surf(X_vals, Y_vals, Z_vals, sigma_vals, 'EdgeColor', 'none');
colorbar;
title('3D Surface Plot with Contour of Charge Distribution');
xlabel('X');
ylabel('Y');
zlabel('Z');
view(3);
axis equal;
shading interp;

% Add contour plot on the surface
hold on;
contour3(X_vals, Y_vals, Z_vals, sigma_vals, 20, 'LineColor', 'black');
hold off;

% Enhance the plot with lighting
% camlight;
% lighting phong;