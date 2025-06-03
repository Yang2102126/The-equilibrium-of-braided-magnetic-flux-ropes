function [B_r,B_theta,B_z]=double_helix_mag(k,a,b,r,theta,z)
% Compute magnetic field components (B_r, B_theta, B_z) from a double-helix current structure
% based on the model described in Yang Zhang's "Magnetic Double Helix" paper.
%
% Inputs:
%   - a: Major radius of the helix
%   - b: Minor radius (cross-sectional radius)
%   - k: Helical wavenumber
%
% The magnetic field is computed using cylindrical harmonic expansions in three distinct radial regions:
%   1. Core:        r < a - b       (vacuum)
%   2. Inside:      a - b < r < a + b   (current-carrying shell)
%   3. Outside:     r > a + b       (vacuum)
%
% A dimensionless normalization is used for the magnetic field:
%     B̄ = B / (μ₀ * J_z * b² / 4a)
%
% The method evaluates Fourier-Bessel integrals and switches to asymptotic expansion for large harmonic numbers
% when convergence issues or numerical warnings are detected.

mu_0=1;      % Normalized vacuum permeability
J_z=4*a/b^2; % Normalized axial current density

% Numerical settings
n_total = 1e4;  % Maximum harmonic order
relerror = 1e-6;  % Relative error tolerance for convergence

% Clear last warning messages, useful for checking numerical issues.
lastwarn('');

% Region-specific evaluation of cylindrical components
if r < a - b
    [B_r, B_theta, B_z] = region_core(k, a, b, r, theta, z, mu_0, J_z, n_total, relerror);
elseif r > a + b
    [B_r, B_theta, B_z] = region_outside(k, a, b, r, theta, z, mu_0, J_z, n_total, relerror);
else
    [B_r, B_theta, B_z] = region_inside(k, a, b, r, theta, z, mu_0, J_z, n_total, relerror);
end

% Convert cylindrical (B_r, B_theta, B_z) to Cartesian (B_x, B_y, B_z)
B_x = B_r * cos(theta) - B_theta * sin(theta);
B_y = B_r * sin(theta) + B_theta * cos(theta);

end
%% region for r<a-b
function [B_r, B_theta, B_z] = region_core(k, a, b, r, theta, z, mu_0, J_z, n_total, relerror)
% Clear last warning messages, useful for checking numerical issues.
lastwarn('');
B_r=0;
B_theta=0;
B_z=2*mu_0*J_z*k/pi*integral(@(x) x.*acos((x.^2+(a^2-b^2))./(2*a*x)),a-b,a+b);
for n=2:2:n_total
fun_B_r=@(u) 4*mu_0*J_z/(pi*k)/n^3*sin(n*(theta-k*z)).*...
    (1/2*(besseli(n-1,u)+besseli(n+1,u)).*...
    integral(@(x) x.^2*(-1/2).*((besselk(n-1,x)+besselk(n+1,x)).*sin(n*acos((x.^2+(n*k)^2*(a^2-b^2))./(2*a*n*k*x)))),n*k*(a-b),n*k*(a+b)));
fun_B_theta=@(u) 4*mu_0*J_z./(pi*k*u)/n^2*cos(n*(theta-k*z)).*...
    (besseli(n,u).*...
    integral(@(x) x.^2*(-1/2).*((besselk(n-1,x)+besselk(n+1,x)).*sin(n*acos((x.^2+(n*k)^2*(a^2-b^2))./(2*a*n*k*x)))),n*k*(a-b),n*k*(a+b)));
fun_B_z=@(u) -4*mu_0*J_z./(pi*k)/n^3*cos(n*(theta-k*z)).*...
    (besseli(n,u).*...
    integral(@(x) x.^2*(-1/2).*((besselk(n-1,x)+besselk(n+1,x)).*sin(n*acos((x.^2+(n*k)^2*(a^2-b^2))./(2*a*n*k*x)))),n*k*(a-b),n*k*(a+b)));

delta_B_r=fun_B_r(n*k*r);
delta_B_theta=fun_B_theta(n*k*r);
delta_B_z=fun_B_z(n*k*r);

[warnMsg, warnId] = lastwarn;    
if length(warnMsg)>0
    n_break=n;
    disp('Warning detected. Switching to asymptotic expansion for integration.');
    break
end

B_r=B_r+delta_B_r;
B_theta=B_theta+delta_B_theta;
B_z=B_z+delta_B_z;

% Check convergence
epsilon=1e-8;
if (abs(delta_B_r / B_r) < relerror || abs(B_r) < epsilon) && ...
   (abs(delta_B_theta / B_theta) < relerror || abs(B_theta) < epsilon) && ...
   (abs(delta_B_z / B_z) < relerror || abs(B_z) < epsilon)
    break;  % Convergence achieved for all components, exit the loop
end

end


if length(warnMsg)>0
for n=n_break:2:n_total
fun_B_r_1=@(u) -4*mu_0*J_z/(pi*k)/n^2*sin(n*(theta-k*z)).*...
    (integral(@(x) (x./u/2.*((1+(u/n).^2).*(1+(x/n).^2)).^(1/4).*...
    exp(n*((1+(u/n).^2).^(1/2)-(1+(x/n).^2).^(1/2)+log(u./x.*(1+(1+(x/n).^2).^(1/2))./(1+(1+(u/n).^2).^(1/2))))).*...
    (1+1/24/n*(-9*(1+(u/n).^2).^(-1/2)+7*(1+(u/n).^2).^(-3/2))).*(1-1/24/n*(-9*(1+(x/n).^2).^(-1/2)+7*(1+(x/n).^2).^(-3/2))).*...
    sin(n*acos((x.^2+(n*k)^2*(a^2-b^2))./(2*a*n*k*x)))),n*k*(a-b),n*k*(a+b)));
fun_B_theta_1=@(u) 4*mu_0*J_z./(pi*k*u)/n^2*cos(n*(theta-k*z)).*...
    (-integral(@(x) (x./2.*((1+(x/n).^2)./(1+(u/n).^2)).^(1/4).*...
    exp(n*((1+(u/n).^2).^(1/2)-(1+(x/n).^2).^(1/2)+log(u./x.*(1+(1+(x/n).^2).^(1/2))./(1+(1+(u/n).^2).^(1/2))))).*...
    (1+1/24/n*(3*(1+(u/n).^2).^(-1/2)-5*(1+(u/n).^2).^(-3/2))).*(1-1/24/n*(-9*(1+(x/n).^2).^(-1/2)+7*(1+(x/n).^2).^(-3/2))).*...
    sin(n*acos((x.^2+(n*k)^2*(a^2-b^2))./(2*a*n*k*x)))),n*k*(a-b),n*k*(a+b)));
fun_B_z_1=@(u) -4*mu_0*J_z./(pi*k)/n^3*cos(n*(theta-k*z)).*...
    (-integral(@(x) (x./2.*((1+(x/n).^2)./(1+(u/n).^2)).^(1/4).*...
    exp(n*((1+(u/n).^2).^(1/2)-(1+(x/n).^2).^(1/2)+log(u./x.*(1+(1+(x/n).^2).^(1/2))./(1+(1+(u/n).^2).^(1/2))))).*...
    (1+1/24/n*(3*(1+(u/n).^2).^(-1/2)-5*(1+(u/n).^2).^(-3/2))).*(1-1/24/n*(-9*(1+(x/n).^2).^(-1/2)+7*(1+(x/n).^2).^(-3/2))).*...
    sin(n*acos((x.^2+(n*k)^2*(a^2-b^2))./(2*a*n*k*x)))),n*k*(a-b),n*k*(a+b)));

delta_B_r=fun_B_r_1(n*k*r);
delta_B_theta=fun_B_theta_1(n*k*r);
delta_B_z=fun_B_z_1(n*k*r);

B_r=B_r+delta_B_r;
B_theta=B_theta+delta_B_theta;
B_z=B_z+delta_B_z;

% Check convergence
epsilon=1e-8;
if (abs(delta_B_r / B_r) < relerror || abs(B_r) < epsilon) && ...
   (abs(delta_B_theta / B_theta) < relerror || abs(B_theta) < epsilon) && ...
   (abs(delta_B_z / B_z) < relerror || abs(B_z) < epsilon)
    break;  % Convergence achieved for all components, exit the loop
end

end
end
end

%% region for r>a+b
function [B_r, B_theta, B_z] = region_outside(k, a, b, r, theta, z, mu_0, J_z, n_total, relerror)
B_r=0;
B_theta=2*mu_0*J_z/pi/r*integral(@(x) x.*acos((x.^2+(a^2-b^2))./(2*a*x)),a-b,a+b);
B_z=0;
for n=2:2:n_total
fun_B_r=@(u) 4*mu_0*J_z/(pi*k)/n^3*sin(n*(theta-k*z)).*...
    ((-1/2)*(besselk(n-1,u)+besselk(n+1,u)).*...
    integral(@(x) x.^2*(1/2).*((besseli(n-1,x)+besseli(n+1,x)).*sin(n*acos((x.^2+(n*k)^2*(a^2-b^2))./(2*a*n*k*x)))),n*k*(a-b),n*k*(a+b)));
fun_B_theta=@(u) 4*mu_0*J_z./(pi*k*u)/n^2*cos(n*(theta-k*z)).*...
    (besselk(n,u).*...
    integral(@(x) x.^2*(1/2).*((besseli(n-1,x)+besseli(n+1,x)).*sin(n*acos((x.^2+(n*k)^2*(a^2-b^2))./(2*a*n*k*x)))),n*k*(a-b),n*k*(a+b)));
fun_B_z=@(u) -4*mu_0*J_z./(pi*k)/n^3*cos(n*(theta-k*z)).*...
    (besselk(n,u).*...
    integral(@(x) x.^2*(1/2).*((besseli(n-1,x)+besseli(n+1,x)).*sin(n*acos((x.^2+(n*k)^2*(a^2-b^2))./(2*a*n*k*x)))),n*k*(a-b),n*k*(a+b)));

delta_B_r=fun_B_r(n*k*r);
delta_B_theta=fun_B_theta(n*k*r);
delta_B_z=fun_B_z(n*k*r);

[warnMsg, warnId] = lastwarn;    
if length(warnMsg)>0
    n_break=n;
    disp('Warning detected. Switching to asymptotic expansion for integration.');
    break
end

B_r=B_r+delta_B_r;
B_theta=B_theta+delta_B_theta;
B_z=B_z+delta_B_z;

% Check convergence
epsilon=1e-8;
if (abs(delta_B_r / B_r) < relerror || abs(B_r) < epsilon) && ...
   (abs(delta_B_theta / B_theta) < relerror || abs(B_theta) < epsilon) && ...
   (abs(delta_B_z / B_z) < relerror || abs(B_z) < epsilon)
    break;  % Convergence achieved for all components, exit the loop
end

end

if length(warnMsg)>0
for n=n_break:2:n_total
fun_B_r_1=@(u) -4*mu_0*J_z/(pi*k)/n^2*sin(n*(theta-k*z)).*...
    (integral(@(x) (x./u/2.*((1+(u/n).^2).*(1+(x/n).^2)).^(1/4).*...
    exp(n*((1+(x/n).^2).^(1/2)-(1+(u/n).^2).^(1/2)+log(x./u.*(1+(1+(u/n).^2).^(1/2))./(1+(1+(x/n).^2).^(1/2))))).*...
    (1-1/24/n*(-9*(1+(u/n).^2).^(1/2)+7*(1+(u/n).^2).^(-3/2))).*(1+1/24/n*(-9*(1+(x/n).^2).^(1/2)+7*(1+(x/n).^2).^(-3/2))).*...
    sin(n*acos((x.^2+(n*k)^2*(a^2-b^2))./(2*a*n*k*x)))),n*k*(a-b),n*k*(a+b)));
fun_B_theta_1=@(u) 4*mu_0*J_z./(pi*k*u)/n^2*cos(n*(theta-k*z)).*...
    (integral(@(x) (x./2.*((1+(x/n).^2)./(1+(u/n).^2)).^(1/4).*...
    exp(n*((1+(x/n).^2).^(1/2)-(1+(u/n).^2).^(1/2)+log(x./u.*(1+(1+(u/n).^2).^(1/2))./(1+(1+(x/n).^2).^(1/2))))).*...
    (1-1/24/n*(3*(1+(u/n).^2).^(-1/2)-5*(1+(u/n).^2).^(-3/2))).*(1+1/24/n*(-9*(1+(x/n).^2).^(-12)+7*(1+(x/n).^2).^(-3/2))).*...
    sin(n*acos((x.^2+(n*k)^2*(a^2-b^2))./(2*a*n*k*x)))),n*k*(a-b),n*k*(a+b)));
fun_B_z_1=@(u) -4*mu_0*J_z./(pi*k)/n^3*cos(n*(theta-k*z)).*...
    (integral(@(x) (x./2.*((1+(x/n).^2)./(1+(u/n).^2)).^(1/4).*...
    exp(n*((1+(x/n).^2).^(1/2)-(1+(u/n).^2).^(1/2)+log(x./u.*(1+(1+(u/n).^2).^(1/2))./(1+(1+(x/n).^2).^(1/2))))).*...
    (1-1/24/n*(3*(1+(u/n).^2).^(-1/2)-5*(1+(u/n).^2).^(-3/2))).*(1+1/24/n*(-9*(1+(x/n).^2).^(-12)+7*(1+(x/n).^2).^(-3/2))).*...
    sin(n*acos((x.^2+(n*k)^2*(a^2-b^2))./(2*a*n*k*x)))),n*k*(a-b),n*k*(a+b)));
delta_B_r=fun_B_r_1(n*k*r);
delta_B_theta=fun_B_theta_1(n*k*r);
delta_B_z=fun_B_z_1(n*k*r);


B_r=B_r+delta_B_r;
B_theta=B_theta+delta_B_theta;
B_z=B_z+delta_B_z;

% Check convergence
epsilon=1e-8;
if (abs(delta_B_r / B_r) < relerror || abs(B_r) < epsilon) && ...
   (abs(delta_B_theta / B_theta) < relerror || abs(B_theta) < epsilon) && ...
   (abs(delta_B_z / B_z) < relerror || abs(B_z) < epsilon)
    break;  % Convergence achieved for all components, exit the loop
end

end
end
end

%% region for a-b<r<a+b
function [B_r, B_theta, B_z] = region_inside(k, a, b, r, theta, z, mu_0, J_z, n_total, relerror)
B_r=0;
B_theta=2*mu_0*J_z/pi/r*integral(@(x) x.*acos((x.^2+(a^2-b^2))./(2*a*x)),a-b,r);
B_z=2*mu_0*J_z*k/pi*integral(@(x) x.*acos((x.^2+(a^2-b^2))./(2*a*x)),r,a+b);
for n=2:2:n_total
fun_B_r=@(u) 4*mu_0*J_z/(pi*k)/n^3*sin(n*(theta-k*z)).*...
    (1/2*(besseli(n-1,u)+besseli(n+1,u)).*...
    integral(@(x) x.^2*(-1/2).*((besselk(n-1,x)+besselk(n+1,x)).*sin(n*acos((x.^2+(n*k)^2*(a^2-b^2))./(2*a*n*k*x)))),u,n*k*(a+b))+...
    (-1/2)*(besselk(n-1,u)+besselk(n+1,u)).*...
    integral(@(x) x.^2*(1/2).*((besseli(n-1,x)+besseli(n+1,x)).*sin(n*acos((x.^2+(n*k)^2*(a^2-b^2))./(2*a*n*k*x)))),n*k*(a-b),u));
fun_B_theta=@(u) 4*mu_0*J_z./(pi*k*u)/n^2*cos(n*(theta-k*z)).*...
    (besseli(n,u).*...
    integral(@(x) x.^2*(-1/2).*((besselk(n-1,x)+besselk(n+1,x)).*sin(n*acos((x.^2+(n*k)^2*(a^2-b^2))./(2*a*n*k*x)))),u,n*k*(a+b))+...
    besselk(n,u).*...
    integral(@(x) x.^2*(1/2).*((besseli(n-1,x)+besseli(n+1,x)).*sin(n*acos((x.^2+(n*k)^2*(a^2-b^2))./(2*a*n*k*x)))),n*k*(a-b),u));
fun_B_z=@(u) -4*mu_0*J_z./(pi*k)/n^3*cos(n*(theta-k*z)).*...
    (besseli(n,u).*...
    integral(@(x) x.^2*(-1/2).*((besselk(n-1,x)+besselk(n+1,x)).*sin(n*acos((x.^2+(n*k)^2*(a^2-b^2))./(2*a*n*k*x)))),u,n*k*(a+b))+...
    besselk(n,u).*...
    integral(@(x) x.^2*(1/2).*((besseli(n-1,x)+besseli(n+1,x)).*sin(n*acos((x.^2+(n*k)^2*(a^2-b^2))./(2*a*n*k*x)))),n*k*(a-b),u));
delta_B_r=fun_B_r(n*k*r);
delta_B_theta=fun_B_theta(n*k*r);
delta_B_z=fun_B_z(n*k*r);

[warnMsg, warnId] = lastwarn;    
if length(warnMsg)>0
    n_break=n;
    disp('Warning detected. Switching to asymptotic expansion for integration.');
    break
end

B_r=B_r+delta_B_r;
B_theta=B_theta+delta_B_theta;
B_z=B_z+delta_B_z;

% Check convergence
epsilon=1e-8;
if (abs(delta_B_r / B_r) < relerror || abs(B_r) < epsilon) && ...
   (abs(delta_B_theta / B_theta) < relerror || abs(B_theta) < epsilon) && ...
   (abs(delta_B_z / B_z) < relerror || abs(B_z) < epsilon)
    break;  % Convergence achieved for all components, exit the loop
end


end

if length(warnMsg)>0
for n=n_break:2:n_total
fun_B_r_1=@(u) -4*mu_0*J_z/(pi*k)/n^2*sin(n*(theta-k*z)).*...
    (integral(@(x) (x./u/2.*((1+(u/n).^2).*(1+(x/n).^2)).^(1/4).*...
    exp(n*((1+(u/n).^2).^(1/2)-(1+(x/n).^2).^(1/2)+log(u./x.*(1+(1+(x/n).^2).^(1/2))./(1+(1+(u/n).^2).^(1/2))))).*...
    (1+1/24/n*(-9*(1+(u/n).^2).^(-1/2)+7*(1+(u/n).^2).^(-3/2))).*(1-1/24/n*(-9*(1+(x/n).^2).^(-12)+7*(1+(x/n).^2).^(-3/2))).*...
    sin(n*acos((x.^2+(n*k)^2*(a^2-b^2))./(2*a*n*k*x)))),u,n*k*(a+b))+...
    integral(@(x) (x./u/2.*((1+(u/n).^2).*(1+(x/n).^2)).^(1/4).*...
    exp(n*((1+(x/n).^2).^(1/2)-(1+(u/n).^2).^(1/2)+log(x./u.*(1+(1+(u/n).^2).^(1/2))./(1+(1+(x/n).^2).^(1/2))))).*...
    (1-1/24/n*(-9*(1+(u/n).^2).^(-1/2)+7*(1+(u/n).^2).^(-3/2))).*(1+1/24/n*(-9*(1+(x/n).^2).^(-1/2)+7*(1+(x/n).^2).^(-3/2))).*...
    sin(n*acos((x.^2+(n*k)^2*(a^2-b^2))./(2*a*n*k*x)))),n*k*(a-b),u));
fun_B_theta_1=@(u) 4*mu_0*J_z./(pi*k*u)/n^2*cos(n*(theta-k*z)).*...
    (-integral(@(x) (x./2.*((1+(x/n).^2)./(1+(u/n).^2)).^(1/4).*...
    exp(n*((1+(u/n).^2).^(1/2)-(1+(x/n).^2).^(1/2)+log(u./x.*(1+(1+(x/n).^2).^(1/2))./(1+(1+(u/n).^2).^(1/2))))).*...
    (1+1/24/n*(3*(1+(u/n).^2).^(-1/2)-5*(1+(u/n).^2).^(-3/2))).*(1-1/24/n*(-9*(1+(x/n).^2).^(-1/2)+7*(1+(x/n).^2).^(-3/2))).*...
    sin(n*acos((x.^2+(n*k)^2*(a^2-b^2))./(2*a*n*k*x)))),u,n*k*(a+b))+...
    integral(@(x) (x./2.*((1+(x/n).^2)./(1+(u/n).^2)).^(1/4).*...
    exp(n*((1+(x/n).^2).^(1/2)-(1+(u/n).^2).^(1/2)+log(x./u.*(1+(1+(u/n).^2).^(1/2))./(1+(1+(x/n).^2).^(1/2))))).*...
    (1-1/24/n*(3*(1+(u/n).^2).^(-1/2)-5*(1+(u/n).^2).^(-3/2))).*(1+1/24/n*(-9*(1+(x/n).^2).^(-1/2)+7*(1+(x/n).^2).^(-3/2))).*...
    sin(n*acos((x.^2+(n*k)^2*(a^2-b^2))./(2*a*n*k*x)))),n*k*(a-b),u));
fun_B_z_1=@(u) -4*mu_0*J_z./(pi*k)/n^3*cos(n*(theta-k*z)).*...
    (-integral(@(x) (x./2.*((1+(x/n).^2)./(1+(u/n).^2)).^(1/4).*...
    exp(n*((1+(u/n).^2).^(1/2)-(1+(x/n).^2).^(1/2)+log(u./x.*(1+(1+(x/n).^2).^(1/2))./(1+(1+(u/n).^2).^(1/2))))).*...
    (1+1/24/n*(3*(1+(u/n).^2).^(-1/2)-5*(1+(u/n).^2).^(-3/2))).*(1-1/24/n*(-9*(1+(x/n).^2).^(-1/2)+7*(1+(x/n).^2).^(-3/2))).*...
    sin(n*acos((x.^2+(n*k)^2*(a^2-b^2))./(2*a*n*k*x)))),u,n*k*(a+b))+...
    integral(@(x) (x./2.*((1+(x/n).^2)./(1+(u/n).^2)).^(1/4).*...
    exp(n*((1+(x/n).^2).^(1/2)-(1+(u/n).^2).^(1/2)+log(x./u.*(1+(1+(u/n).^2).^(1/2))./(1+(1+(x/n).^2).^(1/2))))).*...
    (1-1/24/n*(3*(1+(u/n).^2).^(-1/2)-5*(1+(u/n).^2).^(-3/2))).*(1+1/24/n*(-9*(1+(x/n).^2).^(-1/2)+7*(1+(x/n).^2).^(-3/2))).*...
    sin(n*acos((x.^2+(n*k)^2*(a^2-b^2))./(2*a*n*k*x)))),n*k*(a-b),u));
delta_B_r=fun_B_r_1(n*k*r);
delta_B_theta=fun_B_theta_1(n*k*r);
delta_B_z=fun_B_z_1(n*k*r);


B_r=B_r+delta_B_r;
B_theta=B_theta+delta_B_theta;
B_z=B_z+delta_B_z;

% Check convergence
epsilon=1e-8;
if (abs(delta_B_r / B_r) < relerror || abs(B_r) < epsilon) && ...
   (abs(delta_B_theta / B_theta) < relerror || abs(B_theta) < epsilon) && ...
   (abs(delta_B_z / B_z) < relerror || abs(B_z) < epsilon)
    break;  % Convergence achieved for all components, exit the loop
end

end
end

end





