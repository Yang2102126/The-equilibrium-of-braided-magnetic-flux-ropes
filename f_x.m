function [f_x]=f_x(k,a,b)
% Computes the function f_x corresponding to Equation (A18) in
% Yang Zhang's "Magnetic Double Helix" paper.
%
% This function evaluates a series expression involving integrals of
% Bessel functions (I_n and K_n) and trigonometric terms, which arise
% in the analytic model for the magnetic field structure of braided flux ropes.
%
% Inputs:
%   k - dimensionless wavenumber parameter (e.g., k = r0 / λ)
%   a - major radius of the braid envelope
%   b - minor radius of the braid filament
%
% Output:
%   f_x - numerical evaluation of the analytical function defined in Eq. (A18)
%
% The series is computed term-by-term up to a maximum n, with convergence
% checked against a relative error threshold. If the integral evaluation
% encounters numerical instability at large n, an asymptotic expansion
% of the Bessel functions is used to ensure accurate evaluation.

% To avoid division by zero when k is zero, add a small epsilon.
epsilon=1e-6;
if k==0
    k=k+epsilon;
end
% Clear last warning messages, useful for checking numerical issues.
lastwarn('');
% Set maximum number of series terms and relative error threshold.
n_total=10^4;
relerror=10^(-4);
fun_1=@(r) (-integral(@(x) x.*acos((x.^2+(a^2-b^2))./(2*a*x)),a-b,r)+k^2*r.^2*integral(@(x) x.*acos((x.^2+(a^2-b^2))./(2*a*x)),r,a+b)).*sin(acos((r.^2+(a^2-b^2))./(2*a*r)));
f_x=16*a/(pi^2*b^4)*integral(fun_1,a-b,a+b,'ArrayValued',true);
for n=2:2:n_total
fun_2=@(u) -32*a/(pi^2*k^3*b^4)/n^4*...
    (((n*cos(acos((u.^2+(n*k)^2*(a^2-b^2))./(2*a*n*k*u))).*sin(n*acos((u.^2+(n*k)^2*(a^2-b^2))./(2*a*n*k*u)))-sin(acos((u.^2+(n*k)^2*(a^2-b^2))./(2*a*n*k*u))).*cos(n*acos((u.^2+(n*k)^2*(a^2-b^2))./(2*a*n*k*u))))/(n^2-1)).*(1+(u/n).^2).*...
    (besseli(n,u).*...
    integral(@(x) x.^2*(-1/2).*((besselk(n-1,x)+besselk(n+1,x)).*sin(n*acos((x.^2+(n*k)^2*(a^2-b^2))./(2*a*n*k*x)))),u,n*k*(a+b))+...
    besselk(n,u).*...
    integral(@(x) x.^2*(1/2).*((besseli(n-1,x)+besseli(n+1,x)).*sin(n*acos((x.^2+(n*k)^2*(a^2-b^2))./(2*a*n*k*x)))),n*k*(a-b),u))+...
    ((cos(acos((u.^2+(n*k)^2*(a^2-b^2))./(2*a*n*k*u))).*sin(n*acos((u.^2+(n*k)^2*(a^2-b^2))./(2*a*n*k*u)))-n*sin(acos((u.^2+(n*k)^2*(a^2-b^2))./(2*a*n*k*u))).*cos(n*acos((u.^2+(n*k)^2*(a^2-b^2))./(2*a*n*k*u))))/(n^2-1)).*(u/n).*...
    (1/2*(besseli(n-1,u)+besseli(n+1,u)).*...
    integral(@(x) x.^2*(-1/2).*((besselk(n-1,x)+besselk(n+1,x)).*sin(n*acos((x.^2+(n*k)^2*(a^2-b^2))./(2*a*n*k*x)))),u,n*k*(a+b))+...
    (-1/2)*(besselk(n-1,u)+besselk(n+1,u)).*...
    integral(@(x) x.^2*(1/2).*((besseli(n-1,x)+besseli(n+1,x)).*sin(n*acos((x.^2+(n*k)^2*(a^2-b^2))./(2*a*n*k*x)))),n*k*(a-b),u)));
    delta_f_x=integral(fun_2,n*k*(a-b),n*k*(a+b),'ArrayValued',true);
% When n becomes large, the modified Bessel function I_n(nkr) grows exponentially,
% while K_n(nkr) decays exponentially. Although their product remains finite,
% evaluating them separately in numerical integration can cause overflow,
% leading to numerical instability and warning messages.
% Therefore, when such issues are detected for large n, we switch to the
% asymptotic expansion of the Bessel functions to stabilize and accurately compute the integral.
[warnMsg, warnId] = lastwarn;    
if length(warnMsg)>0
    break
end
% Check convergence
if abs(delta_f_x/f_x)<relerror
break
end
f_x=f_x+delta_f_x;
n
end
n_break=n

% Use asymptotic expansion for large n terms beyond n_break
for n=n_break:2:n_total
fun_2_1=@(u) -32*a/(pi^2*k^3*b^4)/n^4*...
    (((n*cos(acos((u.^2+(n*k)^2*(a^2-b^2))./(2*a*n*k*u))).*sin(n*acos((u.^2+(n*k)^2*(a^2-b^2))./(2*a*n*k*u)))-sin(acos((u.^2+(n*k)^2*(a^2-b^2))./(2*a*n*k*u))).*cos(n*acos((u.^2+(n*k)^2*(a^2-b^2))./(2*a*n*k*u))))/(n^2-1)).*(1+(u/n).^2).*...
    (-integral(@(x) (x./2.*((1+(x/n).^2)./(1+(u/n).^2)).^(1/4).*...
    exp(n*((1+(u/n).^2).^(1/2)-(1+(x/n).^2).^(1/2)+log(u./x.*(1+(1+(x/n).^2).^(1/2))./(1+(1+(u/n).^2).^(1/2))))).*...
    (1+1/24/n*(3*(1+(u/n).^2).^(-1/2)-5*(1+(u/n).^2).^(-3/2))).*(1-1/24/n*(-9*(1+(x/n).^2).^(-1/2)+7*(1+(x/n).^2).^(-3/2))).*...
    sin(n*acos((x.^2+(n*k)^2*(a^2-b^2))./(2*a*n*k*x)))),u,n*k*(a+b))+...
    integral(@(x) (x./2.*((1+(x/n).^2)./(1+(u/n).^2)).^(1/4).*...
    exp(n*((1+(x/n).^2).^(1/2)-(1+(u/n).^2).^(1/2)+log(x./u.*(1+(1+(u/n).^2).^(1/2))./(1+(1+(x/n).^2).^(1/2))))).*...
    (1-1/24/n*(3*(1+(u/n).^2).^(-1/2)-5*(1+(u/n).^2).^(-3/2))).*(1+1/24/n*(-9*(1+(x/n).^2).^(-1/2)+7*(1+(x/n).^2).^(-3/2))).*...
    sin(n*acos((x.^2+(n*k)^2*(a^2-b^2))./(2*a*n*k*x)))),n*k*(a-b),u))+...
    ((cos(acos((u.^2+(n*k)^2*(a^2-b^2))./(2*a*n*k*u))).*sin(n*acos((u.^2+(n*k)^2*(a^2-b^2))./(2*a*n*k*u)))-n*sin(acos((u.^2+(n*k)^2*(a^2-b^2))./(2*a*n*k*u))).*cos(n*acos((u.^2+(n*k)^2*(a^2-b^2))./(2*a*n*k*u))))/(n^2-1)).*(u/n).*...
    (-(integral(@(x) (n*x./u/2.*((1+(u/n).^2).*(1+(x/n).^2)).^(1/4).*...
    exp(n*((1+(u/n).^2).^(1/2)-(1+(x/n).^2).^(1/2)+log(u./x.*(1+(1+(x/n).^2).^(1/2))./(1+(1+(u/n).^2).^(1/2))))).*...
    (1+1/24/n*(-9*(1+(u/n).^2).^(-1/2)+7*(1+(u/n).^2).^(-3/2))).*(1-1/24/n*(-9*(1+(x/n).^2).^(-1/2)+7*(1+(x/n).^2).^(-3/2))).*...
    sin(n*acos((x.^2+(n*k)^2*(a^2-b^2))./(2*a*n*k*x)))),u,n*k*(a+b))+...
    integral(@(x) (n*x./u/2.*((1+(u/n).^2).*(1+(x/n).^2)).^(1/4).*...
    exp(n*((1+(x/n).^2).^(1/2)-(1+(u/n).^2).^(1/2)+log(x./u.*(1+(1+(u/n).^2).^(1/2))./(1+(1+(x/n).^2).^(1/2))))).*...
    (1-1/24/n*(-9*(1+(u/n).^2).^(-1/2)+7*(1+(u/n).^2).^(-3/2))).*(1+1/24/n*(-9*(1+(x/n).^2).^(-1/2)+7*(1+(x/n).^2).^(-3/2))).*...
    sin(n*acos((x.^2+(n*k)^2*(a^2-b^2))./(2*a*n*k*x)))),n*k*(a-b),u))));
    % delta_f_x=integral(fun_2_1,n*k*(a-b),n*k*(a+b),'ArrayValued',true,'AbsTol',1e-8,'RelTol',1e-5);
    delta_f_x=integral(fun_2_1,n*k*(a-b),n*k*(a+b),'ArrayValued',true);
    n
    % Check convergence
    if abs(delta_f_x/f_x)<relerror
    break
    end
    f_x=f_x+delta_f_x;
end
end

