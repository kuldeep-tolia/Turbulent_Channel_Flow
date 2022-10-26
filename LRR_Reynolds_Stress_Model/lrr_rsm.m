% DNS data at Re_delta=7890, Re_tau=395 (Moin, Kim & Mansour, PoF, 1999).
% All quantites are normalized by u_tau and nu unless stated otherwise.
% Delta denotes the channel half-width.

format long
clc; close all; clear all

% Read DNS data [half-channel is given (till centerline)]

load y_dns.dat
load u_dns.dat
load u2_dns.dat
load v2_dns.dat
load w2_dns.dat
load uv_dns.dat
load dns_data.dat

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nu = 1 / 395;        % Fluid viscosity
rho = 1;             % Fluid density
kappa = 0.41;        % Von Karmann constant 
dPdx = -1.00;        % Pressure gradient
delta = 1.0;         % Channel half-width
E = 9.0;
ny = 100;            % number of grid points after first node placement 

% RSM constants

c_mu = 0.09;
c1 = 1.8;
c2 = 0.6;
c1prime = 0.5;
c2prime = 0.3;
c1_eps = 1.44;
c2_eps = 1.92;
sigma_k = 1.0;
sigma_eps = 1.3; 

% Residual error limit

residue_limit = 10^(-4);

% Under relaxation factor

urf = 0.6;

% Initialising by running k-omega model  

[y_node, y_face, n, U, k, epsilon, uu, vv, ww, uv, ustar] = k_epsilon(nu, rho, dPdx, delta, ny);

% Additional variables
  
numerator = 0.0;
denominator = 0.0;
residue = 1.0;
U_old(:) = U(:);
k_old(:) = k(:);
uu_old(:) = uu(:);
vv_old(:) = vv(:);
ww_old(:) = ww(:);
uv_old(:) = uv(:);
epsilon_old(:) = epsilon(:);
counter = 1;

% Initialise more variables

dUdy = 0.1 * ones(n, 1);
duvdy = 0.1 * ones(n, 1);
P11 = 0.1 * ones(n, 1);
nu_t = rand(n, 1);
nu_t_face = 0.1 * ones(n-1, 1);
f = 0.1 * ones(n, 1);
tau_w = ustar^2;
 
while (residue > residue_limit)
   
   % Compute nu_t
 
   for i = 2:n
 
       nu_t(i) = c_mu * k(i) * k(i) / epsilon(i);
 
   end
 
   nu_t(1) = 0.0;

   % Compute nu_t at faces using harmonic mean approach

   nu_t_face(1) = nu_t(1);
   nu_t_face(n-1) = nu_t(n);

   for i = 2:n-2

       nu_t_face(i) = (2.0 * (y_node(i+1) - y_node(i))) / (((y_face(i) - y_face(i-1)) / nu_t(i)) + ((y_face(i+1) - y_face(i)) / nu_t(i+1)));    

   end

   % Compute U
  
   for i = 2:n-1         % Gauss-Seidel loop
    
     if (i == 2)         % Bottom wall cell treatment
      
       aN = nu / (y_node(i+1) - y_node(i));
       aS = nu / (y_node(i) - y_node(i-1));
       Sc = -ustar^2 - duvdy(i);
       b = Sc * (y_face(i) - y_face(i-1));
       aP = aN + aS;
       U(i) = ((aN * U(i+1)) + b) / aP;
      
     elseif (i == (n-1)) % Top wall cell treatment
      
       aN = 0.0;
       aS = nu / (y_node(i) - y_node(i-1));
       Sc = -dPdx - duvdy(i);
       b = Sc * (y_face(i) - y_face(i-1));
       aP = aS;
       U(i) = ((aS * U(i-1)) + b) / aP;
      
     else
      
       aN = nu / (y_node(i+1) - y_node(i));
       aS = nu / (y_node(i) - y_node(i-1));
       Sc = -dPdx - duvdy(i);
       b = Sc * (y_face(i) - y_face(i-1));
       aP = aN + aS;
       U(i) = ((aN * U(i+1)) + (aS * U(i-1)) + b) / aP;
       
     end
    
   end
  
   U(1) = 0.0;        % velocity wall bc
   U(n) = U(n-1);     % velocity symmetry bc

   % Compute gradients 
  
   dUdy(1) = (-U(3) + (4 * U(2)) - (3 * U(1))) / (y_node(3) - y_node(1));          % 3 point forward difference scheme
   dUdy(n) = 0.0;                                                                  % velocity symmetry condition
  
   duvdy(1) = (-uv(3) + (4 * uv(2)) - (3 * uv(1))) / (y_node(3) - y_node(1));      % 3 point forward difference scheme
   duvdy(n) = ((3 * uv(n)) - (4 * uv(n-1)) + uv(n-2)) / (y_node(n) - y_node(n-2)); % 3 point backward difference scheme
  
   for i = 2:n-1
    
     dUdy(i) = (U(i+1) - U(i-1)) / (y_node(i+1) - y_node(i-1));    % central difference scheme
    
     duvdy(i) = (uv(i+1) - uv(i-1)) / (y_node(i+1) - y_node(i-1)); % central difference scheme
    
   end
  
   % Compute ustar and tau_w
  
   ustar = (kappa * U(2)) / log(E * ustar * y_node(2) / nu);
   tau_w = ustar^2;
  
   % Compute k 
  
   k = 0.50 * (uu + vv + ww);
  
   % Compute P11
  
   for i = 2:n
    
     P11(i) = -2.0 * uv(i) * dUdy(i); 
 
   end

   P11(1) = 0.0;
  
   % Compute Pk
  
   Pk(:) = 0.50 * P11(:);
  
   % Compute f
  
   for i = 2:n
    
     f(i) = (k(i) ^ 1.5) / (2.55 * y_node(i) * epsilon(i));
    
   end
   
   % Compute uu

   uu(1) = 0.0;               % fluctuations wall bc
   uu(2) = 3.67 * ustar^2;    % wall function condition
      
   for i = 3:n-1
     
     if (i == (n-1))
       
       aN = 0;
       aS = (nu + (nu_t_face(i-1) / sigma_k)) / (y_node(i) - y_node(i-1));
       Sc = (2 * uv(i) * dUdy(i) * (c2 - 1)) + (2 * epsilon(i) * (c1 - 1) / 3) + (2 * c2 * Pk(i) * (1 + c2prime * f(i)) / 3) + (c1prime * epsilon(i) * vv(i) * f(i) / k(i));
       Sp = -c1 * epsilon(i) / k(i);
       b = Sc * (y_face(i) - y_face(i-1));
       aP = aS - (Sp * (y_face(i) - y_face(i-1)));
       uu(i) = ((aS * uu(i-1)) + b) / aP;
       
     else
       
       aN = (nu + (nu_t_face(i) / sigma_k)) / (y_node(i+1) - y_node(i));
       aS = (nu + (nu_t_face(i-1) / sigma_k)) / (y_node(i) - y_node(i-1));
       Sc = (2 * uv(i) * dUdy(i) * (c2 - 1)) + (2 * epsilon(i) * (c1 - 1) / 3) + (2 * c2 * Pk(i) * (1 + c2prime * f(i)) / 3) + (c1prime * epsilon(i) * vv(i) * f(i) / k(i));
       Sp = -c1 * epsilon(i) / k(i);
       b = Sc * (y_face(i) - y_face(i-1));
       aP = aN + aS - (Sp * (y_face(i) - y_face(i-1)));
       uu(i) = ((aN * uu(i+1)) + (aS * uu(i-1)) + b) / aP;
       
     end
   end
   
   uu(n) = uu(n-1);           % symmetry bc
   
   
   % Compute vv

   vv(1) = 0.0;               % fluctuations wall bc
   vv(2) = 0.83 * ustar^2;    % wall function condition
   
   for i = 3:n-1
     
     if (i == (n-1))
       
       aN = 0;
       aS = (nu + (nu_t_face(i-1) / sigma_k)) / (y_node(i) - y_node(i-1));
       Sc = (2 * epsilon(i) * (c1 - 1) / 3) + (2 * c2 * Pk(i) * (1 - 2 * c2prime * f(i)) / 3);
       Sp = -((c1 + 2 * c1prime * f(i)) * epsilon(i) / k(i));
       b = Sc * (y_face(i) - y_face(i-1));
       aP = aS - (Sp * (y_face(i) - y_face(i-1)));
       vv(i) = ((aS * vv(i-1)) + b) / aP;
       
     else
       
       aN = (nu + (nu_t_face(i) / sigma_k)) / (y_node(i+1) - y_node(i));
       aS = (nu + (nu_t_face(i-1) / sigma_k)) / (y_node(i) - y_node(i-1));
       Sc = (2 * epsilon(i) * (c1 - 1) / 3) + (2 * c2 * Pk(i) * (1 - 2 * c2prime * f(i)) / 3);
       Sp = -((c1 + 2 * c1prime * f(i)) * epsilon(i) / k(i));
       b = Sc * (y_face(i) - y_face(i-1));
       aP = aN + aS - (Sp * (y_face(i) - y_face(i-1)));
       vv(i) = ((aN * vv(i+1)) + (aS * vv(i-1)) + b) / aP;
       
     end
   
   end
   
   vv(n) = vv(n-1);           % symmetry bc
   
   % Compute ww

   ww(1) = 0.0;               % fluctuations wall bc
   ww(2) = 2.17 * ustar^2;    % wall function condition
   
   for i = 3:n-1
     
     if (i == (n-1))
       
       aN = 0.0;
       aS = (nu + (nu_t_face(i-1) / sigma_k)) / (y_node(i) - y_node(i-1));
       Sc = (2 * epsilon(i) * (c1 - 1) / 3) + (2 * c2 * Pk(i) * (1 + c2prime * f(i)) / 3) + (c1prime * epsilon(i) * vv(i) * f(i) / k(i));
       Sp = -c1 * epsilon(i) / k(i);
       b = Sc * (y_face(i) - y_face(i-1));
       aP = aS - (Sp * (y_face(i) - y_face(i-1)));
       ww(i) = ((aS * ww(i-1)) + b) / aP;
       
     else
       
       aN = (nu + (nu_t_face(i) / sigma_k)) / (y_node(i+1) - y_node(i));
       aS = (nu + (nu_t_face(i-1) / sigma_k)) / (y_node(i) - y_node(i-1));
       Sc = (2 * epsilon(i) * (c1 - 1) / 3) + (2 * c2 * Pk(i) * (1 + c2prime * f(i)) / 3) + (c1prime * epsilon(i) * vv(i) * f(i) / k(i));
       Sp = -c1 * epsilon(i) / k(i);
       b = Sc * (y_face(i) - y_face(i-1));
       aP = aN + aS - (Sp * (y_face(i) - y_face(i-1)));
       ww(i) = ((aN * ww(i+1)) + (aS * ww(i-1)) + b) / aP;
       
     end
   
   end
   
   ww(n) = ww(n-1);         % symmetry bc
   
   % Compute uv

   uv(1) = 0.0;               % fluctuations wall bc
   uv(n) = 0.0;               % zero shear stress at symmetry
   uv(2) = -1 * ustar^2;          % wall function condition
   
   for i = 3:n-1
     
     if (i == (n-1))
       
       aN = (nu + (nu_t_face(i) / sigma_k)) / (y_node(i+1) - y_node(i));
       aS = (nu + (nu_t_face(i-1) / sigma_k)) / (y_node(i) - y_node(i-1));
       Sc = (vv(i) * dUdy(i) * (c2 - 1 - 1.5 * c2 * c2prime * f(i))) - (uv(i) * epsilon(i) * (c1 + 1.5 * c1prime *f(i)) / k(i));
       b = Sc * (y_face(i) - y_face(i-1));
       aP = aN + aS;
       uv(i) = ((aS * uv(i-1)) + b) / aP;
       
     else
       
       aN = (nu + (nu_t_face(i) / sigma_k)) / (y_node(i+1) - y_node(i));
       aS = (nu + (nu_t_face(i-1) / sigma_k)) / (y_node(i) - y_node(i-1));
       Sc = (vv(i) * dUdy(i) * (c2 - 1 - 1.5 * c2 * c2prime * f(i))) - (uv(i) * epsilon(i) * (c1 + 1.5 * c1prime *f(i)) / k(i));
       b = Sc * (y_face(i) - y_face(i-1));
       aP = aN + aS;
       uv(i) = ((aN * uv(i+1)) + (aS * uv(i-1)) + b) / aP;
       
     end
     
   end
   
   % Compute epsilon

   epsilon(2) = ustar^3 / (kappa * y_node(2));    % epsilon wall function condition
      
   for i = 3:n-1
     
     if (i == (n-1))
       
       aN = 0.0;
       aS = (nu + (nu_t_face(i-1) / sigma_eps)) / (y_node(i) - y_node(i-1));
       Sc = c1_eps * epsilon(i) * Pk(i) / k(i);
       Sp = -c2_eps * epsilon(i) / k(i);
       b = Sc * (y_face(i) - y_face(i-1));
       aP = aS - (Sp * (y_face(i) - y_face(i-1)));
       epsilon(i) = ((aS * epsilon(i-1)) + b) / aP;
       
     else 
       
       aN = (nu + (nu_t_face(i) / sigma_eps)) / (y_node(i+1) - y_node(i));
       aS = (nu + (nu_t_face(i-1) / sigma_eps)) / (y_node(i) - y_node(i-1));
       Sc = c1_eps * epsilon(i) * Pk(i) / k(i);
       Sp = -c2_eps * epsilon(i) / k(i);
       b = Sc * (y_face(i) - y_face(i-1));
       aP = aN + aS - (Sp * (y_face(i) - y_face(i-1)));
       epsilon(i) = ((aN * epsilon(i+1)) + (aS * epsilon(i-1)) + b) / aP;
       
     end
     
   end
   
   epsilon(n) = epsilon(n-1);                      % epsilon symmetry bc

   % Under relax variables

   U(:) = (urf * U(:)) + ((1 - urf) * U_old(:));
   k(:) = (urf * k(:)) + ((1 - urf) * k_old(:));
   uu(:) = (urf * uu(:)) + ((1 - urf) * uu_old(:));
   vv(:) = (urf * vv(:)) + ((1 - urf) * vv_old(:));
   ww(:) = (urf * ww(:)) + ((1 - urf) * ww_old(:));
   uv(:) = (urf * uv(:)) + ((1 - urf) * uv_old(:));
   epsilon(:) = (urf * epsilon(:)) + ((1 - urf) * epsilon_old(:));

   numerator = 0.0;
   denominator = 0.0;

   for i =1:n

       numerator = numerator + (U(i) - U_old(i))^2;
       denominator = denominator + U(i)^2;

   end

   residue = sqrt(numerator / denominator);

   U_old(:) = U(:);
   k_old(:) = k(:);
   uu_old(:) = uu(:);
   vv_old(:) = vv(:);
   ww_old(:) = ww(:);
   uv_old(:) = uv(:);
   epsilon_old(:) = epsilon(:);

   counter = counter + 1;

   t = mod(counter, 1000);

   if (t == 0)

       fprintf('iter = %i\n', counter)
       fprintf('U = %.7f\n', U(2))
       fprintf('residual = %.7f\n', residue)
       fprintf('nut = %.7f\n', nu_t(2))
       fprintf('Pk = %.7f\n', Pk(2))
       fprintf('k = %.7f\n', k(2))
       fprintf('uu = %.7f\n', uu(2))
       fprintf('vv = %.7f\n', vv(2));
       fprintf('ww = %.7f\n', ww(2));
       fprintf('uv = %.7f\n', uv(2));
       fprintf('epsilon = %.7f\n', epsilon(2));
       fprintf('\n')

   end
   
end

% Plotting and Post-Processing

eps_dns = dns_data(:,2) / nu;         % eps is normalized by ustar^4/nu
y_plus = y_dns * 1 /nu;
y_plus(y_plus < 30) = [];
U_plus_log_law = (log(y_plus) / kappa) + 5.2;

fprintf('Wall Shear Stress = %.7f\n\n', tau_w);

figure(1)
plot(y_dns, u_dns, 'bo', y_node(2:n), U(2:n), 'k-', 'LineWidth', 2);
xlabel('y (m)'); ylabel('U (m / s)'); title('U-Velocity');
legend('DNS data', 'RSM data','Location','northeast'); legend boxoff;
ax = gca;
ax.FontSize = 24;

figure(2)
plot(y_dns, 0.5 * (u2_dns + v2_dns + w2_dns), 'bo', y_node(2:n), k(2:n), 'k-', 'LineWidth', 2);
xlabel('y (m)'); ylabel('k (m^2 / s^2)'); title('Turbulent Kinetic Energy');
legend('DNS data', 'RSM data', 'Location', 'northeast'); legend boxoff;
ax = gca;
ax.FontSize = 24;

figure(3)
plot(y_dns, u2_dns, 'bo', y_node(2:n), uu(2:n), 'k-', 'LineWidth', 2);
xlabel('y (m)'); ylabel('<uu> (m^2 / s^2)'); title('Reynolds Normal Stress <uu>');
legend('DNS data', 'RSM data', 'Location', 'northeast'); legend boxoff;
ax = gca;
ax.FontSize = 24;

figure(4)
plot(y_dns, v2_dns, 'bo', y_node(2:n), vv(2:n), 'k-', 'LineWidth', 2);
xlabel('y (m)'); ylabel('<vv> (m^2 / s^2)'); title('Reynolds Normal Stress <vv>');
legend('DNS data', 'RSM data', 'Location', 'northeast'); legend boxoff;
ax = gca;
ax.FontSize = 24;

figure(5)
plot(y_dns, w2_dns, 'bo', y_node(2:n), ww(2:n), 'k-', 'LineWidth', 2);
xlabel('y (m)'); ylabel('<ww> (m^2 / s^2)'); title('Reynolds Normal Stress <ww>');
legend('DNS data', 'RSM data', 'Location', 'northeast'); legend boxoff;
ax = gca;
ax.FontSize = 24;

figure(6)
plot(y_dns, uv_dns, 'bo', y_node(2:n), uv(2:n), 'k-', 'LineWidth', 2);
xlabel('y (m)'); ylabel('<uv> (m^2 / s^2)'); title('Reynolds Shear Stress <uv>');
legend('DNS data', 'RSM data', 'Location', 'southeast'); legend boxoff;
ax = gca;
ax.FontSize = 24;

figure(7)
plot(y_dns, eps_dns, 'bo', y_node(2:n), epsilon(2:n), 'k-', 'LineWidth', 2);
xlabel('y (m)'); ylabel([char(949) ' (m^2 / s^3)']); title('Dissipation Rate');
legend('DNS data', 'RSM data', 'Location', 'northeast'); legend boxoff;
ax = gca;
ax.FontSize = 24;

figure(8)
plot(y_dns, dns_data(:, 3) / nu, 'bo', y_node(2:n), Pk(2:n), 'k-', 'LineWidth', 2);
xlabel('y (m)'); ylabel('P_k (m^2 / s^3)'); title('Production Rate');
legend('DNS data', 'RSM data', 'Location', 'northeast'); legend boxoff;
ax = gca;
ax.FontSize = 24;

figure(9)
plot(y_dns * 1 / nu, u_dns / 1, 'bo', y_node(2:n) / nu, U(2:n) / ustar, 'k-', y_plus, U_plus_log_law, 'r^', 'LineWidth', 2);
xlabel('y^{+}'); ylabel('U^{+}'); title('U-Velocity plotted in inner coordinates');
legend('DNS data', 'RSM data', 'Log-Law formulation', 'Location', 'southeast'); legend boxoff;
ax = gca;
ax.FontSize = 24;