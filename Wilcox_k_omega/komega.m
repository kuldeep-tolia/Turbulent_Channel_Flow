% DNS data at Re_delta=7890, Re_tau=395 (Moin, Kim & Mansour, PoF, 1999).
% All quantites are normalized by u_tau and nu unless stated otherwise.
% Delta denotes the channel half-width.

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
ustar = 1;           % Wall friction velocity
rho = 1;             % Fluid density
kappa = 0.41;        % Von Karmann constant 
dPdx = -1.00;        % Pressure gradient
delta = 1.0;         % Channel half-width

iter = 1;            % Iteration counter
m = length(y_dns);   % # of points in array = 97

% k-omega model constants

beta = 0.09;
c1 = 5.0 / 9.0;
c2 = 3.0 / 40.0;
sigma_k = 2.0;
sigma_w = 2.0;

% Residual error limit

residue_limit = 10^(-4);

% Under relaxation factor

urf = 0.7;

% Grid (based on DNS data)
% y_node stores location of cell nodal values
% y_face stores location of cell face values 

y_node = [y_dns; (2.0 - flip(y_dns))];
y_node(m) = [];
n = length(y_node);

y_face(1) = y_node(1);
y_face(2) = 2 * (y_node(2) - y_node(1));

for i = 3:n-1
  
  y_face(i) = 2 * y_node(i) - y_face(i-1);
  
end

y_face = y_face';

% Initialise variables

k = ones(n, 1);
k_old = ones(n, 1);
w = ones(n, 1);
w_old = ones(n, 1);
U = ones(n, 1);
U_old = ones(n, 1);
nu_t = ones(n, 1);
nu_t_face = ones(n-1, 1);
dUdy = ones(n, 1);
lm = ones(n, 1);
Pk = ones(n, 1);
dkdy = ones(n, 1);
d2kdy2 = ones(n, 1);

% Calculating mixing length based on Prandtl's Mixing Length model
% Adapted from H. Versteeg and W. Malalasekera 1995

lm = delta * (0.14 - (0.08 * (1 - y_node / delta) .^ 2) - ...
      (0.06 * (1 - y_node / delta) .^ 4));
      
% Initial conditions for old & new variables

k(:) = 5 * 10^-3;
k_old(:) = k(:);
U(:) = rand(n, 1);
U_old(:) = U(:);
w(:) = 10^-5;
w_old(:) = w(:);
residue = 1.0;
num1 = 0.0;
den1 = 0.0;

% Boundary condition

U(1) = 0.0;                                     % velocity wall bc
U(n) = 0.0;                                     % velocity wall bc
k(1) = 0.0;                                     % TKE wall bc
k(n) = 0.0;                                     % TKE wall bc
w(2) = (6.0 * nu) / (c2 * y_node(2) ^ 2);       % omega bc for near wall cell
w(1) = w(2);                                    % interpolating omega to wall
w(n-1) = (6.0 * nu) / (c2 * (2.0 - y_node(n-1)) ^ 2);     % omega bc for near wall cell
w(n) = w(n-1);                                  % interpolating omega to wall

while residue > residue_limit
  
  % Compute eddy viscosity
  
  if (iter <= 300)
    
    dUdy(1) = (U(2) - U(1)) / (y_node(2) - y_node(1));          % forward difference scheme
    dUdy(n) = (U(n-1) - U(n)) / (y_node(n) - y_node(n-1));      % backward difference scheme
    
    for i = 2:n-1
      
      dUdy(i) = (U(i+1) - U(i-1)) / (y_node(i+1) - y_node(i-1));    % central difference scheme
    
    end
  
    for i = 1:n
        
      nu_t(i) = lm(i) * lm(i) * abs(dUdy(i));       % eddy viscosity based on mixing length model
    
    end
  
  else
  
    for i = 1:n
    
      nu_t(i) = k(i) / w(i);        % eddy viscosity coupled with k and omega
   
    end
  
  end
  
  %----------------------------- NUMERICAL TIP --------------------------
  % Often it can be tricky to start the simulations. They often diverge.
  % It can then be useful to compute the turbulence viscosity from the 
  % mixing-length model for the first few iterations (say 100?).
  % In this way the k & omega model equations are de-coupled from U for the first few iterations.

  % Computing nu_t at faces using harmonic mean approach
  
  nu_t_face(1) = nu_t(1);
  nu_t_face(n-1) = nu_t(n);
  
  for i = 2:n-2
    
    nu_t_face(i) = (2.0 * (y_node(i+1) - y_node(i))) / (((y_face(i) - y_face(i-1)) / nu_t(i)) + ((y_face(i+1) - y_face(i)) / nu_t(i+1)));
    
  end
  
  % Compute U
  
  for i = 2:n-1         % Gauss-Seidel loop
    
    if (i == 2)         % Bottom wall cell treatment
      
      aN = (nu + nu_t_face(i)) / (y_node(i+1) - y_node(i));
      aS = (nu + nu_t_face(i-1)) / (y_node(i) - y_node(i-1));
      b = 1.0 * (y_face(i) - y_face(i-1));
      aP = aN + aS;
      U(i) = ((aN * U(i+1)) + b) / aP;
      
    elseif (i == (n-1)) % Top wall cell treatment
      
      aN = (nu + nu_t_face(i)) / (y_node(i+1) - y_node(i));
      aS = (nu + nu_t_face(i-1)) / (y_node(i) - y_node(i-1));
      b = 1.0 * (y_face(i) - y_face(i-1));
      aP = aN + aS;
      U(i) = ((aS * U(i-1)) + b) / aP;
      
    else                % Interior cell treatment
      
      aN = (nu + nu_t_face(i)) / (y_node(i+1) - y_node(i));
      aS = (nu + nu_t_face(i-1)) / (y_node(i) - y_node(i-1));
      b = 1.0 * (y_face(i) - y_face(i-1));
      aP = aN + aS;
      U(i) = ((aN * U(i+1)) + (aS * U(i-1)) + b) / aP;
      
    end
    
  end
  
  U(1) = 0.0;           % velocity wall bc
  U(n) = 0.0;           % velocity wall bc
  
  % Compute dUdy
  
  dUdy(1) = (U(2) - U(1)) / (y_node(2) - y_node(1));          % forward difference scheme
  dUdy(n) = (U(n-1) - U(n)) / (y_node(n) - y_node(n-1));      % backward difference scheme
    
  for i = 2:n-1
      
    dUdy(i) = (U(i+1) - U(i-1)) / (y_node(i+1) - y_node(i-1));     % central difference scheme
    
  end
  
  % Compute Pk
  
  for i = 1:n
    
    Pk(i) = nu_t(i) * (dUdy(i) ^ 2);
    
  end
  
  % Compute k
  
  for i = 2:n-1         % Gauss-Seidel loop
    
    if (i == 2)         % Bottom wall cell treatment
      
      aN = (nu + (nu_t_face(i) / sigma_k)) / (y_node(i+1) - y_node(i));
      aS = (nu + (nu_t_face(i-1) / sigma_k)) / (y_node(i) - y_node(i-1));
      b = Pk(i) * (y_face(i) - y_face(i-1));
      aP = aN + aS - ((-1.0 * beta * w(i)) * (y_face(i) - y_face(i-1)));
      k(i) = ((aN * k(i+1)) + b) / aP;
      
    elseif (i == (n-1)) % Top wall cell treatment
      
      aN = (nu + (nu_t_face(i) / sigma_k)) / (y_node(i+1) - y_node(i));
      aS = (nu + (nu_t_face(i-1) / sigma_k)) / (y_node(i) - y_node(i-1));
      b = Pk(i) * (y_face(i) - y_face(i-1));
      aP = aN + aS - ((-1.0 * beta * w(i)) * (y_face(i) - y_face(i-1)));
      k(i) = ((aS * k(i-1)) + b) / aP;
      
    else                % Interior cell treatment
      
      aN = (nu + (nu_t_face(i) / sigma_k)) / (y_node(i+1) - y_node(i));
      aS = (nu + (nu_t_face(i-1) / sigma_k)) / (y_node(i) - y_node(i-1));
      b = Pk(i) * (y_face(i) - y_face(i-1));
      aP = aN + aS - ((-1.0 * beta * w(i)) * (y_face(i) - y_face(i-1)));
      k(i) = ((aN * k(i+1)) + (aS * k(i-1)) + b) / aP;
      
    end
    
  end
  
  k(1) = 0.0;         % TKE wall bc
  k(n) = 0.0;         % TKE wall bc
  
  % Compute omega
  
  for i = 3:n-2        % Interior cell treatment
    
    aN = (nu + (nu_t_face(i) / sigma_w)) / (y_node(i+1) - y_node(i));
    aS = (nu + (nu_t_face(i-1) / sigma_w)) / (y_node(i) - y_node(i-1));
    b = (c1 * Pk(i) * w(i) / k(i)) * (y_face(i) - y_face(i-1));
    aP = aN + aS - ((-c2 * w(i)) * (y_face(i) - y_face(i-1)));
    w(i) = ((aN * w(i+1)) + (aS * w(i-1)) + b) / aP;
      
  end
  
  w(2) = (6.0 * nu) / (c2 * y_node(2) ^ 2);                 % omega bc for near wall cell 
  w(1) = w(2);                                              % interpolating omega to wall
  w(n-1) = (6.0 * nu) / (c2 * (2.0 - y_node(n-1)) ^ 2);     % omega bc for near wall cell
  w(n) = w(n-1);                                            % interpolating omega to wall
  
  % Under-relaxing variables
  
  U(:) = (urf * U(:)) + ((1.0 - urf) * U(:));
  k(:) = (urf * k(:)) + ((1.0 - urf) * k(:));
  w(:) = (urf * w(:)) + ((1.0 - urf) * w(:));
  
  for i = 1:n         % Calculate l2norm for velocity
    
    num1 = num1 + (U(i) - U_old(i)) ^ 2;
    den1 = den1 + (U(i) ^ 2);
    
  end
  
  residue = sqrt(num1) / sqrt(den1);
    
  U_old(:) = U(:);
  k_old(:) = k(:);
  w_old(:) = w(:);
  
  iter = iter + 1;
  
  t = mod(iter, 100);
  
  if (t == 0)         % Print residuals and other quantities at cell(10)
    
    fprintf('iter = %i\n', iter)
    fprintf('residue = %.7f\n', residue)
    fprintf('U = %.7f\n', U(10))
    fprintf('nut = %.7f\n', nu_t(10))
    fprintf('dUdy = %.7f\n', dUdy(10))
    fprintf('k = %.7f\n', k(10))
    fprintf('Pk = %.7f\n', Pk(10))
    fprintf('w = %.7f\n', w(10));
    fprintf('\n')
  
  end

end
      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 6 columns in dns_data.dat as below:
%
%      y+         Diss        prod     vel_p_grad   Turb_diff   Visc_diff
%
% Please note that all terms are normalized by ustar^4/nu

k_dns = 0.5 * (u2_dns + v2_dns + w2_dns);         % Calculate TKE from DNS data
eps_dns = dns_data(:,2) * ustar ^ 4 / nu;         % eps is normalized by ustar^4/nu

% Calculate dkdy

dkdy(1) = (k(2) - k(1)) / (y_node(2) - y_node(1));          % forward difference scheme
dkdy(n) = (k(n-1) - k(n)) / (y_node(n) - y_node(n-1));      % backward difference scheme

for i = 2:n-1
  
  dkdy(i) = (k(i+1) - k(i-1)) / (y_node(i+1) - y_node(i-1));     % calculate difference scheme
  
end

% Calculate d2kdy2

d2kdy2(1) = (dkdy(2) - dkdy(1)) / (y_node(2) - y_node(1));        % forward difference scheme
d2kdy2(n) = (dkdy(n-1) - dkdy(n)) / (y_node(n) - y_node(n-1));    % backward difference scheme

for i = 2:n-1
  
  d2kdy2(i) = (dkdy(i+1) - dkdy(i-1)) / (y_node(i+1) - y_node(i-1));     % central difference scheme
  
end

dns_diff = (dns_data(:, 5) / nu) + (dns_data(:, 6) / nu) + (dns_data(:, 4) / nu);   % Total diffusion rate for DNS data
model_diff = (-nu_t(1:m) .* d2kdy2(1:m) / sigma_k) + (nu * d2kdy2(1:m));            % Modelled total diffusion rate
P_omega = c1 * (dUdy(:).^2);          % Production rate of omega

% Plotting quantities

figure(1)
plot(u_dns, y_node(1:m), 'bo', U(1:m), y_node(1:m), 'k-', 'LineWidth', 2);
xlabel('U'); ylabel('y / h'); title('U-velocity');
legend('DNS data', 'k-\omega model','Location','northeast'); legend boxoff;
h=get(gcf, "currentaxes");
set(h, "fontsize", 20);

figure(2)
plot(y_node(1:m), k_dns, 'bo', y_node(1:m), k(1:m), 'k-', 'LineWidth', 2);
xlabel('y / h'); ylabel('k'); title('Turbulence kinetic energy');
legend('DNS data', 'k-\omega model', 'Location', 'northeast'); legend boxoff;
h=get(gcf, "currentaxes");
set(h, "fontsize", 20);

figure(3)
plot(y_node(1:m), eps_dns, 'bo', y_node(1:m), (beta .* k(1:m) .* w(1:m)), 'k-', 'LineWidth', 2);
xlabel('y / h'); ylabel('\epsilon'); title('Dissipation rate of k');
legend('DNS data', 'k-\omega model', 'Location', 'northeast'); legend boxoff;
h=get(gcf, "currentaxes");
set(h, "fontsize", 20);

figure(4)
plot(y_dns, -uv_dns, 'bo', y_node(1:m), (Pk(1:m) ./ abs(dUdy(1:m))), 'k-', 'LineWidth', 2)
xlabel('y / h'); ylabel('-<uv>'); title('Turbulence shear stress');
legend('DNS data', 'k-\omega model', 'Location', 'northeast'); legend boxoff;
h=get(gcf, "currentaxes");
set(h, "fontsize", 20);

figure(5)
plot(y_dns, dns_data(:,3) / nu, 'bo', y_node(1:m), (Pk(1:m)), 'k-', 'LineWidth', 2);
xlabel('y / h'); ylabel('P_k'); title('Production rate of k');
legend('DNS data', 'k-\omega model', 'Location', 'northeast'); legend boxoff;
h=get(gcf, "currentaxes");
set(h, "fontsize", 20);

figure(6)
plot(y_dns, dns_data(:,5) / nu, 'bo', y_node(1:m), (nu_t(1:m) .* d2kdy2(1:m) / sigma_k), 'k-', 'Linewidth', 2.0);
xlabel('y / h'); ylabel('Turbulent diffusion rate'); title('Turbulent diffusion of k');
legend('DNS data', 'k-\omega model', 'Location', 'northeast'); legend boxoff;
h=get(gcf, "currentaxes");
set(h, "fontsize", 20);

figure(7)
plot(y_dns, dns_data(:,6) / nu, 'bo', y_node(1:m), nu * d2kdy2(1:m), 'k-', 'Linewidth', 2.0);
xlabel('y / h'); ylabel('Viscous diffusion rate'); title('Viscous diffusion of k');
legend('DNS data', 'k-\omega model', 'Location', 'northeast'); legend boxoff;
h=get(gcf, "currentaxes");
set(h, "fontsize", 20);

figure(8)
plot(y_dns, dns_diff, 'bo', y_node(1:m), model_diff, 'k-', 'Linewidth', 2.0);
xlabel('y / h'); ylabel('Total diffusion rate'); title('Total diffusion of k');
legend('DNS data', 'k-\omega model', 'Location', 'northeast'); legend boxoff;
h=get(gcf, "currentaxes");
set(h, "fontsize", 20);

figure(9)
plot(y_node(1:m), nu_t(1:m), 'k-', 'LineWidth', 2);
xlabel('y / h'); ylabel('\nu_t'); title('Eddy viscosity');
legend('k-\omega model', 'Location', 'northeast'); legend boxoff;
h=get(gcf, "currentaxes");
set(h, "fontsize", 20);

figure(10)
plot(y_node(1:m), P_omega(1:m), 'k-','Linewidth', 2);
xlabel('y / H'); ylabel('Production rate'); title('Production rate of \omega');
h=get(gcf, "currentaxes");
set(h, "fontsize", 20);

% Calculate wall shear stress

tau_w = dUdy(1) * nu;
fprintf('\n\n wall shear stress = %.7f\n\n', tau_w);
