function [y_node, y_face, n, U, k, epsilon, uu, vv, ww, uv, ustar] = k_epsilon(nu, rho, dPdx, delta, ny)

% Read DNS data [half-channel is given (till centerline)]

load y_dns.dat
load u_dns.dat
load u2_dns.dat
load v2_dns.dat
load w2_dns.dat
load uv_dns.dat
load dns_data.dat

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

kappa = 0.41;        % Von Karmann constant 
E = 9.0;

% k-epsilon model constants

c_mu = 0.09;                 
c1_ke = 1.44;
c2_ke = 1.92;
sigma_k = 1.0;         
sigma_eps = 1.3;  
c1_eps = 1.44;
c2_eps = 1.92;

% Grid
% y_dns(26) corresponds to a point yplus >= 30

y_node(1, 1) = 0.0;
y_node(2, 1) = 0.1;                         % placement of first node near wall @ y = 0.1
dy = (delta - (2.0 * y_node(2))) / ny;
y_node(3) = (2.0 * y_node(2)) + (0.5 * dy);

for i = 4:ny+2
  
  y_node(i) = y_node(i-1) + dy;
  
end

y_node(ny+3) = 1.0;
y_face(1, 1) = 0.0;
y_face(2, 1) = 2.0 * y_node(2);

for i = 3:ny+2
  
  y_face(i) = y_face(i-1) + dy;
  
end

n = length(y_node);

% Initialise variables

k = ones(n, 1);
epsilon = ones(n, 1);
U = ones(n, 1);
U_old = ones(n, 1);
nu_t = ones(n, 1);
nu_t_face = ones(n-1, 1);
dUdy = ones(n, 1);
lm = ones(n, 1);
Pk = ones(n, 1);
uu = ones(n, 1);
vv = ones(n, 1);
ww = ones(n, 1);
uv = ones(n, 1);
residue = 1.0;
tolerance = 10^-4;
numerator = 0.0;
denominator = 0.0;

% Calculating mixing length based on Prandtl's Mixing Length model
% Adapted from H. Versteeg and W. Malalasekera 1995

lm = delta * (0.14 - (0.08 * (1 - y_node / delta) .^ 2) - ...
      (0.06 * (1 - y_node / delta) .^ 4));
      
% Initial conditions for old & new variables

k(:) = 5 * 10^-3;
U(:) = rand(n, 1);
U_old(:) = U(:);
epsilon(:) = 10^-5;
ustar = rand(1);                % guess value for ustar between [0, 1]

% Boundary condition

U(1) = 0.0;                                     % velocity wall bc
U(n) = U(n-1);                                  % velocity symmetry bc
k(1) = 0.0;                                     % tke wall bc
k(2) = (c_mu ^ (-0.5)) * (ustar ^ 2);           % tke wall function condition
k(n) = k(n-1);                                  % tke symmetry bc
epsilon(2) = (ustar ^ 3) / (kappa * y_node(2)); % epsilon wall function condition
epsilon(n) = epsilon(n-1);                      % epsilon symmetry bc

%while (residue >= tolerance)
 for iter = 1:200
    
    for i = 2:n
      
      nu_t(i) = c_mu * k(i) * k(i) / epsilon(i);    % eddy viscosity coupled with k and epsilon
     
    end
    
    nu_t(1) = 0.0;
  
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
      b = -dPdx * (y_face(i) - y_face(i-1));
      aP = aN + aS;
      U(i) = ((aN * U(i+1)) + b) / aP;
      
    elseif (i == (n-1)) % Top wall cell treatment
      
      aN = 0.0;
      aS = (nu + nu_t_face(i-1)) / (y_node(i) - y_node(i-1));
      b = -dPdx * (y_face(i) - y_face(i-1));
      aP = aS;
      U(i) = ((aS * U(i-1)) + b) / aP;
      
    else
      
      aN = (nu + nu_t_face(i)) / (y_node(i+1) - y_node(i));
      aS = (nu + nu_t_face(i-1)) / (y_node(i) - y_node(i-1));
      b = -dPdx * (y_face(i) - y_face(i-1));
      aP = aN + aS;
      U(i) = ((aN * U(i+1)) + (aS * U(i-1)) + b) / aP;
      
    end
    
  end
  
  U(1) = 0.0;         % velocity wall bc
  U(n) = U(n-1);      % velocity symmetry bc
  
  % Compute dUdy
  
  dUdy(1) = (U(2) - U(1)) / (y_node(2) - y_node(1));          % forward difference scheme
  dUdy(n) = 0.0;
  
  for i = 2:n-1
    
    dUdy(i) = (U(i+1) - U(i-1)) / (y_node(i+1) - y_node(i-1));     % central difference scheme
    
  end
  
  % Compute Pk
  
  for i = 2:n
    
    Pk(i) = nu_t(i) * (dUdy(i) ^ 2);
    
  end
  
  % Compute ustar
  
  ustar = (kappa * U(2)) / log(E * ustar * y_node(2) / nu);
  
  % Compute k
  
  for i = 3:n-1          % Gauss-Seidel loop
    
    if (i == (n-1))
      
      aN = 0.0;
      aS = (nu + (nu_t_face(i-1) / sigma_k)) / (y_node(i) - y_node(i-1));
      b = (Pk(i) - epsilon(i)) * (y_face(i) - y_face(i-1));
      aP = aS;
      k(i) = ((aS * k(i-1)) + b) / aP;
      
    else
      
      aN = (nu + (nu_t_face(i) / sigma_k)) / (y_node(i+1) - y_node(i));
      aS = (nu + (nu_t_face(i-1) / sigma_k)) / (y_node(i) - y_node(i-1));
      b = (Pk(i) - epsilon(i)) * (y_face(i) - y_face(i-1));
      aP = aN + aS;
      k(i) = ((aN * k(i+1)) + (aS * k(i-1)) + b) / aP;
      
    end
    
  end
  
  k(1) = 0.0;                           % tke wall bc
  k(2) = (c_mu ^ -0.5) * (ustar ^ 2);   % tke wall function condition
  k(n) = k(n-1);                        % tke symmetry bc
  
  % Compute epsilon
  
  for i = 3:n-1          % Gauss-Seidel loop
    
    if (i == (n-1))
      
      aN = 0.0;
      aS = (nu + (nu_t_face(i-1) / sigma_eps)) / (y_node(i) - y_node(i-1));
      b = (c1_eps * Pk(i) * epsilon(i) / k(i)) * (y_face(i) - y_face(i-1));
      aP = aS - ((-c2_eps * epsilon(i) / k(i)) * (y_face(i) - y_face(i-1)));
      epsilon(i) = ((aS * epsilon(i-1)) + b) / aP;

    else

      aN = (nu + (nu_t_face(i) / sigma_eps)) / (y_node(i+1) - y_node(i));
      aS = (nu + (nu_t_face(i-1) / sigma_eps)) / (y_node(i) - y_node(i-1));
      b = (c1_eps * Pk(i) * epsilon(i) / k(i)) * (y_face(i) - y_face(i-1));
      aP = aN + aS - ((-c2_eps * epsilon(i) / k(i)) * (y_face(i) - y_face(i-1)));
      epsilon(i) = ((aN * epsilon(i+1)) + (aS * epsilon(i-1)) + b) / aP;
      
    end
    
  end
  
  epsilon(2) = (ustar ^ 3) / (kappa * y_node(2)); % epsilon wall function condition
  epsilon(n) = epsilon(n-1);                      % epsilon symmetry bc
  epsilon(1) = epsilon(2);

end

% Compute initial field for Reynold stresses

uu = (2 / 3) * k;
vv = (2 / 3) * k;
ww = (2 / 3) * k;
uv = -nu_t .* abs(dUdy);

end
