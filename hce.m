clear;
clc;

N = 20; % resolution (increase for finer grid)
L = 0.1; % size
h = L / N; % size of cell

x = linspace(0, L, N);
y = linspace(0, L, N);

[X, Y] = meshgrid(x, y);

% Basic conditions
T = zeros(size(X));
T0 = 300; % initial temperature

T = T + T0;

% Solver settings
t = 0; % start time
tmax = 60; % simulation time
Nt = 1000; % number of time steps
tau = (tmax - t) / Nt; % time step

% Boundary conditions (fixed)
Tl = 400; % left boundary temperature
Tb = 600; % bottom boundary temperature
T(:, 1) = Tl; % set left boundary
T(1, :) = Tb; % set bottom boundary

% Physical properties
rho = 7800; % density
cp = 460; % specific heat capacity
k = 46; % thermal conductivity

A = k * tau / (rho * cp * h^2);

% Helper arrays for the Thomas algorithm (tridiagonal matrix)
a = -A * ones(N-2, 1);
b = (1 + 2*A) * ones(N-2, 1);
c = -A * ones(N-2, 1);

% Set the range of temperatures for the color scale
Tmin = min([T0, Tl, Tb]);
Tmax = max([T0, Tl, Tb]);

% Time loop
while t < tmax
    Tprev = T;  % save the previous temperature distribution
    
    % Step 1: Update temperature along rows (horizontal)
    for iy = 2:N-1
        % Construct the right-hand side (known terms from previous time step)
        d = Tprev(iy, 2:N-1);
        d(1) = d(1) + A * T(iy, 1);   % left boundary condition
        d(end) = d(end) + A * T(iy, end); % right boundary (adiabatic)

        % Solve the tridiagonal system using Thomas algorithm
        T(iy, 2:N-1) = thomas_algorithm(a, b, c, d);
    end

    % Step 2: Update temperature along columns (vertical)
    for ix = 2:N-1
        % Construct the right-hand side (known terms from previous time step)
        d = Tprev(2:N-1, ix);
        d(1) = d(1) + A * T(1, ix);   % bottom boundary condition
        d(end) = d(end) + A * T(end, ix); % top boundary (adiabatic)

        % Solve the tridiagonal system using Thomas algorithm
        T(2:N-1, ix) = thomas_algorithm(a, b, c, d);
    end
    
    % Reapply boundary conditions (fix the temperatures on the boundaries)
    T(:, 1) = Tl;   % left boundary fixed at Tl
    T(1, :) = Tb;   % bottom boundary fixed at Tb
    T(:, end) = T(:, end-1); % right boundary: adiabatic condition (no heat flux)
    T(end, :) = T(end-1, :); % top boundary: adiabatic condition (no heat flux)
    
    % Visualization (fix the color scale to Tmin and Tmax)
    contourf(X, Y, T, 20, 'EdgeColor', 'none');
    caxis([Tmin Tmax]);  % fix the color axis
    xlabel('x, m');
    ylabel('y, m');
    title(sprintf('Температура пластинки, K, t = %f', t));
    colorbar;
    drawnow;
    
    % Advance time
    t = t + tau;
end

% Function for solving tridiagonal system using the Thomas algorithm
function x = thomas_algorithm(a, b, c, d)
    N = length(b); % size of the system
    % Forward sweep
    for i = 2:N
        w = a(i-1) / b(i-1);
        b(i) = b(i) - w * c(i-1);
        d(i) = d(i) - w * d(i-1);
    end
    % Backward substitution
    x = zeros(N, 1);
    x(N) = d(N) / b(N);
    for i = N-1:-1:1
        x(i) = (d(i) - c(i) * x(i+1)) / b(i);
    end
end
