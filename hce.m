clear;
clc;

N = 20; % resolution
L = 0.1; % size
h = L / N; % size of cell

x = linspace(0, L, N);
y = linspace(0, L, N);

[X, Y] = meshgrid(x, y);

% Basic conditions

T = zeros(size(X));
T0 = 300; % basic temp

T = T + T0;

% Solver settings

t = 0; % basic time
tmax = 60; % end time
Nt = 1000; % step quantity by time
tau = (tmax - t) / Nt;


% Border conditions

Tl = 400;
Tb = 600;
T(:, 1) = Tl;
T(1, :) = Tb;

Tprev = T;

rho = 7800; % density
cp = 460; % heat capacity
k = 46; % thermal conductivity

A = k * tau / (rho * cp * h^2);

while t < tmax

    for iy = 2:N-1
       for ix = 2:N-1 
           sum4 = T(iy + 1, ix) + T(iy - 1, ix) + T(iy, ix + 1) + T(iy, ix - 1);

           T(iy, ix) = A * sum4 + (1 - 4*A) * Tprev(iy, ix);
       end
    end

    Tprev = T;

    % Result visualization
    contourf(X, Y, T, 20, 'EdgeColor','none');
    xlabel('x, m');
    ylabel('y, m');
    %title('Температура пластинки, K');
    title(sprintf('Температура пластинки, K, t = %f', t));
    colorbar;
    drawnow;

    t = t + tau;
end






