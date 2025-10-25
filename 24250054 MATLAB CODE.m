% Clear workspace and command window
clear; clc; close all;

% Parameters for simulation
params.rs = 0.5;    % Growth rate of trees
params.rE = 0.3;    % Energy recovery rate
params.KS = 100;    % Carrying capacity for tree height
params.KE = 150;     % Carrying capacity for energy
params.P = 0.1;     % Impact parameter of budworms
params.B = 150;      % Budworm population (high value as requested)

% Simulation parameters
S0_range = linspace(20, 80, 4);  % Initial tree heights
E0_range = linspace(20, 60, 4);  % Initial energy levels
tspan = [0, 100];
dt = 0.1;

% Generate phase portrait
disp('Generating phase portrait...');
figure('Position', [100, 100, 800, 600]);
hold on;

for S0 = S0_range
    for E0 = E0_range
        [~, S, E] = simulate_system(params, S0, E0, tspan, dt);
        plot(S, E, '-', 'LineWidth', 1);
    end
end

xlabel('Tree Height (S)');
ylabel('Forest Energy Reserve (E)');
title('Phase Portrait of Forest-Budworm System');
grid on;
hold off;

% Single trajectory simulation
disp('Simulating single trajectory...');
[t, S, E] = simulate_system(params, 50, 40, tspan, dt);

% Plot time series
figure('Position', [100, 100, 1000, 400]);

% Tree height plot
subplot(1, 2, 1);
plot(t, S, 'LineWidth', 1.5);
xlabel('Time');
ylabel('Height');
title('Tree Height Over Time');
grid on;

% Energy reserve plot
subplot(1, 2, 2);
plot(t, E, 'LineWidth', 1.5);
xlabel('Time');
ylabel('Energy');
title('Forest Energy Over Time');
grid on;

% Adjust subplot spacing
set(gcf, 'Position', get(0, 'Screensize').*[1 1 0.8 0.8]);

% Display some analysis results
disp('Analysis Results:');
fprintf('Final tree height: %.2f\n', S(end));
fprintf('Final energy reserve: %.2f\n', E(end));
fprintf('Maximum tree height: %.2f\n', max(S));
fprintf('Minimum energy reserve: %.2f\n', min(E));

% Save figures if needed
saveas(1, 'phase_portrait.png');
saveas(2, 'time_series.png');

% Helper functions
function [t, S, E] = simulate_system(params, S0, E0, tspan, dt)
    % Simulate the system using RK4
    t = tspan(1):dt:tspan(2);
    n_steps = length(t);
    
    % Initialize arrays
    S = zeros(1, n_steps);
    E = zeros(1, n_steps);
    
    % Set initial conditions
    S(1) = S0;
    E(1) = E0;
    
    % Integrate using RK4
    for i = 2:n_steps
        state = [S(i-1); E(i-1)];
        new_state = runge_kutta_4(params, state, t(i-1), dt);
        S(i) = new_state(1);
        E(i) = new_state(2);
    end
end

function new_state = runge_kutta_4(params, state, t, dt)
    % 4th order Runge-Kutta method
    k1 = derivatives(params, t, state);
    k2 = derivatives(params, t + dt/2, state + dt*k1/2);
    k3 = derivatives(params, t + dt/2, state + dt*k2/2);
    k4 = derivatives(params, t + dt, state + dt*k3);
    
    new_state = state + dt*(k1 + 2*k2 + 2*k3 + k4)/6;
end

function dstate = derivatives(params, ~, state)
    % Calculate derivatives for the system
    S = state(1);
    E = state(2);
    
    dSdt = params.rs * S * (1 - (params.KE * S)/(params.KS * E));
    dEdt = params.rE * E * (1 - E/params.KE) - params.P * params.B/S;
    
    dstate = [dSdt; dEdt];
end