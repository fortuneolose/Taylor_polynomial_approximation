% Taylor Polynomial Approximation of sin(x)^2
% This script compares different orders of Taylor approximations
% with the actual function sin(x)^2

clear; clc; close all;

% Define the range for x
x = linspace(-2*pi, 2*pi, 1000);

% Calculate the actual function
y_actual = sin(x).^2;

% Taylor series for sin(x) around x=0:
% sin(x) = x - x^3/3! + x^5/5! - x^7/7! + x^9/9! - ...
% 
% sin(x)^2 can be approximated by squaring the Taylor series
% We'll compute several orders of approximation

% Order 2: Using sin(x) ≈ x
sin_approx_2 = x;
y_taylor_2 = sin_approx_2.^2;

% Order 4: Using sin(x) ≈ x - x^3/6
sin_approx_4 = x - x.^3/6;
y_taylor_4 = sin_approx_4.^2;

% Order 6: Using sin(x) ≈ x - x^3/6 + x^5/120
sin_approx_6 = x - x.^3/6 + x.^5/120;
y_taylor_6 = sin_approx_6.^2;

% Order 8: Using sin(x) ≈ x - x^3/6 + x^5/120 - x^7/5040
sin_approx_8 = x - x.^3/6 + x.^5/120 - x.^7/5040;
y_taylor_8 = sin_approx_8.^2;

% Alternative: Direct Taylor expansion of sin(x)^2 using trig identity
% sin(x)^2 = (1 - cos(2x))/2
% cos(2x) = 1 - (2x)^2/2! + (2x)^4/4! - (2x)^6/6! + ...
% So sin(x)^2 = (1 - (1 - 2x^2 + 2x^4/3 - 4x^6/45 + ...))/2 
%             = x^2 - x^4/3 + 2x^6/45 - ...

y_direct_2 = x.^2;
y_direct_4 = x.^2 - x.^4/3;
y_direct_6 = x.^2 - x.^4/3 + 2*x.^6/45;
y_direct_8 = x.^2 - x.^4/3 + 2*x.^6/45 - x.^8/315;

% Create plots
figure('Position', [100, 100, 1200, 800]);

% Plot 1: Approximations by squaring sin(x) Taylor series
subplot(2, 2, 1);
plot(x, y_actual, 'k-', 'LineWidth', 2, 'DisplayName', 'Actual sin(x)^2');
hold on;
plot(x, y_taylor_2, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Order 2');
plot(x, y_taylor_4, 'g--', 'LineWidth', 1.5, 'DisplayName', 'Order 4');
plot(x, y_taylor_6, 'b--', 'LineWidth', 1.5, 'DisplayName', 'Order 6');
plot(x, y_taylor_8, 'm--', 'LineWidth', 1.5, 'DisplayName', 'Order 8');
grid on;
xlabel('x (radians)');
ylabel('y');
title('Taylor Approximation: Squaring sin(x) Series');
legend('Location', 'best');
xlim([-2*pi, 2*pi]);
ylim([-0.5, 2]);

% Plot 2: Direct Taylor expansion using trig identity
subplot(2, 2, 2);
plot(x, y_actual, 'k-', 'LineWidth', 2, 'DisplayName', 'Actual sin(x)^2');
hold on;
plot(x, y_direct_2, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Order 2');
plot(x, y_direct_4, 'g--', 'LineWidth', 1.5, 'DisplayName', 'Order 4');
plot(x, y_direct_6, 'b--', 'LineWidth', 1.5, 'DisplayName', 'Order 6');
plot(x, y_direct_8, 'm--', 'LineWidth', 1.5, 'DisplayName', 'Order 8');
grid on;
xlabel('x (radians)');
ylabel('y');
title('Direct Taylor Expansion (using cos identity)');
legend('Location', 'best');
xlim([-2*pi, 2*pi]);
ylim([-0.5, 2]);

% Plot 3: Error analysis for squared method
subplot(2, 2, 3);
semilogy(x, abs(y_actual - y_taylor_2), 'r-', 'LineWidth', 1.5, 'DisplayName', 'Order 2');
hold on;
semilogy(x, abs(y_actual - y_taylor_4), 'g-', 'LineWidth', 1.5, 'DisplayName', 'Order 4');
semilogy(x, abs(y_actual - y_taylor_6), 'b-', 'LineWidth', 1.5, 'DisplayName', 'Order 6');
semilogy(x, abs(y_actual - y_taylor_8), 'm-', 'LineWidth', 1.5, 'DisplayName', 'Order 8');
grid on;
xlabel('x (radians)');
ylabel('Absolute Error (log scale)');
title('Error: Squaring sin(x) Series');
legend('Location', 'best');
xlim([-2*pi, 2*pi]);

% Plot 4: Error analysis for direct method
subplot(2, 2, 4);
semilogy(x, abs(y_actual - y_direct_2), 'r-', 'LineWidth', 1.5, 'DisplayName', 'Order 2');
hold on;
semilogy(x, abs(y_actual - y_direct_4), 'g-', 'LineWidth', 1.5, 'DisplayName', 'Order 4');
semilogy(x, abs(y_actual - y_direct_6), 'b-', 'LineWidth', 1.5, 'DisplayName', 'Order 6');
semilogy(x, abs(y_actual - y_direct_8), 'm-', 'LineWidth', 1.5, 'DisplayName', 'Order 8');
grid on;
xlabel('x (radians)');
ylabel('Absolute Error (log scale)');
title('Error: Direct Taylor Expansion');
legend('Location', 'best');
xlim([-2*pi, 2*pi]);

% Display maximum errors
fprintf('Maximum Absolute Errors (Squaring Method):\n');
fprintf('  Order 2: %.6e\n', max(abs(y_actual - y_taylor_2)));
fprintf('  Order 4: %.6e\n', max(abs(y_actual - y_taylor_4)));
fprintf('  Order 6: %.6e\n', max(abs(y_actual - y_taylor_6)));
fprintf('  Order 8: %.6e\n', max(abs(y_actual - y_taylor_8)));
fprintf('\nMaximum Absolute Errors (Direct Method):\n');
fprintf('  Order 2: %.6e\n', max(abs(y_actual - y_direct_2)));
fprintf('  Order 4: %.6e\n', max(abs(y_actual - y_direct_4)));
fprintf('  Order 6: %.6e\n', max(abs(y_actual - y_direct_6)));
fprintf('  Order 8: %.6e\n', max(abs(y_actual - y_direct_8)));