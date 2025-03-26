clear all; close all; clc;

% 定义参数
m = 0.999;  % 测试接近1的情况
L = 1;
x = linspace(0, L, 1000);

% 计算 theta(x)，确保定义域合法
arg = (m*L/2) * (2*(x/L) - (x/L).^2);
if any(abs(arg) > 1)
    error('theta(x) is complex: adjust m or x.');
end
theta_x = asin(arg);

% 方法1：直接调用内置函数（默认双精度）
F = ellipticF(theta_x, m);
E = ellipticE(theta_x, m);
y_default = (F - E) / sqrt(m);

% 方法2：符号计算高精度结果（可选）
syms theta m_sym;
F_sym = ellipticF(theta, m_sym);
E_sym = ellipticE(theta, m_sym);
y_sym = (F_sym - E_sym) / sqrt(m_sym);
y_high_precision = double(subs(y_sym, {theta, m_sym}, {theta_x, m}));

% 比较两种方法的差异
error = abs(y_default - y_high_precision);
fprintf('Max error between default and high-precision: %.2e\n', max(error));

% 绘图
plot(x, y_default, 'b-', x, y_high_precision, 'r--');
legend('Default precision', 'High precision');
xlabel('x'); ylabel('y(x)');
title('Comparison of Calculation Precision');