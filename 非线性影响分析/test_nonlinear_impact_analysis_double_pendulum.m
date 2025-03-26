dbstop if error; % 设置断点，当发生错误时进入调试模式

% 定义双摆参数
params.m1 = 1.0;  % 第一个摆锤质量
params.m2 = 1.0;  % 第二个摆锤质量
params.l1 = 1.0;  % 第一个摆长
params.l2 = 1.0;  % 第二个摆长
params.k = 1.0;   % 弹簧刚度
params.c = 0.1;   % 阻尼系数

% 定义初始条件（二维向量）
q0 = [pi/2; 0];   % 初始角度：第一个摆垂直，第二个摆水平
dq0 = [0; 0];     % 初始角速度为零
tspan = [0, 10];

% 定义势能函数和力函数
potential = @(t, q, params) 0.5*params.k*q(1)^2 + 0.25*params.k*q(1)^4;
forces = @(t, q, dq, params) [...
    -params.k*q(1) - params.c*dq(1);  % 第一自由度广义力
    -params.k*q(2) - params.c*dq(2)   % 第二自由度广义力
];

% 双摆质量矩阵（2×2对角矩阵）
mass_matrix = @(t, q, dq, params) deal(...
    diag([...
        params.m1 * params.l1^2 + params.m2 * (params.l1^2 + params.l2^2), ... % M(1,1)
        params.m2 * params.l2^2 ...                                            % M(2,2)
    ]), ...
    zeros(2) ... % Mdot (质量矩阵是常数，导数为零)
);

% 定义约束函数
constraints = @(t, q, dq, params) double_pendulum_constraints(t, q, dq, params);

% 调用函数
[t_nonlin, z_nonlin, t_lin, z_lin] = nonlinear_impact_analysis(...
    params, tspan, q0, dq0, [], ...
    'Potential', potential, ...
    'Forces', forces, ...
    'MassMatrix', mass_matrix, ...
    'Constraints', constraints);

% 绘制结果
if ~isempty(t_nonlin)
    figure;
    plot(t_nonlin, z_nonlin(:,1)); hold on; % 非线性解 - θ1
    plot(t_nonlin, z_nonlin(:,2));         % 非线性解 - θ2
    plot(t_lin, z_lin(:,1), '--');         % 线性解 - θ1
    plot(t_lin, z_lin(:,2), '--');         % 线性解 - θ2
    legend('非线性-θ1', '非线性-θ2', '线性-θ1', '线性-θ2');
    xlabel('时间 (s)');
    ylabel('角度 (rad)');
    title('双摆角度响应');
end

% 定义双摆约束函数
function [phi, J] = double_pendulum_constraints(~, q, ~, params)
    % 双摆位置约束（末端固定在原点）
    x2 = params.l1 * sin(q(1)) + params.l2 * sin(q(1) + q(2));
    y2 = params.l1 * cos(q(1)) + params.l2 * cos(q(1) + q(2));
    phi = [x2; y2];  % 约束条件应为 phi = [0; 0]
    
    % 雅可比矩阵
    J = [...
        params.l1 * cos(q(1)) + params.l2 * cos(q(1)+q(2)), params.l2 * cos(q(1)+q(2)); 
        -params.l1 * sin(q(1)) - params.l2 * sin(q(1)+q(2)), -params.l2 * sin(q(1)+q(2))...
    ];
end
