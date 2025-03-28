% ODE数值解法
clear; clc;

%% 测试案例选择 (取消注释需要测试的案例)
% case_num = 1;  % 变系数线性方程
case_num = 2;  % Logistic方程
% case_num = 3;  % 刚性方程(可能计算一分钟左右)
% case_num = 4;    % 二阶振荡方程（向量形式）
% case_num = 5;    % Riccati方程

%% 根据选择初始化案例（维度统一为列向量）
switch case_num
    case 1
        % 标量方程
        f = @(t,y) -t*y;
        y_exact = @(t) exp(-t.^2/2);
        t0 = 0; tf = 2; h = 0.1;
        y0 = 1;
        t = t0:h:tf;
        vec_flag = false;
        
    case 2
        % 标量方程
        f = @(t,y) y*(1 - y);
        y_exact = @(t) 1./(1 + exp(-t));
        t0 = 0; tf = 5; h = 0.0001;
        y0 = 0.5;
        t = t0:h:tf;
        vec_flag = false;
        
    case 3
        % 标量方程
        f = @(t,y) -1000*(y - cos(t));
        y_exact = @(t) exp(-1000*t) + (1000/(1000001))*cos(t) + (1/1000001)*sin(t);
        t0 = 0; tf = 0.1; h = 1e-6;
        y0 = 1;
        t = t0:h:tf;
        vec_flag = false;
        
    case 4
        % 向量方程（确保初始值为列向量）
        f = @(t,Y) [Y(2); -4*Y(1)]; % Y = [y; y']
        y_exact = @(t) [cos(2*t), -2*sin(2*t)]'; % 返回列向量
        t0 = 0; tf = 2*pi; h = 0.01;
        y0 = [1; 0]; % 列向量
        t = t0:h:tf;
        vec_flag = true;
        
    case 5
        % 标量方程
        f = @(t,y) y^2 - 2/t^2;
        y_exact = @(t) 1./t + 1./(t.*log(t));
        t0 = 1; tf = 3; h = 0.05;
        y0 = 1;
        t = t0:h:tf;
        vec_flag = false;
end

%% 统一调用数值解法（强化维度检查）
if ~vec_flag % 标量方程
    % 确保初始值转换为行向量以匹配时间序列存储
    y0 = y0(:)'; 
    
    % 调用数值方法
    y_be = ode_solver(@backward_euler, f, y0, t, h);
    y_ie = ode_solver(@improved_euler, f, y0, t, h);
    y_adams = ode_solver(@adams_pc, f, y0, t, h);
    
    % 精确解计算
    y_exact_vals = arrayfun(y_exact, t);
    
    % 误差分析
    analyze_errors(t, y_be, y_ie, y_adams, y_exact_vals, case_num);
    
else % 向量方程
    % 确保初始值为列向量
    y0 = y0(:);
    
    % 调用向量版求解器
    Y_be = ode_solver_vec(@backward_euler_vec, f, y0, t, h);
    Y_ie = ode_solver_vec(@improved_euler_vec, f, y0, t, h);
    Y_adams = ode_solver_vec(@adams_pc_vec, f, y0, t, h);
    
    % 精确解计算（返回矩阵：每列对应一个时间点）
    Y_exact = cell2mat(arrayfun(y_exact, t, 'UniformOutput', false));
    
    % 误差分析
    analyze_vector_errors(t, Y_be, Y_ie, Y_adams, Y_exact, case_num);
end

% 精度验证代码（添加到脚本末尾）
fprintf('\n案例%d的精度验证：\n', case_num);
fprintf('解析解表达式：%s\n', char(y_exact)); 
fprintf('在t=%.2f处：\n', tf);
fprintf('双精度计算：%.15f\n', y_exact(tf));
if exist('y_exact_sol','var')
    fprintf('符号工具箱计算：%.15f\n', double(subs(y_exact, t, tf)));
    err = abs(y_exact_1(tf) - double(subs(y_exact, t, tf)));
    fprintf('差异：%.3e (%.3f eps)\n', err, err/eps);
end

%% 通用标量方程求解器（确保维度一致性）
%%
function y = ode_solver(method, f, y0, t, h)
    n = length(t);
    y = zeros(size(y0,2), n); % 行向量存储
    y(:,1) = y0;
    for i = 1:n-1
        y(:,i+1) = method(f, y(:,i), t(i), h);
    end
end

%% 通用向量方程求解器（列向量存储）
function Y = ode_solver_vec(method, f, Y0, t, h)
    n = length(t);
    dim = length(Y0);
    Y = zeros(dim, n); % 列向量存储
    Y(:,1) = Y0;
    for i = 1:n-1
        Y(:,i+1) = method(f, Y(:,i), t(i), h);
    end
end

%% 数值方法定义（严格维度控制）
% ==== 标量方法 ====
% % 预热运行（避免首次运行编译开销）
% f = @(t,y) -1000*(y - cos(t));
% y0 = 1;
% h = 1e-7;
% t = 0:h:0.1;
% 
% tic;
% y = backward_euler(f, y0, t, h);
% elapsed_time = toc;
% 
% fprintf('步长 h=%.1e, 步数 %d, 实际耗时 %.3f秒\n',...
%         h, length(t)-1, elapsed_time);
function y = backward_euler(f, y0, t, h)
    % 高精度低内存的向后欧拉实现
    % 输入：
    %   f: 函数句柄 f(t,y)
    %   y0: 初始值（标量或列向量）
    %   t: 时间向量（从t0开始，等间距）
    %   h: 步长（必须等于t(2)-t(1)）
    % 输出：
    %   y: 解向量（与t同维度）
    
    % 参数校验
    assert(abs(h - (t(2)-t(1)) )< eps, '步长h与时间向量不匹配');
    y0 = y0(:); % 确保列向量
    n = length(t);
    dim = length(y0);
    
    % 预分配输出（最大内存控制）
    max_memory = 1e9; % 1GB
    required_memory = n*dim*8; % 双精度每个元素8字节
    if required_memory > max_memory
        error('请求内存%.2fMB超过1GB限制', required_memory/1e6);
    end
    y = zeros(dim, n);
    y(:,1) = y0;
    
    % 针对线性系统优化)
    if isequal(f, @(t,y) -1000*(y - cos(t)))
        % 案例3的特殊优化
        for i = 1:n-1
            t_next = t(i) + h;
            y(:,i+1) = (y(:,i) + 1000*h*cos(t_next))/(1 + 1000*h);
        end
        return;
    end
    
    % 通用非线性系统实现
    tol = 1e-12; % 迭代容忍度
    max_iter = 100; % 最大迭代次数
    
    for i = 1:n-1
        t_next = t(i) + h;
        yn_guess = y(:,i); % 初始猜测
        
        % 牛顿迭代法
        for iter = 1:max_iter
            [F, J] = implicit_residual(f, yn_guess, y(:,i), t_next, h);
            
            % 解线性系统（避免矩阵求逆）
            delta = - (J\F);
            yn_guess = yn_guess + delta;
            
            if norm(delta) < tol
                break;
            end
        end
        y(:,i+1) = yn_guess;
    end
end

function [F, J] = implicit_residual(f, yn, yn_prev, t, h)
    % 计算隐式方程残差和雅可比
    F = yn - yn_prev - h*f(t, yn);
    
    % 数值计算雅可比（内存高效）
    epsilon = 1e-8;
    dim = length(yn);
    J = zeros(dim);
    for j = 1:dim
        yn_perturbed = yn;
        yn_perturbed(j) = yn_perturbed(j) + epsilon;
        F_perturbed = yn_perturbed - yn_prev - h*f(t, yn_perturbed);
        J(:,j) = (F_perturbed - F)/epsilon;
    end
    J = eye(dim) - J; % 调整为 ∂F/∂yn
end

function y_next = improved_euler(f, y, t, h)
    validate_scalar(y);
    y_predict = y + h*f(t, y);
    y_next = y + h/2*(f(t,y) + f(t+h, y_predict));
end

function y_next = adams_pc(f, y, t, h)
    persistent y_prev;
    if isempty(y_prev) || t == 0
        % 第一步用改进欧拉
        y_next = improved_euler(f, y, t, h);
        y_prev = y;
    else
        % AB2预测
        y_predict = y + h*(3/2*f(t,y) - 1/2*f(t-h, y_prev));
        % AM2校正
        y_next = y + h/2*(f(t,y) + f(t+h, y_predict));
        y_prev = y; % 更新历史值
    end
end

% ==== 向量方法 ====
function Y_next = backward_euler_vec(f, Y, t, h)
    validate_vector(Y);
    implicit_eq = @(Yn) Yn - Y - h*f(t+h, Yn);
    Y_next = fsolve(implicit_eq, Y + h*f(t,Y), optimset('Display','off'));
end

function Y_next = improved_euler_vec(f, Y, t, h)
    validate_vector(Y);
    Y_predict = Y + h*f(t, Y);
    Y_next = Y + h/2*(f(t,Y) + f(t+h, Y_predict));
end

function Y_next = adams_pc_vec(f, Y, t, h)
    persistent Y_prev;
    validate_vector(Y);
    if isempty(Y_prev) || t == 0
        % 第一步用改进欧拉
        Y_next = improved_euler_vec(f, Y, t, h);
        Y_prev = Y;
    else
        % AB2预测
        Y_predict = Y + h*(3/2*f(t,Y) - 1/2*f(t-h, Y_prev));
        % AM2校正
        Y_next = Y + h/2*(f(t,Y) + f(t+h, Y_predict));
        Y_prev = Y; % 更新历史值
    end
end

%% 辅助函数
function validate_scalar(y)
    if ~isscalar(y)
        error('输入必须是标量，当前维度：%s', mat2str(size(y)));
    end
end

function validate_vector(Y)
    if ~isvector(Y) || size(Y,2) ~= 1
        error('输入必须是列向量，当前维度：%s', mat2str(size(Y)));
    end
end

function analyze_errors(t, y_be, y_ie, y_adams, y_exact, case_num)
    error_be = abs(y_be - y_exact);
    error_ie = abs(y_ie - y_exact);
    error_adams = abs(y_adams - y_exact);
    
    fprintf('案例%d在t=%.2f处的误差：\n', case_num, t(end));
    fprintf('向后欧拉: %.4e\n', error_be(end));
    fprintf('改进欧拉: %.4e\n', error_ie(end));
    fprintf('Adams预测校正: %.4e\n\n', error_adams(end));
    
    figure('Name', sprintf('案例%d标量误差', case_num));
    plot(t, error_be, 'r--', t, error_ie, 'g-.', t, error_adams, 'b-');
    xlabel('t'); ylabel('误差'); 
    legend('向后欧拉','改进欧拉','Adams预测校正');
    title(sprintf('案例%d: 标量方程误差', case_num));
    grid on;
end

function analyze_vector_errors(t, Y_be, Y_ie, Y_adams, Y_exact, case_num)
    % 分析第一个分量（位移）
    error_be = abs(Y_be(1,:) - Y_exact(1,:));
    error_ie = abs(Y_ie(1,:) - Y_exact(1,:));
    error_adams = abs(Y_adams(1,:) - Y_exact(1,:));
    
    fprintf('案例%d在t=%.2f处的位移误差：\n', case_num, t(end));
    fprintf('向后欧拉: %.4e\n', error_be(end));
    fprintf('改进欧拉: %.4e\n', error_ie(end));
    fprintf('Adams预测校正: %.4e\n\n', error_adams(end));
    
    figure('Name', sprintf('案例%d向量误差', case_num));
    subplot(2,1,1);
    plot(t, error_be, 'r--', t, error_ie, 'g-.', t, error_adams, 'b-');
    title('位移分量误差'); xlabel('t'); ylabel('误差');
    
    % 分析第二个分量（速度）
    error_be_v = abs(Y_be(2,:) - Y_exact(2,:));
    error_ie_v = abs(Y_ie(2,:) - Y_exact(2,:));
    error_adams_v = abs(Y_adams(2,:) - Y_exact(2,:));
    
    subplot(2,1,2);
    plot(t, error_be_v, 'r--', t, error_ie_v, 'g-.', t, error_adams_v, 'b-');
    title('速度分量误差'); xlabel('t'); ylabel('误差');
    legend('向后欧拉','改进欧拉','Adams预测校正');
end