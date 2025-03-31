%% ODE数据生成脚本
clear; clc; close all;
format longG;
rng(42, 'twister');

% ==================== 配置部分 ====================
config.DoublePrecision = true;
config.GPUMemThreshold = 0.8;
config.TestDimensions = [2, 5, 10];
config.TestSteps = [1e4, 5e4];

% ==================== 主执行部分 ====================
try
    % 1. GPU初始化
    [gpu_available, gpu_info] = init_gpu_r2022a(config);
    if ~gpu_available
        error('终止执行：GPU不可用');
    end
    fprintf('>> [状态] %s (%.1fGB显存)\n', gpu_info.Name, gpu_info.TotalMemory/1e9);

    % 2. 并行池管理（修正线程池创建方式）
    p = gcp('nocreate');
    if isempty(p)
        % 关键修改：移除workers数量参数
        p = parpool('threads'); 
    end
    fprintf('>> 并行池启动 (线程池模式)\n');
    
    % 3. 扩展测试系统定义（新增8种类型）
test_systems = {
    % 基础类型
    '线性一阶',    @linear_1st_order,   false, false, '1st', false, false;
    '非线性一阶',  @nonlinear_1st_order, false, true,  '1st', false, true;
    '刚性系统',    @stiff_system,        true,  false, '1st', false, false;
    '延迟系统',    @delay_system,       false, true,  '1st', true,  true;
    '二阶线性',    @linear_2nd_order,   false, false, '2nd', false, false;
    '二阶非线性',  @nonlinear_2nd_order, false, true,  '2nd', false, true;
    
    % 新增类型（机器学习关键扩展）
    '周期系统',    @periodic_system,    false, true,  '1st', false, true;  % 新增1
    '混沌系统',    @lorenz_system,      false, true,  '1st', false, true;  % 新增2  
    '脉冲系统',    @impulse_system,     false, true,  '1st', true,  false; % 新增3
    '分段系统',    @piecewise_system,   false, true,  '1st', false, true;  % 新增4
    '高振荡系统',  @high_osc_system,    true,  false, '1st', false, false; % 新增5
    '参数时变系统',@time_varying_system,false, false, '1st', false, true;  % 新增6
    '守恒系统',    @conservative_system,false, true,  '2nd', false, false; % 新增7
    '随机驱动系统',@stochastic_system,  false, true,  '1st', false, true   % 新增8
};

    % 4. 主计算循环
    results = table();
    for sys_idx = 1:size(test_systems,1)
        sys_data = test_systems(sys_idx, :);
        parfor (dim_idx = 1:length(config.TestDimensions), 0) % 使用默认并行设置
            d = config.TestDimensions(dim_idx);
            for steps = config.TestSteps
                [result, valid] = run_enhanced_test(sys_data, d, steps, config, gpu_info);
                if valid, results = [results; result]; end
            end
        end
    end

    % 5. 结果处理
    if ~isempty(results)
        results.Properties.VariableNames = {
            'Dimension','Steps','System','ElapsedTime',...
            'IsStiff','IsNonlinear','Order','IsDelay','HasParams'
        };
        save('ode_timing_data.mat', 'results');
        disp('==== 增强测试结果 ====');
        disp(sortrows(results, 'ElapsedTime', 'descend'));
    end

catch ME
    fprintf('\n!! 错误: %s\n', ME.message);
end

% 确保并行池关闭
delete(gcp('nocreate'));
fprintf('>> 计算完成，并行池已释放\n');

% ==================== 核心函数 ====================
function [result, valid] = run_enhanced_test(sys_data, base_dim, steps, config, gpu_info)
    sys_name = sys_data{1};
    ode_fun = sys_data{2};
    is_stiff = sys_data{3};
    is_nonlin = sys_data{4};
    order = sys_data{5};
    is_delay = sys_data{6};
    has_params = sys_data{7};
    
    % 维度转换（二阶系统）
    actual_dim = base_dim * (1 + strcmp(order,'2nd'));
    
    result = table();
    valid = false;
    
    fprintf('计算 %s dim=%d steps=%d...\n', sys_name, actual_dim, steps);
    
    try
        % 显存检查
        mem_need = 32 * actual_dim^2 * steps / 1e6;
        if mem_need > config.GPUMemThreshold * gpu_info.TotalMemory/1e6
            fprintf('  跳过：需要%.1fMB显存\n', mem_need);
            return;
        end
        
        % 维度预验证
        y_test = gpuArray(rand(actual_dim,1,'double'));
        if is_delay
            Z_test = zeros(actual_dim, 1, 'double');
            test_output = ode_fun(0, y_test, Z_test);
        else
            test_output = ode_fun(0, y_test);
        end
        if numel(test_output) ~= actual_dim
            error('维度不匹配: 输出应为%d×1',actual_dim);
        end
        
        % 执行计算
        tic;
        if is_delay
            lags = 1; % 延迟时间
            sol = dde23(ode_fun, lags, @(t) zeros(actual_dim,1), [0, 5]);
        else
            [~,~] = ode45(@(t,y) ode_wrapper(t,y,ode_fun), ...
                     linspace(0,5,steps), rand(actual_dim,1));
        end
        elapsed = toc;
        
        % 构建结果
        result = table(...
            actual_dim, steps, {sys_name}, elapsed,...
            is_stiff, is_nonlin, {order}, is_delay, has_params);
        valid = true;
        fprintf('  成功: %.3fs\n', elapsed);
        
    catch ME
        fprintf('!! %s失败(dim=%d): %s\n', sys_name, actual_dim, ME.message);
        reset(gpuDevice);
    end
    
    function dy = ode_wrapper(t, y, fun)
        y_gpu = gpuArray(colvec(y,'double'));
        if is_delay
            Z_test = zeros(size(y_gpu,1), 1, 'double');
            dy_gpu = fun(t, y_gpu, Z_test);
        else
            dy_gpu = fun(t, y_gpu);
        end
        dy = gather(dy_gpu);
    end
end

% ==================== 增强系统函数 ====================
function dy = linear_1st_order(~, y)
    d = numel(y);
    A = gpuArray.randn(d,d,'double')/(2*sqrt(double(d)));
    A = A - diag(1.1*abs(diag(A))) - 0.1*eye(d,'double');
    dy = A * y;
end

function dy = nonlinear_1st_order(~, y)
    d = numel(y);
    norm_factor = sum(y.^2)/double(d);
    dy = y .* (1 - norm_factor);
end

function dy = stiff_system(~, y)
    d = numel(y);
    main_diag = -100 * ones(d,1);
    off_diag = 0.1 * randn(d, d);
    off_diag = off_diag - diag(diag(off_diag));
    J = diag(main_diag) + off_diag;
    dy = J * y;
end

function dy = delay_system(t, y, Z)
    ylag = Z(:,1);
    dy = 0.5*ylag - y.^3./(1+y.^2);
end

function dy = linear_2nd_order(~, y)
    n = length(y)/2;
    dy = zeros(2*n,1);
    dy(1:n) = y(n+1:end);
    dy(n+1:end) = -0.1*y(n+1:end) - 0.5*y(1:n);
end

function dy = nonlinear_2nd_order(~, y)
    n = length(y)/2;
    dy = zeros(2*n,1);
    dy(1:n) = y(n+1:end);
    dy(n+1:end) = -0.1*y(n+1:end) - y(1:n).^3;
end

% ==================== 新增系统函数 ====================
function dy = periodic_system(~, y)
    d = numel(y);
    dy = zeros(d, 1);
    for i = 1:2:d-1
        dy(i) = y(i+1);
        dy(i+1) = -0.1*y(i+1) + sin(y(i));
    end
    if mod(d,2) == 1
        dy(end) = sin(y(end)) - 0.1*y(end);
    end
end

function dy = lorenz_system(~, y)
    d = numel(y);
    dy = zeros(d,1);
    if d >= 3
        sigma = 10; rho = 28; beta = 8/3;
        dy(1) = sigma*(y(2)-y(1));
        dy(2) = y(1)*(rho-y(3))-y(2);
        dy(3) = y(1)*y(2)-beta*y(3);
    end
    if d > 3
        dy(4:end) = -0.1 * y(4:end);
    end
end

function dy = impulse_system(t, y, Z)
    ylag = Z(:,1);
    dy = -0.5*y + 1.2*ylag.*(mod(floor(t),2)==0);
end

function dy = piecewise_system(~, y)
    d = numel(y);
    if y(1) > 0
        A = -0.2*eye(d) + diag(ones(d-1,1),1) - diag(ones(d-1,1),-1);
    else
        A = -0.2*eye(d) - diag(ones(d-1,1),1) + diag(ones(d-1,1),-1);
    end
    dy = A * y;
end

function dy = high_osc_system(~, y)
    d = numel(y);
    dy = zeros(d,1);
    omega = 100;
    for i = 1:2:d-1
        dy(i) = y(i+1);
        dy(i+1) = -omega^2*y(i) - 0.01*y(i+1);
    end
    if mod(d,2) == 1
        dy(end) = -0.01*y(end);
    end
end

function dy = time_varying_system(t, y)
    d = numel(y);
    dy = zeros(d,1);
    a = 0.1 + 0.05*sin(t);
    b = 1 + 0.5*cos(t/2);
    for i = 1:2:d-1
        dy(i) = y(i+1);
        dy(i+1) = -a*y(i+1) - b*y(i);
    end
    if mod(d,2) ==1
        dy(end) = -a*y(end);
    end
end

function dy = conservative_system(~, y)
    n = length(y)/2;
    q = y(1:n); p = y(n+1:end);
    dq = p;
    dp = -q./(sqrt(sum(q.^2))^3);
    dy = [dq; dp];
end

function dy = stochastic_system(t, y)
    persistent w
    d = numel(y);
    if isempty(w)
        w = 0.1*randn(d,1);
    end
    dy = zeros(d,1);
    for i = 1:2:d-1
        dy(i) = y(i+1);
        dy(i+1) = -0.1*y(i+1) - y(i)^3 + 0.05*sin(t);
    end
    if mod(d,2) ==1
        dy(end) = -0.1*y(end) + 0.05*sin(t);
    end
    dy = dy + w;
end

% ==================== 辅助函数 ====================
function [success, info] = init_gpu_r2022a(~)
    success = false;
    info = struct();
    try
        if gpuDeviceCount == 0, error('无GPU设备'); end
        device = gpuDevice(1);
        reset(device);
        
        test = gpuArray.rand(10,'double') * gpuArray.rand(10,'double')';
        if ~isa(test,'gpuArray') || ~strcmp(classUnderlying(test),'double')
            error('双精度测试失败');
        end
        
        success = true;
        info = struct('Name',device.Name, 'TotalMemory',device.TotalMemory, ...
                     'ComputeCapability',device.ComputeCapability);
    catch ME
        fprintf('GPU初始化失败: %s\n', ME.message);
    end
end

function v = colvec(v, precision)
    if isrow(v), v = v'; end
    if nargin>1 && ~strcmp(class(v),precision)
        v = cast(v,precision);
    end
end