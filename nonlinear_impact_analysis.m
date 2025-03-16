function [t_nonlin, z_nonlin, t_lin, z_lin] = nonlinear_impact_analysis_v3(params, tspan, q0, dq0, options, varargin)
    %% 增强版非完全约束系统分析框架
    % 专注于稳定性和准确性，适用于复杂力学系统
    
    %% 输入参数处理 (增强验证)
    p = inputParser;
    addRequired(p, 'params', @(x) validate_params(x));
    addRequired(p, 'tspan', @(x) validateattributes(x, {'numeric'}, {'vector', 'increasing'}));
    addRequired(p, 'q0', @(x) validateattributes(x, {'numeric'}, {'vector'}));
    addRequired(p, 'dq0', @(x) validateattributes(x, {'numeric'}, {'vector', 'numel', numel(q0)}));
    addRequired(p, 'options', @(x) isempty(x) || isstruct(x) || isa(x, 'odeset'));
    addParameter(p, 'Potential', [], @(x) isa(x, 'function_handle'));
    addParameter(p, 'Constraints', [], @(x) isa(x, 'function_handle'));
    addParameter(p, 'ConstraintJacobian', [], @(x) isa(x, 'function_handle')); % 新增约束雅可比
    addParameter(p, 'Forces', [], @(x) isa(x, 'function_handle'));
    addParameter(p, 'MassMatrix', [], @(x) isa(x, 'function_handle'));
    addParameter(p, 'BaumgarteAlpha', 1.0, @(x) validateattributes(x, {'numeric'}, {'scalar', '>=', 0}));
    addParameter(p, 'BaumgarteBeta', 1.0, @(x) validateattributes(x, {'numeric'}, {'scalar', '>=', 0}));
    addParameter(p, 'SVDTol', 1e-8, @(x) validateattributes(x, {'numeric'}, {'scalar', '>', 0})); % SVD截断容差
    addParameter(p, 'ConstraintTol', 1e-6, @(x) validateattributes(x, {'numeric'}, {'scalar', '>', 0})); % 约束容差
    parse(p, params, tspan, q0, dq0, options, varargin{:});
    opts = p.Results;
    
    % 必要函数检查
    assert(~isempty(opts.Constraints), '必须提供约束函数');
    assert(~isempty(opts.Forces), '必须提供力函数');
    assert(~isempty(opts.MassMatrix), '必须提供质量矩阵函数');

    %% 系统维度验证
    n = numel(q0);  % 系统自由度
    [~, test_M] = opts.MassMatrix(0, q0, dq0, params);
    validate_dimensions(test_M, n, '质量矩阵');
    
    %% 求解器配置 (自适应策略)
    if isempty(options)
        options = odeset(...
            'RelTol', 1e-10, ... % 更高的相对容差
            'AbsTol', 1e-12, ... % 更高的绝对容差
            'Stats', 'on', ...
            'Mass', @(t,z) mass_matrix_wrapper(t,z,opts), ...
            'MStateDependence', 'strong', ...
            'MaxStep', diff(tspan)/100, ...
            'InitialStep', 1e-6); % 更小的初始步长
    end

    %% 主求解过程 (带故障恢复)
    try
        % 非线性系统求解
        [t_nonlin, z_nonlin] = ode15s(@(t,z) full_dynamics(t,z,opts), tspan, [q0; dq0], options);
        
        % 线性化系统求解
        [t_lin, z_lin] = ode15s(@(t,z) linearized_dynamics(t,z,opts), tspan, [q0; dq0], options);
    catch ME
%         handle_solver_error(ME, tspan, [q0; dq0], opts);
    end
    
    %% 结果后处理
    [t_nonlin, z_nonlin] = gather_results(t_nonlin, z_nonlin);
    [t_lin, z_lin] = gather_results(t_lin, z_lin);
    
    %% 系统分析
    if nargout == 0
        analyze_results(t_nonlin, z_nonlin, t_lin, z_lin, opts);
    end
end

%% 核心动力学函数 --------------------------------------------------------
function dz = full_dynamics(t, z, opts)
    [q, dq, n] = unpack_state(z);
    [M, A, F, b] = evaluate_system(t, q, dq, opts);
    
    % Baumgarte 约束稳定化
    if ~isempty(A)
        [phi, ~] = opts.Constraints(t, q, dq, opts.params);
        b = b - opts.BaumgarteAlpha*phi - opts.BaumgarteBeta*(A*dq);
    end
    
    % 鲁棒求解策略
    ddq = solve_constrained_dynamics(M, A, F, b, dq, n, opts);
    
    dz = pack_derivatives(dq, ddq);
end

function dz = linearized_dynamics(t, z, opts)
    [q, dq, n] = unpack_state(z);
    [M, A, F, b] = evaluate_system(t, q, dq, opts);
    
    % 数值线性化约束
    if ~isempty(A) && isempty(opts.ConstraintJacobian)
        A = numerical_jacobian(@(q) opts.Constraints(t, q, dq, opts.params), q);
    end
    
    ddq = solve_constrained_dynamics(M, A, F, b, dq, n, opts);
    dz = pack_derivatives(dq, ddq);
end

%% 辅助函数 ------------------------------------------------------------
function [q, dq, n] = unpack_state(z)
    n = length(z)/2;
    q = z(1:n);
    dq = z(n+1:end);
end

function dz = pack_derivatives(dq, ddq)
    dz = [dq; ddq];
end

function [M, A, F, b] = evaluate_system(t, q, dq, opts)
    % 计算系统各部分
    [M, F] = deal(opts.MassMatrix(t, q, dq, opts.params), ...
                 opts.Forces(t, q, dq, opts.params));
    
    % 约束计算
    if ~isempty(opts.Constraints)
        [phi, J] = opts.Constraints(t, q, dq, opts.params);
        A = J;
        b = -opts.ConstraintJacobian(t, q, dq, opts.params)*dq; % 假设提供约束雅可比
    else
        A = []; b = [];
    end
    
    % 维度一致性验证
    validate_dimensions(M, length(q), '质量矩阵');
    if ~isempty(A)
        validate_dimensions(A, length(q), '约束矩阵');
        validate_dimensions(b, size(A,1), '约束向量');
    end
end

function ddq = solve_constrained_dynamics(M, A, F, b, dq, n, opts)
    % 鲁棒约束处理策略
    if isempty(A)
        ddq = M \ F;
    else
        % 构建增广矩阵
        K = [M, A'; 
             A, zeros(size(A,1))];
        rhs = [F; b];
        
        % SVD求解 (处理秩亏情况)
        [U,S,V] = svd(K, 'econ');
        s = diag(S);
        valid_s = s > opts.SVDTol * max(s);
        inv_S = diag(1./s(valid_s));
        solution = V(:,valid_s) * inv_S * U(:,valid_s)' * rhs;
        
        % 验证解的质量
        residual = norm(K*solution - rhs);
        if residual > opts.ConstraintTol
            warning('约束求解残差过大: %.2e', residual);
        end
        
        ddq = solution(1:n);
    end
end

function J = numerical_jacobian(fun, q)
    % 自适应步长的数值雅可比
    epsilon = max(1e-8, 1e-6*norm(q));
    n = length(q);
    f0 = fun(q);
    m = length(f0);
    J = zeros(m, n);
    
    for k = 1:n
        q_plus = q; q_plus(k) = q_plus(k) + epsilon;
        q_minus = q; q_minus(k) = q_minus(k) - epsilon;
        J(:,k) = (fun(q_plus) - fun(q_minus)) / (2*epsilon);
    end
end

%% 验证函数 ------------------------------------------------------------
function validate_params(params)
    % 参数结构体验证
    required_fields = {'masses', 'lengths', 'stiffness'}; % 根据实际情况修改
    for fn = required_fields
        if ~isfield(params, fn{1})
            error('参数结构体缺少必要字段: %s', fn{1});
        end
    end
end

function validate_dimensions(mat, expected_dim, name)
    % 矩阵维度验证
    if isempty(mat), return; end
    if size(mat,1) ~= expected_dim || (ndims(mat)>=2 && size(mat,2) ~= expected_dim)
        error('%s 维度不匹配: 期望 %dx%d，实际 %dx%d', ...
              name, expected_dim, expected_dim, size(mat,1), size(mat,2));
    end
end

function handle_solver_error(ME, tspan, z0, opts)
    % 求解器错误处理
    fprintf('[错误] ODE求解失败: %s\n', ME.message);
    fprintf('尝试调整求解器参数...\n');
    
    % 自动调整策略
    new_options = odeset(opts.options, ...
        'RelTol', 1e-8, ...
        'AbsTol', 1e-10, ...
        'MaxStep', diff(tspan)/50);
    
    try
        % 非线性系统求解
        [t_nonlin, z_nonlin] = ode15s(@(t,z) full_dynamics(t,z,opts), tspan, z0, new_options);
        
        % 线性化系统求解
        [t_lin, z_lin] = ode15s(@(t,z) linearized_dynamics(t,z,opts), tspan, z0, new_options);
    catch ME2
        error('无法收敛: %s', ME2.message);
    end
end

%% 结果处理函数 --------------------------------------------------------
function [t, z] = gather_results(t, z)
    % 数据回传
    t = t;
    z = z;
end

function analyze_results(t_nonlin, z_nonlin, t_lin, z_lin, opts)
    % 结果分析可视化
    figure('Name', '动力学分析', 'Position', [100 100 1200 800])
    
    % 状态轨迹比较
    subplot(2,1,1)
    plot(t_nonlin, z_nonlin(:,1:end/2)), hold on
    plot(t_lin, z_lin(:,1:end/2), '--')
    title('广义坐标对比'), xlabel('时间'), ylabel('状态')
    
    % 能量分析
    subplot(2,1,2)
    plot_energy(t_nonlin, z_nonlin, opts)
    title('系统能量'), xlabel('时间'), ylabel('能量')
end

function plot_energy(t, z, opts)
    % 能量计算与绘图
    n = size(z,2)/2;
    KE = zeros(size(t));
    PE = zeros(size(t));
    
    for i = 1:length(t)
        [q, dq] = unpack_state(z(i,:)');
        M = opts.MassMatrix(t(i), q, dq, opts.params);
        KE(i) = 0.5 * dq' * M * dq;
        if ~isempty(opts.Potential)
            PE(i) = opts.Potential(t(i), q, opts.params);
        end
    end
    
    plot(t, KE, 'b', t, PE, 'r', t, KE+PE, 'g--')
    legend('动能', '势能', '总能量')
end