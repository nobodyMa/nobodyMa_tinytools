%经典四阶龙格-库塔法（RK4）

function [t, y] = rk4(f, tspan, y0, h)
    % 输入参数：
    %   f: 函数句柄，定义ODE
    %   tspan: 时间区间向量，格式为 [起始时间, 终止时间]
    %   y0: 初始条件，列向量（如 [1; 2] 表示多变量系统）
    %   h: 固定步长（标量）
    % 输出：
    %   t: 时间点列向量（格式为 N×1）
    %   y: 解矩阵，格式为 N×M，N为时间点数，M为变量维度

    t = tspan(1):h:tspan(end);          % 生成离散时间点向量，语法：start:step:end
    if t(end) ~= tspan(end)             % 检查是否覆盖终止时间
        t = [t, tspan(end)];           % 若不覆盖，添加终止时间点
    end
    n = length(t);                      % 时间点总数（语法：length(向量) 返回元素个数）
    m = length(y0);                     % 变量维度（如单变量为1，多变量为2+）
    y = zeros(n, m);                    % 预分配解矩阵（语法：zeros(行数, 列数)）
    y(1,:) = y0';                       % 初始条件赋值，转置为行向量存入第一行

    for i = 1:n-1                       % 时间循环，从第一个到倒数第二个时间点
        tn = t(i);                      % 当前时间
        yn = y(i,:)';                   % 当前状态列向量（从行转置为列）
        
        % --- 计算四个斜率 ---
        k1 = f(tn, yn);                 % 斜率k1，基于当前点(tn, yn)
        k2 = f(tn + h/2, yn + h*k1/2);  % 斜率k2，基于中点(tn + h/2, yn + k1*h/2)
        k3 = f(tn + h/2, yn + h*k2/2);  % 斜率k3，基于修正后的中点
        k4 = f(tn + h, yn + h*k3);      % 斜率k4，基于下一个时间点(tn + h, yn + k3*h)
        
        % --- 更新下一步的解 ---
        y(i+1,:) = yn' + h*(k1 + 2*k2 + 2*k3 + k4)/6; % 加权平均，转置为行向量存储
    end
    t = t';                             % 将时间向量转为列向量（语法：' 表示转置）
end

%显式Adams-Bashforth四阶方法

function [t, y] = adams_bashforth_4(f, tspan, y0, h)
    % 输入参数同RK4
    t = tspan(1):h:tspan(end);
    if t(end) ~= tspan(end)
        t = [t, tspan(end)];
    end
    n = length(t);
    m = length(y0);
    y = zeros(n, m);
    
    % --- 用RK4初始化前4步 ---
    [t_init, y_init] = rk4(f, [t(1), t(4)], y0, h); % 调用RK4生成前4步
    y(1:4,:) = y_init(1:4,:);                       % 将前4步结果存入y
    
    % --- Adams-Bashforth主循环 ---
    for i = 4:n-1
        tn = t(i);
        yn = y(i,:)';                   % 当前状态列向量
        
        % 提取前四步的斜率（历史数据）
        f_n   = f(tn, yn);              % 当前步斜率
        f_n_1 = f(t(i-1), y(i-1,:)');   % 前1步斜率
        f_n_2 = f(t(i-2), y(i-2,:)');   % 前2步斜率
        f_n_3 = f(t(i-3), y(i-3,:)');   % 前3步斜率
        
        % Adams-Bashforth四阶公式
        y_next = yn + h*(55*f_n - 59*f_n_1 + 37*f_n_2 - 9*f_n_3)/24;
        y(i+1,:) = y_next';             % 转置为行向量存储
    end
    t = t';                             % 时间向量转为列向量
end

%隐式梯形法（Crank-Nicolson）

function [t, y] = implicit_trapezoid(f, tspan, y0, h, tol, max_iter)
    % 新增输入参数：
    %   tol: 牛顿迭代收敛容差（默认1e-6）
    %   max_iter: 最大迭代次数（默认100）
    if nargin < 6, max_iter = 100; end    % 检查输入参数个数（语法：nargin）
    if nargin < 5, tol = 1e-6; end
    
    t = tspan(1):h:tspan(end);
    if t(end) ~= tspan(end)
        t = [t, tspan(end)];
    end
    n = length(t);
    m = length(y0);
    y = zeros(n, m);
    y(1,:) = y0';
    
    for i = 1:n-1
        tn = t(i);
        yn = y(i,:)';                   % 当前状态列向量
        tn1 = t(i+1);                   % 下一步时间
        
        % --- 牛顿迭代求解隐式方程 ---
        y_guess = yn;                   % 初始猜测为当前值
        for iter = 1:max_iter
            % 计算残差 F = y_guess - yn - h/2*(f(tn,yn) + f(tn1,y_guess))
            F = y_guess - yn - h/2*(f(tn, yn) + f(tn1, y_guess));
            
            % 计算雅可比矩阵 J = I - h/2 * df/dy(tn1, y_guess)
            J = eye(m) - h/2 * jacobian(f, tn1, y_guess); % 需自定义jacobian函数
            
            % 求解线性系统 J * delta = -F
            delta = -J \ F;             % 语法：反斜杠\ 表示解线性方程
            
            % 更新猜测值
            y_guess = y_guess + delta;
            
            % 检查收敛条件
            if norm(delta) < tol
                break;
            end
        end
        
        y(i+1,:) = y_guess';            % 存储结果
    end
    t = t';
end

% --- 自定义雅可比矩阵计算函数 ---
function J = jacobian(f, t, y)
    % 数值计算雅可比矩阵 df/dy
    epsilon = 1e-8;                     % 扰动步长
    m = length(y);
    J = zeros(m, m);
    f0 = f(t, y);                       % 原始函数值
    
    for k = 1:m
        y_perturbed = y;
        y_perturbed(k) = y_perturbed(k) + epsilon;  % 对第k个变量加扰动
        f_perturbed = f(t, y_perturbed);
        J(:,k) = (f_perturbed - f0) / epsilon;      % 有限差分近似导数
    end
end

%辛积分方法（Verlet算法）

function [t, y, v] = verlet(f, tspan, y0, v0, h)
    % 输入参数：
    %   f: 加速度函数句柄，a = f(t, y)，接受标量t和位置y，返回标量加速度
    %   y0: 初始位置（标量）
    %   v0: 初始速度（标量）
    % 输出：
    %   y: 位置向量（N×1）
    %   v: 速度向量（N×1）
    
    t = tspan(1):h:tspan(end);
    if t(end) ~= tspan(end)
        t = [t, tspan(end)];
    end
    n = length(t);
    y = zeros(n, 1);
    v = zeros(n, 1);
    y(1) = y0;
    v(1) = v0;
    
    % --- 第一步：用速度Verlet初始化 ---
    a_prev = f(t(1), y(1));             % 初始加速度
    y(2) = y(1) + h*v(1) + 0.5*h^2*a_prev;
    
    % --- 后续步骤：位置Verlet ---
    for i = 2:n-1
        a_current = f(t(i), y(i));      % 当前加速度
        y(i+1) = 2*y(i) - y(i-1) + h^2*a_current;  % Verlet位置更新公式
        
        % 计算速度（中心差分）
        v(i) = (y(i+1) - y(i-1)) / (2*h);          % 语法：向量索引从1开始
    end
    
    % --- 补全最后一个速度 ---
    v(end) = (y(end) - y(end-1)) / h;  % 前向差分近似
end

