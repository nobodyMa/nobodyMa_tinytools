vibration_visualizer_enhanced

%%
function vibration_visualizer_enhanced
    % 创建主界面
    fig = figure('Name','振动现象ODE可视化工具', 'Position',[100 100 1000 650], 'NumberTitle','off',...
        'MenuBar','none', 'ToolBar','none', 'Color',[0.95 0.95 0.95]);
    
    % 振动类型选择 (顶部居中)
    uicontrol('Style','text', 'Position',[350 610 300 20], 'String','选择振动类型:',...
        'FontSize',10, 'BackgroundColor',[0.95 0.95 0.95]);
    vibration_type = uicontrol('Style','popup', 'Position',[350 580 300 30],...
        'String',{'自由振动(无阻尼)','自由振动(有阻尼)','强迫振动(正弦外力)','强迫振动(动力型外力)','共振现象','静态外力作用','自定义外力函数'},...
        'FontSize',10, 'Callback',@update_parameters, 'BackgroundColor','white');
    
    % 参数输入区域 (左侧面板)
    param_panel = uipanel('Title','系统参数', 'Position',[0.02 0.4 0.35 0.25],...
        'FontSize',11, 'BackgroundColor',[0.95 0.95 0.95], 'BorderType','etchedin');
    
    % 质量参数
    uicontrol('Parent',param_panel, 'Style','text', 'Position',[10 90 80 20],...
        'String','质量 m:', 'FontSize',10, 'HorizontalAlignment','left',...
        'BackgroundColor',[0.95 0.95 0.95]);
    m_edit = uicontrol('Parent',param_panel, 'Style','edit', 'Position',[100 90 100 25],...
        'String','1', 'FontSize',10, 'Tag','m', 'BackgroundColor','white');
    
    % 刚度参数
    uicontrol('Parent',param_panel, 'Style','text', 'Position',[10 60 80 20],...
        'String','刚度 k:', 'FontSize',10, 'HorizontalAlignment','left',...
        'BackgroundColor',[0.95 0.95 0.95]);
    k_edit = uicontrol('Parent',param_panel, 'Style','edit', 'Position',[100 60 100 25],...
        'String','1', 'FontSize',10, 'Tag','k', 'BackgroundColor','white');
    
    % 阻尼参数 (默认隐藏)
    c_label = uicontrol('Parent',param_panel, 'Style','text', 'Position',[10 30 80 20],...
        'String','阻尼 c:', 'FontSize',10, 'HorizontalAlignment','left',...
        'Visible','off', 'Tag','c_label', 'BackgroundColor',[0.95 0.95 0.95]);
    c_edit = uicontrol('Parent',param_panel, 'Style','edit', 'Position',[100 30 100 25],...
        'String','0.1', 'FontSize',10, 'Tag','c', 'Visible','off', 'BackgroundColor','white');
    
    % 初值条件面板 (左侧中部)
    iv_panel = uipanel('Title','初值条件', 'Position',[0.02 0.2 0.35 0.18],...
        'FontSize',11, 'BackgroundColor',[0.95 0.95 0.95], 'BorderType','etchedin');
    
    uicontrol('Parent',iv_panel, 'Style','text', 'Position',[10 60 80 20],...
        'String','初始位移 x0:', 'FontSize',10, 'HorizontalAlignment','left',...
        'BackgroundColor',[0.95 0.95 0.95]);
    x0_edit = uicontrol('Parent',iv_panel, 'Style','edit', 'Position',[100 60 100 25],...
        'String','1', 'FontSize',10, 'Tag','x0', 'BackgroundColor','white');
    
    uicontrol('Parent',iv_panel, 'Style','text', 'Position',[10 30 80 20],...
        'String','初始速度 v0:', 'FontSize',10, 'HorizontalAlignment','left',...
        'BackgroundColor',[0.95 0.95 0.95]);
    v0_edit = uicontrol('Parent',iv_panel, 'Style','edit', 'Position',[100 30 100 25],...
        'String','0', 'FontSize',10, 'Tag','v0', 'BackgroundColor','white');
    
    % 外力参数面板 (默认隐藏)
    force_panel = uipanel('Title','外力参数', 'Position',[0.02 0.02 0.35 0.16],...
        'FontSize',11, 'Visible','off', 'Tag','force_panel',...
        'BackgroundColor',[0.95 0.95 0.95], 'BorderType','etchedin');
    
    % 外力幅值
    uicontrol('Parent',force_panel, 'Style','text', 'Position',[10 90 80 20],...
        'String','幅值 F0:', 'FontSize',10, 'HorizontalAlignment','left',...
        'BackgroundColor',[0.95 0.95 0.95]);
    F0_edit = uicontrol('Parent',force_panel, 'Style','edit', 'Position',[100 90 100 25],...
        'String','1', 'FontSize',10, 'Tag','F0', 'BackgroundColor','white');
    
    % 外力频率
    uicontrol('Parent',force_panel, 'Style','text', 'Position',[10 60 80 20],...
        'String','频率 ω:', 'FontSize',10, 'HorizontalAlignment','left',...
        'BackgroundColor',[0.95 0.95 0.95]);
    omega_edit = uicontrol('Parent',force_panel, 'Style','edit', 'Position',[100 60 100 25],...
        'String','1', 'FontSize',10, 'Tag','omega', 'BackgroundColor','white');
    
    % 外力类型选择 (用于动力型外力)
    force_type_label = uicontrol('Parent',force_panel, 'Style','text', 'Position',[10 30 80 20],...
        'String','外力类型:', 'FontSize',10, 'HorizontalAlignment','left',...
        'Visible','off', 'Tag','force_type_label', 'BackgroundColor',[0.95 0.95 0.95]);
    force_type = uicontrol('Parent',force_panel, 'Style','popup', 'Position',[100 30 100 25],...
        'String',{'正弦','阶跃','脉冲','随机'}, 'FontSize',10, 'Visible','off',...
        'Tag','force_type', 'BackgroundColor','white');
    
    % 自定义外力函数面板 (默认隐藏)
    custom_force_panel = uipanel('Title','自定义外力函数', 'Position',[0.02 0.02 0.35 0.16],...
        'FontSize',11, 'Visible','off', 'Tag','custom_force_panel',...
        'BackgroundColor',[0.95 0.95 0.95], 'BorderType','etchedin');
    
    uicontrol('Parent',custom_force_panel, 'Style','text', 'Position',[10 60 280 20],...
        'String','输入外力函数 F(t,x,x'',x"):', 'FontSize',10, 'HorizontalAlignment','left',...
        'BackgroundColor',[0.95 0.95 0.95]);
    custom_force_edit = uicontrol('Parent',custom_force_panel, 'Style','edit', 'Position',[10 30 280 30],...
        'String','sin(t) + 0.5*x''', 'FontSize',10, 'Tag','custom_force', 'BackgroundColor','white',...
        'Tooltip','可以使用t(时间)、x(位移)、xd(速度)、xdd(加速度)变量，例如: sin(t) + 0.5*xd - 0.1*x');
    
    % 时间参数 (右侧面板)
    time_panel = uipanel('Title','时间参数', 'Position',[0.02 0.02 0.35 0.16],...
        'FontSize',11, 'BackgroundColor',[0.95 0.95 0.95], 'BorderType','etchedin');
    
    uicontrol('Parent',time_panel, 'Style','text', 'Position',[10 60 80 20],...
        'String','时间范围:', 'FontSize',10, 'HorizontalAlignment','left',...
        'BackgroundColor',[0.95 0.95 0.95]);
    t_edit = uicontrol('Parent',time_panel, 'Style','edit', 'Position',[100 60 100 25],...
        'String','0:0.01:20', 'FontSize',10, 'Tag','t', 'BackgroundColor','white');
    
    uicontrol('Parent',time_panel, 'Style','text', 'Position',[10 30 80 20],...
        'String','求解器:', 'FontSize',10, 'HorizontalAlignment','left',...
        'BackgroundColor',[0.95 0.95 0.95]);
    solver_popup = uicontrol('Parent',time_panel, 'Style','popup', 'Position',[100 30 100 25],...
        'String',{'自动选择','ode45','ode23','ode113','ode15s','ode23s'}, 'FontSize',10,...
        'BackgroundColor','white', 'Value',1);
    
    % 绘图按钮 (底部居中)
    plot_btn = uicontrol('Style','pushbutton', 'Position',[400 20 120 40],...
        'String','绘制振动曲线', 'FontSize',11, 'Callback',@plot_vibration,...
        'BackgroundColor',[0.4 0.8 0.6], 'ForegroundColor','white',...
        'FontWeight','bold');
    
    % 状态提示文本 (底部)
    status_text = uicontrol('Style','text', 'Position',[20 5 960 20],...
        'String','准备就绪', 'FontSize',10, 'ForegroundColor',[0.2 0.5 0.2],...
        'BackgroundColor',[0.95 0.95 0.95], 'HorizontalAlignment','left');
    
    % 创建坐标轴但初始隐藏
ax = axes('Parent',fig, 'Position',[0.4 0.1 0.58 0.8], 'Tag','main_ax',...
    'Box','on', 'FontSize',10, 'Color','white', 'Visible','off');
grid(ax, 'on');

    % 操作指南面板
guide_panel = uipanel('Parent',fig, 'Title','操作指南', 'Position',[0.4 0.1 0.58 0.8],...
    'FontSize',11, 'BackgroundColor',[0.96 0.96 0.96], 'BorderType','etchedin');

% ===== 完整版操作指南（含求解器指南） =====
% 创建可滚动文本区域 (上半部分)
guide_text = {
    '【核心功能】'
    '1. 自由振动(无/有阻尼)  2. 强迫振动(正弦/阶跃/脉冲/随机)'
    '3. 共振分析  4. 静态外力  5. 自定义非线性外力'
    ''
    '【输入规范】'
    '• 基本参数：m>0, k>0, c≥0 (单位: kg, N/m, N·s/m)'
    '• 初值条件：x0(m), v0(m/s)'
    '• 时间范围：起始:步长:结束 或 linspace(起始,结束,点数)'
    '• 自定义外力：可用t(时间),x(位移),xd(速度),xdd(加速度)变量'
    '  示例：5*exp(-t)*sin(xd)'};

uicontrol('Parent',guide_panel, 'Style','edit', 'Position',[10 120 560 300],...
    'String',guide_text, 'FontSize',9, 'HorizontalAlignment','left',...
    'BackgroundColor','white', 'Max',2, 'Enable','inactive');

% 创建表格数据
solver_data = {
    'ode45'  '非刚性系统(默认选择)'    '无阻尼自由振动'
    'ode23'  '轻度刚性/快速计算'       '初步结果预览'  
    'ode113' '高精度需求'            '长时间仿真'
    'ode15s' '刚性系统/不连续外力'     '大阻尼/脉冲外力'
    'ode23s' '强刚性系统'            '极高刚度(m/k>1e6)'
    '自动'   '智能识别'              '所有常规情况'};

% 创建表格 (下半部分)
uitable('Parent',guide_panel, 'Position',[10 10 560 150],...
    'Data', solver_data,...
    'ColumnName', {'求解器', '适用场景', '典型案例'},...
    'ColumnWidth', {60 180 180},...
    'FontSize', 8,...
    'RowName', [],...
    'BackgroundColor', [1 1 1; 0.95 0.95 0.95]);

% 添加底部提示
uicontrol('Parent',guide_panel, 'Style','text', 'Position',[10 0 560 20],...
    'String','【提示】右键图形可导出数据/图片 | 表达式支持所有MATLAB数学函数',...
    'FontSize',8, 'HorizontalAlignment','left',...
    'BackgroundColor',[0.98 0.98 0.98]);
    
    % 参数显示
    function update_parameters(~,~)
    type = vibration_type.Value;
    
    % 隐藏所有外力相关面板
    force_panel.Visible = 'off';
    custom_force_panel.Visible = 'off';
    time_panel.Position = [0.02 0.02 0.35 0.16]; % 默认位置
    
    % 根据振动类型显示/隐藏相关参数
    switch type
        case 1 % 自由振动(无阻尼)
            c_label.Visible = 'off';
            c_edit.Visible = 'off';
            time_panel.Position = [0.02 0.02 0.35 0.16]; % 正常位置
            
        case 2 % 自由振动(有阻尼)
            c_label.Visible = 'on';
            c_edit.Visible = 'on';
            time_panel.Position = [0.02 0.02 0.35 0.16]; % 正常位置
            
        case {3, 5} % 强迫振动或共振
            c_label.Visible = 'on';
            c_edit.Visible = 'on';
            force_panel.Visible = 'on';
            force_type_label.Visible = 'off';
            force_type.Visible = 'off';
            time_panel.Position = [0.02 0.02 0.35 0.16]; % 保持原位置
            
            if type == 4 % 动力型外力
                force_type_label.Visible = 'on';
                force_type.Visible = 'on';
            elseif type == 5 % 共振
                % 自动设置外力频率为固有频率
                k = str2double(get(k_edit, 'String'));
                m = str2double(get(m_edit, 'String'));
                if m <= 0
                    set(status_text, 'String', '错误: 质量必须大于0');
                    return;
                end
                omega_n = sqrt(k/m); % 固有频率
                set(omega_edit, 'String', num2str(omega_n));
            end
            
        case 4 % 动力型外力
            c_label.Visible = 'on';
            c_edit.Visible = 'on';
            force_panel.Visible = 'on';
            force_type_label.Visible = 'on';
            force_type.Visible = 'on';
            time_panel.Position = [0.02 0.02 0.35 0.16]; % 保持原位置
            
        case 6 % 静态外力作用
            c_label.Visible = 'on';
            c_edit.Visible = 'on';
            force_panel.Visible = 'on';
            set(omega_edit, 'String', '0'); % 静态外力频率为0
            time_panel.Position = [0.02 0.02 0.35 0.16]; % 保持原位置
            
        case 7 % 自定义外力函数
            c_label.Visible = 'on';
            c_edit.Visible = 'on';
            custom_force_panel.Visible = 'on';
            time_panel.Position = [0.02 0.19 0.35 0.16]; % 上移时间面板
    end
    set(status_text, 'String', '参数已更新，点击"绘制振动曲线"进行计算');
end

    % 主绘图函数
    function plot_vibration(~,~)
    cla(ax);
    set(status_text, 'String', '计算中，请稍候...', 'ForegroundColor',[0.2 0.5 0.2]);
    drawnow;
        
        try
            % 获取参数并进行验证
            m = str2double(get(m_edit, 'String'));
            k = str2double(get(k_edit, 'String'));
            x0 = str2double(get(x0_edit, 'String'));
            v0 = str2double(get(v0_edit, 'String'));
            
            if isnan(m) || m <= 0
                error('质量m必须是正数');
            end
            if isnan(k) || k <= 0
                error('刚度k必须是正数');
            end
            if isnan(x0)
                error('初始位移x0必须是数值');
            end
            if isnan(v0)
                error('初始速度v0必须是数值');
            end
            
            % 解析时间范围
            try
                t_span = eval(get(t_edit, 'String'));
                if numel(t_span) < 2
                    error('时间范围至少需要2个点');
                end
                if any(diff(t_span) <= 0)
                    error('时间点必须严格递增');
                end
            catch ME
                error('时间范围格式错误，示例: 0:0.01:10 或 linspace(0,10,1000)');
            end
            
            % 根据振动类型获取其他参数
            type = vibration_type.Value;
            if type ~= 1 % 不是无阻尼自由振动
                c = str2double(get(c_edit, 'String'));
                if isnan(c) || c < 0
                    error('阻尼c必须是非负数');
                end
            else
                c = 0;
            end
            
            % 计算固有频率和阻尼比
            omega_n = sqrt(k/m); % 固有频率
            zeta = c/(2*sqrt(m*k)); % 阻尼比
            
            % 获取求解器选择
            solver_choice = solver_popup.Value;
            solver_list = {'auto','ode45','ode23','ode113','ode15s','ode23s'};
            selected_solver = solver_list{solver_choice};
            
            % 根据振动类型求解
            switch type
                case {1, 2} % 自由振动
                    [t, x, v] = solve_free_vibration(m, c, k, x0, v0, t_span, selected_solver);
                    plot_results(t, x, v, '自由振动响应', type);
                    
                case {3, 5} % 强迫振动或共振
                    F0 = str2double(get(F0_edit, 'String'));
                    omega = str2double(get(omega_edit, 'String'));
                    
                    if isnan(F0)
                        error('外力幅值F0必须是数值');
                    end
                    if isnan(omega) || omega < 0
                        error('外力频率ω必须是非负数');
                    end
                    
                    [t, x, v] = solve_forced_vibration(m, c, k, F0, omega, x0, v0, t_span, selected_solver);
                    
                    if type == 5
                        title_str = sprintf('共振现象 (ω=%.3f ≈ ω_n=%.3f)', omega, omega_n);
                    else
                        title_str = sprintf('强迫振动响应 (ω=%.3f, ω_n=%.3f)', omega, omega_n);
                    end
                    plot_results(t, x, v, title_str, type);
                    
                case 4 % 动力型外力
                    F0 = str2double(get(F0_edit, 'String'));
                    omega = str2double(get(omega_edit, 'String'));
                    force_type_idx = force_type.Value;
                    
                    if isnan(F0)
                        error('外力幅值F0必须是数值');
                    end
                    
                    [t, x, v] = solve_dynamic_force(m, c, k, F0, omega, force_type_idx, x0, v0, t_span, selected_solver);
                    
                    force_types = {'正弦外力', '阶跃外力', '脉冲外力', '随机外力'};
                    plot_results(t, x, v, sprintf('%s响应', force_types{force_type_idx}), type);
                    
                case 6 % 静态外力
                    F0 = str2double(get(F0_edit, 'String'));
                    
                    if isnan(F0)
                        error('外力幅值F0必须是数值');
                    end
                    
                    [t, x, v] = solve_static_force(m, c, k, F0, x0, v0, t_span, selected_solver);
                    plot_results(t, x, v, '静态外力响应', type);
                    
                case 7 % 自定义外力函数
                    custom_force_str = strtrim(get(custom_force_edit, 'String'));
                    if isempty(custom_force_str)
                        error('请输入自定义外力函数');
                    end
                    
                    [t, x, v] = solve_custom_force(m, c, k, custom_force_str, x0, v0, t_span, selected_solver);
                    plot_results(t, x, v, '自定义外力响应', type);
            end
            
            set(status_text, 'String', sprintf('计算完成 (使用求解器: %s)', selected_solver),...
                'ForegroundColor',[0.2 0.5 0.2]);
            
        catch ME
            set(status_text, 'String', sprintf('错误: %s', ME.message), 'ForegroundColor',[0.8 0 0]);
            disp(getReport(ME, 'extended'));
        end
    end

    % 自由振动求解函数 (优化数值稳定性)
    function [t, x, v] = solve_free_vibration(m, c, k, x0, v0, t_span, solver)
        % 转换为状态空间形式
        ode_fun = @(t,y) [y(2); (-c*y(2) - k*y(1))/m];
        
        % 根据系统特性和用户选择确定求解器
        [ode_solver, opts] = select_solver(m, c, k, solver);
        
        % 求解ODE
        [t, y] = ode_solver(ode_fun, t_span, [x0; v0], opts);
        
        x = y(:,1);
        v = y(:,2);
    end

    % 强迫振动(正弦外力)求解函数
    function [t, x, v] = solve_forced_vibration(m, c, k, F0, omega, x0, v0, t_span, solver)
        % 外力函数 F(t) = F0*sin(ωt)
        F = @(t) F0*sin(omega*t);
        
        % 状态空间方程
        ode_fun = @(t,y) [y(2); (F(t) - c*y(2) - k*y(1))/m];
        
        % 根据系统特性和用户选择确定求解器
        [ode_solver, opts] = select_solver(m, c, k, solver);
        
        % 求解ODE
        [t, y] = ode_solver(ode_fun, t_span, [x0; v0], opts);
        
        x = y(:,1);
        v = y(:,2);
    end

    % 动力型外力求解函数
    function [t, x, v] = solve_dynamic_force(m, c, k, F0, omega, force_type, x0, v0, t_span, solver)
        % 根据外力类型定义不同的外力函数
        switch force_type
            case 1 % 正弦
                F = @(t) F0*sin(omega*t);
            case 2 % 阶跃
                F = @(t) F0*(t>=1); % 在t=1时施加阶跃力
            case 3 % 脉冲
                pulse_width = 0.1;
                F = @(t) F0*(abs(t-1)<pulse_width); % 在t=1附近施加短时脉冲
            case 4 % 随机
                rng(0); % 固定随机种子以便重复
                F = @(t) F0*interp1(linspace(t_span(1), t_span(end), 100), randn(1,100), t, 'pchip');
        end
        
        % 状态空间方程
        ode_fun = @(t,y) [y(2); (F(t) - c*y(2) - k*y(1))/m];
        
        % 根据系统特性和用户选择确定求解器
        [ode_solver, opts] = select_solver(m, c, k, solver, force_type);
        
        % 求解ODE
        [t, y] = ode_solver(ode_fun, t_span, [x0; v0], opts);
        
        x = y(:,1);
        v = y(:,2);
    end

    % 静态外力求解函数
    function [t, x, v] = solve_static_force(m, c, k, F0, x0, v0, t_span, solver)
        % 外力是常数 F(t) = F0
        F = @(t) F0;
        
        % 状态空间方程
        ode_fun = @(t,y) [y(2); (F(t) - c*y(2) - k*y(1))/m];
        
        % 根据系统特性和用户选择确定求解器
        [ode_solver, opts] = select_solver(m, c, k, solver);
        
        % 求解ODE
        [t, y] = ode_solver(ode_fun, t_span, [x0; v0], opts);
        
        x = y(:,1);
        v = y(:,2);
    end

    % 自定义外力函数求解
    function [t, x, v] = solve_custom_force(m, c, k, custom_force_str, x0, v0, t_span, solver)
        % 替换变量名为MATLAB可识别的形式
        custom_force_str = strrep(custom_force_str, 'x''', 'xd');
        custom_force_str = strrep(custom_force_str, 'x"', 'xdd');
        custom_force_str = strrep(custom_force_str, 'x', 'x');
        
        % 创建外力函数
        try
            % 测试表达式是否有效
            test_expr = str2func(['@(t,x,xd,xdd) ' custom_force_str]);
            test_val = test_expr(0,1,0,0); % 测试在t=0,x=1,xd=0,xdd=0时的值
            
            if ~isscalar(test_val) || ~isnumeric(test_val)
                error('自定义外力函数必须返回标量数值');
            end
        catch ME
            error('自定义外力函数格式错误: %s', ME.message);
        end
        
        % 状态空间方程
        ode_fun = @(t,y) [y(2); 
                         (test_expr(t,y(1),y(2),(-c*y(2) - k*y(1))/m) - c*y(2) - k*y(1))/m];
        
        % 根据系统特性和用户选择确定求解器
        [ode_solver, opts] = select_solver(m, c, k, solver, [], custom_force_str);
        
        % 求解ODE
        [t, y] = ode_solver(ode_fun, t_span, [x0; v0], opts);
        
        x = y(:,1);
        v = y(:,2);
    end

    % 智能选择求解器函数
    function [ode_solver, opts] = select_solver(m, c, k, solver_choice, force_type, custom_force_str)
        % 默认选项
        opts = odeset('RelTol',1e-6, 'AbsTol',1e-8);
        
        % 如果用户没有选择"自动"，则直接返回选择的求解器
        if ~strcmp(solver_choice, 'auto')
            switch solver_choice
                case 'ode45'
                    ode_solver = @ode45;
                case 'ode23'
                    ode_solver = @ode23;
                case 'ode113'
                    ode_solver = @ode113;
                case 'ode15s'
                    ode_solver = @ode15s;
                case 'ode23s'
                    ode_solver = @ode23s;
            end
            return;
        end
        
        % 自动选择逻辑
        omega_n = sqrt(k/m); % 固有频率
        zeta = c/(2*sqrt(m*k)); % 阻尼比
        
        % 判断系统特性
        is_stiff = false;
        has_discontinuities = false;
        is_nonlinear = false;
        
        % 检查是否有不连续外力
        if nargin >= 5 && ~isempty(force_type)
            has_discontinuities = (force_type == 2 || force_type == 3); % 阶跃或脉冲
        end
        
        % 检查自定义外力是否非线性
        if nargin >= 6 && ~isempty(custom_force_str)
            is_nonlinear = contains(custom_force_str, {'^','*','/','sin','cos','exp'});
        end
        
        % 判断是否为刚性系统
        if zeta > 1.5 || (zeta > 0.5 && omega_n > 50) || has_discontinuities
            is_stiff = true;
        end
        
        % 根据系统特性选择求解器
        if is_stiff
            % 刚性系统或有不连续外力
            ode_solver = @ode15s;
            opts = odeset(opts, 'MaxStep',0.01);
        elseif has_discontinuities
            % 不连续外力
            ode_solver = @ode15s;
            opts = odeset(opts, 'MaxStep',0.01);
        elseif is_nonlinear
            % 非线性系统
            ode_solver = @ode45;
            opts = odeset(opts, 'RelTol',1e-6, 'AbsTol',1e-8);
        elseif omega_n > 100
            % 高频系统
            ode_solver = @ode113;
            opts = odeset(opts, 'RelTol',1e-6, 'AbsTol',1e-8);
        else
            % 默认情况 (非刚性、低频、线性)
            ode_solver = @ode45;
            opts = odeset(opts, 'RelTol',1e-8, 'AbsTol',1e-10);
        end
    end

    % 结果绘图函数 (优化显示)
    function plot_results(t, x, v, title_str, type)
        figure('Name', title_str, 'NumberTitle','off', 'Position',[100 100 1200 700],...
            'Color','white');
        
        % 时域响应图
        subplot(2,2,[1 2]);
        plot(t, x, 'LineWidth', 1.5, 'Color',[0 0.45 0.74]);
        xlabel('时间 t (s)', 'FontSize',11);
        ylabel('位移 x (m)', 'FontSize',11);
        title(title_str, 'FontSize',12, 'FontWeight','bold');
        grid on;
        box on;
        
        % 根据振动类型添加额外信息
        switch type
            case {1, 2} % 自由振动
                % 计算理论解参数
                m = str2double(get(m_edit, 'String'));
                k = str2double(get(k_edit, 'String'));
                c = 0;
                if type == 2
                    c = str2double(get(c_edit, 'String'));
                end
                
                omega_n = sqrt(k/m); % 固有频率
                zeta = c/(2*sqrt(m*k)); % 阻尼比
                
                if zeta < 1 % 欠阻尼
                    omega_d = omega_n*sqrt(1-zeta^2); % 阻尼固有频率
                    legend_str = sprintf('欠阻尼振动 (ζ=%.2f, ω_d=%.2f rad/s)', zeta, omega_d);
                elseif zeta == 1 % 临界阻尼
                    legend_str = sprintf('临界阻尼 (ζ=%.2f)', zeta);
                else % 过阻尼
                    legend_str = sprintf('过阻尼 (ζ=%.2f)', zeta);
                end
                legend(legend_str, 'Location','best', 'FontSize',10);
                
            case {3, 5} % 强迫振动或共振
                m = str2double(get(m_edit, 'String'));
                k = str2double(get(k_edit, 'String'));
                c = str2double(get(c_edit, 'String'));
                omega = str2double(get(omega_edit, 'String'));
                
                omega_n = sqrt(k/m); % 固有频率
                zeta = c/(2*sqrt(m*k)); % 阻尼比
                r = omega/omega_n; % 频率比
                
                % 计算理论稳态振幅
                X = (str2double(get(F0_edit, 'String'))/k) / sqrt((1-r^2)^2 + (2*zeta*r)^2);
                
                % 绘制稳态振幅参考线
                hold on;
                plot([t(1) t(end)], [X X], 'r--', 'LineWidth', 1);
                plot([t(1) t(end)], [-X -X], 'r--', 'LineWidth', 1);
                hold off;
                
                legend_str = sprintf('响应 (ω/ω_n=%.2f, ζ=%.2f)\n稳态振幅: %.3f m', r, zeta, X);
                legend({'瞬态响应', '稳态振幅'}, 'Location','best', 'FontSize',10);
                
            case 4 % 动力型外力
                legend('位移响应', 'Location','best', 'FontSize',10);
                
            case 6 % 静态外力
                % 计算静态位移
                x_static = str2double(get(F0_edit, 'String')) / k;
                
                % 绘制静态位移参考线
                hold on;
                plot([t(1) t(end)], [x_static x_static], 'r--', 'LineWidth', 1);
                hold off;
                
                legend({'动态响应', sprintf('静态平衡位置: %.3f m', x_static)},...
                    'Location','best', 'FontSize',10);
                
            case 7 % 自定义外力
                legend('位移响应', 'Location','best', 'FontSize',10);
                % 显示自定义外力函数
                annotation('textbox',[0.15 0.85 0.7 0.1], 'String',...
                    sprintf('外力函数: F(t,x,x'',x") = %s', get(custom_force_edit, 'String')),...
                    'FitBoxToText','on', 'EdgeColor','none', 'FontSize',10,...
                    'HorizontalAlignment','center');
        end
        
        % 速度-时间图
        subplot(2,2,3);
        plot(t, v, 'LineWidth', 1.5, 'Color',[0.85 0.33 0.1]);
        xlabel('时间 t (s)', 'FontSize',11);
        ylabel('速度 v (m/s)', 'FontSize',11);
        title('速度-时间曲线', 'FontSize',11);
        grid on;
        box on;
        
        % 相图
        subplot(2,2,4);
        plot(x, v, 'LineWidth', 1.5, 'Color',[0.93 0.69 0.13]);
        xlabel('位移 x (m)', 'FontSize',11);
        ylabel('速度 v (m/s)', 'FontSize',11);
        title('相图 (位移-速度)', 'FontSize',11);
        grid on;
        box on;
        
        % 对于无阻尼自由振动，相图应该是闭合曲线
        if type == 1
            axis equal;
        end
        
        % 添加整体标题
        sgtitle(sprintf('振动系统响应分析 - %s', title_str), 'FontSize',13, 'FontWeight','bold');
    end

    % 初始化界面
    update_parameters();
end