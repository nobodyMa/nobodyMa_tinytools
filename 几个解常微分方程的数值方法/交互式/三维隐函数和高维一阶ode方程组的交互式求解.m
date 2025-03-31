equation_solver_guiji
%%
function equation_solver_guiji
    % 创建主界面
    fig = figure('Name','高维方程求解器', 'Position',[100 100 900 600], 'NumberTitle','off');
    
    % 问题类型选择
    uicontrol('Style','text', 'Position',[20 550 200 20], 'String','选择问题类型:', 'FontSize',10);
    problem_type = uicontrol('Style','popup', 'Position',[20 520 200 30],...
        'String',{'隐式方程','ODE方程组'}, 'FontSize',10, 'Callback',@toggle_inputs);
    
    % 方程输入区域
    uicontrol('Style','text', 'Position',[20 470 200 20], 'String','输入方程:', 'FontSize',10);
    eq_edit = uicontrol('Style','edit', 'Position',[20 350 250 120], 'Tag','eq',...
        'HorizontalAlignment','left', 'Max',5, 'FontSize',10, 'Tooltip',...
        ['例如：\n2D: x^2 + y^2 -1\n'...
         '3D: x^2/4 + y^2/9 + z^2 -1\n'...
         'ODE: [y(2); -y(1)]']);
    
    % 初值条件输入（默认隐藏）
    iv_label = uicontrol('Style','text', 'Position',[20 300 200 20], 'String','初值条件:',...
        'Visible','off', 'Tag','iv_label', 'FontSize',10);
    iv_edit = uicontrol('Style','edit', 'Position',[20 250 250 50], 'Tag','iv',...
        'Visible','off', 'FontSize',10, 'Tooltip','列向量格式，如[1; 0; 0.5]');
    
    % 参数输入
    uicontrol('Style','text', 'Position',[20 200 200 20], 'String','参数范围:', 'FontSize',10);
    range_edit = uicontrol('Style','edit', 'Position',[20 170 250 30], 'FontSize',10,...
        'Tooltip',...
        ['3D格式: x范围; y范围; z范围\n'...
         '示例: -2:0.1:2; -3:0.1:3; -1:0.1:1']);
        
    % 绘图按钮
    plot_btn = uicontrol('Style','pushbutton', 'Position',[20 30 100 40],...
        'String','绘制', 'FontSize',11, 'Callback',@plot_solution,...
        'BackgroundColor',[0.8 0.9 0.8]);
    
    % 状态提示文本
    status_text = uicontrol('Style','text', 'Position',[20 80 250 30],...
        'String','', 'FontSize',10, 'ForegroundColor',[0.8 0 0]);
    
    % 坐标轴
    ax = axes('Parent',fig, 'Position',[0.4 0.1 0.55 0.8], 'Tag','main_ax',...
        'Box','on', 'FontSize',10);
    grid(ax, 'on');
    
    % 回调函数 - 切换输入类型
    function toggle_inputs(~,~)
        type = problem_type.Value;
        if type == 2 % ODE系统
            iv_label.Visible = 'on';
            iv_edit.Visible = 'on';
            set(range_edit, 'Tooltip','时间范围，如 0:0.01:50');
        else % 隐式方程
            iv_label.Visible = 'off';
            iv_edit.Visible = 'off';
            set(range_edit, 'Tooltip',...
                ['3D格式: x范围; y范围; z范围\n'...
                 '示例: -2:0.1:2; -3:0.1:3; -1:0.1:1']);
        end
    end

    % 主绘图函数
    function plot_solution(~,~)
        cla(ax); 
        set(status_text, 'String', '计算中，请稍候...');
        drawnow;
        
        try
            % 获取输入
            eq_str = strtrim(get(eq_edit, 'String'));
            range_str = strtrim(get(range_edit, 'String'));
            if isempty(eq_str) || isempty(range_str)
                error('输入不能为空');
            end
            
            % 调试信息
            disp(['=== 调试信息 ===']);
            disp(['方程: ', eq_str]);
            disp(['初值: ', get(iv_edit, 'String')]);
            disp(['范围: ', range_str]);
            
            if problem_type.Value == 1 % 隐式方程
                % 符号变量检测与验证
                vars = symvar(eq_str);
                valid_vars = {'x','y','z'};
                for i = 1:length(vars)
                    var_name = char(vars(i));
                    if ~any(strcmpi(var_name, valid_vars))
                        error('仅允许使用x,y,z变量');
                    end
                end
                dim = numel(vars);
                
                % 解析范围参数
                range_cells = strsplit(range_str, ';');
                range_cells = strtrim(range_cells(~cellfun('isempty', range_cells)));
                
                if dim == 3
                    % 3D隐式方程处理
                    if numel(range_cells) ~= 3
                        error('3D方程需要3个范围参数');
                    end
                    ranges = zeros(3, 3);
                    for i = 1:3
                        parts = strsplit(range_cells{i}, ':');
                        if numel(parts) ~= 3
                            error('范围格式应为start:step:end');
                        end
                        ranges(i,:) = str2double(parts);
                    end
                    x = ranges(1,1):ranges(1,2):ranges(1,3);
                    y = ranges(2,1):ranges(2,2):ranges(2,3);
                    z = ranges(3,1):ranges(3,2):ranges(3,3);
                    
                    % 计算量预估 
                    if total_calcs > 1e14 
                        % 使用隐式方程测试结果进行预估 (1.62e-07 * N^3 + 0.05)
                        est_time = 1.62e-07 * total_calcs + 0.05;
                        answer = questdlg(sprintf('预计计算量%.2e次，耗时约%.1f秒。继续吗?', total_calcs, est_time), ...
                            '计算量警告', '继续', '取消', '取消');
                        if strcmp(answer, '取消')
                            set(status_text, 'String', '计算已取消');
                            return;
                        end
                    end
                    
                    f = str2func(sprintf('@(%s) %s', strjoin(vars,','), eq_str));
                    fimplicit3(ax, f, [x(1) x(end) y(1) y(end) z(1) z(end)],...
                        'MeshDensity',40, 'FaceAlpha',0.7);
                    view(ax, 3);
                    title(ax, '3D隐式方程解曲面', 'FontSize',12);
                    
                elseif dim == 2
                    % 2D隐式方程处理
                    if numel(range_cells) > 1
                        range_cells = range_cells(1);
                    end
                    parts = strsplit(range_cells{1}, ':');
                    x_range = str2double(parts(1)):str2double(parts(2)):str2double(parts(3));
                    
                    f = str2func(sprintf('@(%s) %s', strjoin(vars,','), eq_str));
                    fimplicit(ax, f, [x_range(1) x_range(end)], 'LineWidth', 2);
                    title(ax, '2D隐式方程解曲线', 'FontSize',12);
                end
                
                % 统一图形设置
                xlabel(ax, 'x'); ylabel(ax, 'y'); 
                if dim == 3, zlabel(ax, 'z'); end
                grid(ax, 'on');
                axis(ax, 'tight');
                
            else % ODE方程组
    % ODE输入验证 - 增强灵活性
    y0_str = strtrim(get(iv_edit, 'String'));
    try
        y0 = eval(y0_str); % 支持更灵活的输入格式
        if isempty(y0) || ~isnumeric(y0)
            error('初值必须是数值向量');
        end
        y0 = y0(:); % 转换为列向量
    catch ME
        error('初值格式错误，请输入数值向量，如 [1;0] 或 [1,0]''');
    end
    
    % 时间范围解析增强 - 支持更多格式
    try
        if contains(range_str, ':')
            t_span = eval(range_str); % 支持0:0.01:10格式
        else
            t_span = str2num(range_str); % 支持[0 1 2 3]格式
        end
        if isempty(t_span)
            error('时间范围格式错误');
        end
    catch
        error('时间范围格式错误，示例: 0:0.01:10 或 [0,1,2,3]');
    end
    
    % 改进的ODE计算量预估模型
    if numel(t_span) > 1e4
        % 基于系统复杂度的预估
        dim = numel(y0);
        eq_str_lower = lower(eq_str);
        
        % 复杂度评估
        nonlinear_terms = sum(contains(eq_str_lower, {'^','*','/','sin','cos','exp'}));
        complexity = 1 + nonlinear_terms*0.2;
        
        % 预估时间 (基准时间 * 维度因子 * 复杂度因子)
        base_time = 0.01 * numel(t_span)/1e4;
        est_time = base_time * (1 + 0.1*dim^1.5) * complexity;
        
        if est_time > 5
            answer = questdlg(sprintf([...
                '预估计算时间: %.1f秒\n'...
                '系统维度: %d\n'...
                '非线性项: %d\n'...
                '确认要继续计算吗？'], ...
                est_time, dim, nonlinear_terms), ...
                '长时间计算警告', '继续', '取消', '取消');
            
            if strcmp(answer, '取消')
                set(status_text, 'String', '计算已取消');
                return;
            end
        end
        set(status_text, 'String', sprintf('预估时间: %.2f秒 | 正在计算...', est_time));
        drawnow;
    end
    
    % 动态创建ODE函数 - 增加更多错误检查
    try
        % 替换常见数学函数写法
        eq_str = regexprep(eq_str, '(\d+)e([+-]\d+)', '$1*10^($2)');
        ode_fun = str2func(['@(t,y) ', eq_str]);
        
        % 测试函数有效性
        try
            test_output = ode_fun(0, y0);
        catch ME
            error('方程测试失败: %s\n检查变量是否使用y(1),y(2)...格式', ME.message);
        end
        
        if numel(test_output) ~= numel(y0)
            error('方程输出维度(%d)与初值维度(%d)不匹配',...
                numel(test_output), numel(y0));
        end
    catch ME
        error('ODE函数创建失败: %s', ME.message);
    end
    
    % 求解ODE - 根据系统特性自动选择算法
    opts = odeset('RelTol',1e-6, 'AbsTol',1e-9);
    try
        [t, y] = ode45(ode_fun, t_span, y0, opts);
    catch ME
        % 若ode45失败，尝试使用刚性求解器
        try
            [t, y] = ode15s(ode_fun, t_span, y0, opts);
            set(status_text, 'String', [get(status_text, 'String'), ' (使用刚性求解器)']);
        catch
            rethrow(ME);
        end
    end
    
    % 优化结果可视化
    ode_fig = figure('Name','ODE Solution', 'NumberTitle','off', 'Position',[100 100 900 500]);
    if size(y,2) <= 5
        % 时域图
        subplot(1,2,1);
        plot(t, y, 'LineWidth', 1.5);
        xlabel('时间 t'); ylabel('状态');
        title('时域解');
        grid on;
        
        % 添加图例
legend(arrayfun(@(n)sprintf('y_%d',n), 1:size(y,2), 'UniformOutput',false), 'Location','best');

        % 相空间图
        subplot(1,2,2);
        if size(y,2) >= 3
            plot3(y(:,1), y(:,2), y(:,3), 'LineWidth', 1.5);
            xlabel('y_1'); ylabel('y_2'); zlabel('y_3');
            view(3);
            title('三维相空间');
        elseif size(y,2) == 2
            plot(y(:,1), y(:,2), 'LineWidth', 1.5);
            xlabel('y_1'); ylabel('y_2');
            title('相平面');
        else
            plot(t, y, 'LineWidth', 1.5);
            xlabel('时间 t'); ylabel('y');
            title('标量解');
        end
        grid on;
        axis equal;
        
        % 添加图例
        if size(y,2) <= 5
            legend(arrayfun(@(n)sprintf('y_%d',n), 1:size(y,2), 'UniformOutput',false),...
                'Location','best');
        end
    else
        % 高维系统简化显示
        plot(t, y, 'LineWidth', 1);
        xlabel('时间 t'); ylabel('状态');
        title(sprintf('%d维系统时域解', size(y,2)));
        grid on;
    end
    
    % 更新状态
    set(status_text, 'String', sprintf('计算完成 (结果在Figure %s)', ode_fig.Name));
end
            
        catch ME
            set(status_text, 'String', sprintf('错误: %s', ME.message));
            disp(getReport(ME, 'extended')); % 在命令行显示完整错误堆栈
        end
    end

    % 初始化界面
    toggle_inputs();
end