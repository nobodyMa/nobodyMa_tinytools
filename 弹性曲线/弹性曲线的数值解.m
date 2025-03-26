all close ;all clear;all clc
compare_elliptic_integration_methods()
%%
function compare_elliptic_integration_methods()
    % 参数设置
    L = 1.0;       % 长度尺度参数
    m = 0.5;       % 椭圆积分参数 (0 < m < 1)
    x_values = linspace(0, L, 50);  % x的取值点
    
    % 预分配结果数组
    y_builtin = zeros(size(x_values));
    y_quadgk = zeros(size(x_values));
    y_integral = zeros(size(x_values));
    y_symbolic = zeros(size(x_values));
    
    % 符号计算初始化
    syms t;
    m_sym = sym(m);
    L_sym = sym(L);
    opts = {'RelTol',1e-20,'AbsTol',1e-20};
    
    % 主计算循环
    for i = 1:length(x_values)
        x = x_values(i);
        x_sym = vpa(x,32);
        
        % 计算theta(x)
        theta_expr = (m_sym*L_sym/2)*(2*(x_sym/L_sym)-(x_sym/L_sym)^2);
        theta_val = double(asin(theta_expr));
        theta_sym = asin(theta_expr);
        
        % 方法1: 内置椭圆积分函数
        F_builtin = ellipticF(theta_val, m);
        E_builtin = ellipticE(theta_val, m);
        y_builtin(i) = (F_builtin - E_builtin)/sqrt(m);
        
        % 方法2: 自适应高斯-克朗罗德积分(quadgk)
        integrand_F = @(t) 1./sqrt(1 - m*sin(t).^2);
        integrand_E = @(t) sqrt(1 - m*sin(t).^2);
        F_quadgk = quadgk(integrand_F, 0, theta_val, 'AbsTol',1e-12,'RelTol',1e-12);
        E_quadgk = quadgk(integrand_E, 0, theta_val, 'AbsTol',1e-12,'RelTol',1e-12);
        y_quadgk(i) = (F_quadgk - E_quadgk)/sqrt(m);
        
        % 方法3: 全局自适应积分(integral)
        F_integral = integral(integrand_F, 0, theta_val, 'AbsTol',1e-12,'RelTol',1e-12);
        E_integral = integral(integrand_E, 0, theta_val, 'AbsTol',1e-12,'RelTol',1e-12);
        y_integral(i) = (F_integral - E_integral)/sqrt(m);
        
        % 方法4: 符号计算（基准）
        F_sym = vpaintegral(1/sqrt(1 - m_sym*sin(t)^2), t, 0, theta_sym, opts{:});
        E_sym = vpaintegral(sqrt(1 - m_sym*sin(t)^2), t, 0, theta_sym, opts{:});
        y_symbolic(i) = double((F_sym - E_sym)/sqrt(m_sym));
    end
    
    % ========== 生成LaTeX表达式 ==========
    fprintf(['\n%% ========== LaTeX 表达式输出 ==========\n' ...
        '（使用两个下划线表示LaTeX中的单个下划线）\n']);
    
    % 1. 数据点向量
    fprintf('%% 数据点向量\n');
    fprintf('x__values = [%.4f', x_values(1));
    for i = 2:length(x_values)
        fprintf(', %.4f', x_values(i));
    end
    fprintf('];\n\n');
    
    % 2. 各方法结果
    fprintf('%% 符号计算结果\n');
    fprintf('y__symbolic = [%.8e', y_symbolic(1));
    for i = 2:length(y_symbolic)
        fprintf(', %.8e', y_symbolic(i));
    end
    fprintf('];\n\n');
    
    fprintf('%% 内置椭圆积分函数结果\n');
    fprintf('y__builtin = [%.8e', y_builtin(1));
    for i = 2:length(y_builtin)
        fprintf(', %.8e', y_builtin(i));
    end
    fprintf('];\n\n');
    
    fprintf('%% 自适应高斯-克朗罗德积分结果\n');
    fprintf('y__quadgk = [%.8e', y_quadgk(1));
    for i = 2:length(y_quadgk)
        fprintf(', %.8e', y_quadgk(i));
    end
    fprintf('];\n\n');
    
    fprintf('%% 全局自适应积分结果\n');
    fprintf('y__integral = [%.8e', y_integral(1));
    for i = 2:length(y_integral)
        fprintf(', %.8e', y_integral(i));
    end
    fprintf('];\n\n');
    
    % 3. LaTeX表格
    fprintf('%% LaTeX表格代码\n');
    fprintf('\\begin{tabular}{|c|c|c|c|c|}\n');
    fprintf('\\hline\n');
    fprintf('$x$ & 符号计算 & 内置椭圆积分函数 & 自适应高斯-克朗罗德积分 & 全局自适应积分 \\\\\n');
    fprintf('\\hline\n');
    for i = 1:length(x_values)
        fprintf('%.4f & %.4e & %.4e & %.4e & %.4e \\\\\n', ...
            x_values(i), y_symbolic(i), y_builtin(i), y_quadgk(i), y_integral(i));
    end
    fprintf('\\hline\n');
    fprintf('\\end{tabular}\n');
    
    % 多项式拟合表达式
    fprintf('%% 结果曲线多项式拟合表达式\n');
    
    % 对符号基准结果进行多项式拟合
    p = polyfit(x_values, y_symbolic, 6); % 6次多项式拟合
    syms x;
    poly_expr = poly2sym(p, x);
    fprintf('\n6次多项式拟合表达式:\n ');
    latex_str = latex(vpa(expand(poly_expr), 5));
    fprintf('$$ y_{\\mathrm{fit}}(x) = %s $$\n', latex_str);
    
    % 显示拟合精度
    y_fit = polyval(p, x_values);
    fit_error = max(abs(y_fit - y_symbolic));
    fprintf('\n最大拟合误差: %.2e\n', fit_error);

    % ========== 可视化 ==========
    %结果曲线
    figure('Position',[100 100 1200 600])
    subplot(1,2,1)
    plot(x_values, y_symbolic, 'k-', 'LineWidth',2, 'DisplayName','Symbolic');
    hold on;
    plot(x_values, y_builtin, 'bo', 'DisplayName','内置椭圆积分函数');
    plot(x_values, y_quadgk, 'r^', 'DisplayName','自适应高斯-克朗罗德积分');
    plot(x_values, y_integral, 'gd', 'DisplayName','全局自适应积分');
    xlabel('$x$', 'Interpreter','latex');
    ylabel('$y(x)$', 'Interpreter','latex');
    title('Comparison of Four Methods');
    legend('Location','best');
    grid on;
    
    % 误差分析
    subplot(1,2,2)
    rel_err = @(y) abs(y-y_symbolic)./abs(y_symbolic);
    semilogy(x_values, rel_err(y_builtin), 'bo-', 'DisplayName','内置椭圆积分函数');
    hold on;
    semilogy(x_values, rel_err(y_quadgk), 'r^-', 'DisplayName','自适应高斯-克朗罗德积分');
    semilogy(x_values, rel_err(y_integral), 'gd-', 'DisplayName','全局自适应积分');
    xlabel('$x$', 'Interpreter','latex');
    ylabel('Relative Error');
    title('Relative Error Comparison');
    legend('Location','best');
    grid on;
end