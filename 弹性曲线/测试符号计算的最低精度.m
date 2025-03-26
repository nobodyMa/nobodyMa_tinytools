all clear ;all clc;
test_elliptic_precision()
%%
function test_elliptic_precision()
    % 基本参数设置
    L = 1.0;               % 长度参数
    m = 0.5;               % 椭圆积分参数
    x_test = 0.7*L;        % 测试点位置
    precisions = [16, 32, 64, 128]; % 测试的精度位数
    
    % 计算theta(x)
    syms t;
    theta_expr = (m*L/2)*(2*(x_test/L)-(x_test/L)^2);
    theta_val = double(asin(theta_expr));
    
    % 预存储结果
    results = zeros(length(precisions), 3); % [精度, F值, E值]
    rel_diffs = zeros(length(precisions)-1, 3); % 相邻精度相对差异
    
    % 主测试循环
    fprintf('%-10s %-25s %-25s %-15s\n', 'Precision', 'F(θ,m)', 'E(θ,m)', 'y(x)');
    for i = 1:length(precisions)
        digits(precisions(i));
        % 修正：使用合理的容差值，不能为0
        opts = {'RelTol', 1e-20, 'AbsTol', 10^(-precisions(i))};
        
        % 高精度符号计算
        F_sym = vpaintegral(1/sqrt(1 - m*sin(t)^2), t, 0, theta_val, opts{:});
        E_sym = vpaintegral(sqrt(1 - m*sin(t)^2), t, 0, theta_val, opts{:});
        y_val = double((F_sym - E_sym)/sqrt(m));
        
        % 存储结果
        results(i,:) = [precisions(i), double(F_sym), double(E_sym)];
        
        % 打印当前结果
        fprintf('%-10d %-25.15e %-25.15e %-15.10f\n', ...
                precisions(i), double(F_sym), double(E_sym), y_val);
        
        % 计算与上一个精度的相对差异（从第二个精度开始）
        if i > 1
            rel_diffs(i-1,1) = abs(results(i,2)-results(i-1,2))/abs(results(i-1,2));
            rel_diffs(i-1,2) = abs(results(i,3)-results(i-1,3))/abs(results(i-1,3));
            rel_diffs(i-1,3) = precisions(i);
        end
    end
    
    % ========== 可视化分析 ==========
    figure('Position', [100, 100, 1200, 500]);
    
    % 结果值变化
    subplot(1,2,1);
    semilogy(precisions, abs(results(:,2)), 'bo-', 'DisplayName', 'F(θ,m)');
    hold on;
    semilogy(precisions, abs(results(:,3)), 'r^-', 'DisplayName', 'E(θ,m)');
    xlabel('计算精度（有效位数）');
    ylabel('积分值大小');
    title('不同精度下的积分值');
    legend('Location', 'best');
    grid on;
    
    % 相邻精度相对差异
    subplot(1,2,2);
    semilogy(rel_diffs(:,3), rel_diffs(:,1), 'bo-', 'DisplayName', 'F(θ,m)相对差异');
    hold on;
    semilogy(rel_diffs(:,3), rel_diffs(:,2), 'r^-', 'DisplayName', 'E(θ,m)相对差异');
    xlabel('较高精度（有效位数）');
    ylabel('与前一精度的相对差异');
    title('精度提升带来的改进');
    legend('Location', 'best');
    grid on;
    
    % 控制台输出关键结论
    fprintf('\n===== 关键结论 =====\n');
    fprintf('1. F(θ,m) 的收敛精度: %.2e (在%d位后稳定)\n', ...
            min(rel_diffs(:,1)), find(rel_diffs(:,1)<1e-15,1)+1);
    fprintf('2. E(θ,m) 的收敛精度: %.2e (在%d位后稳定)\n', ...
            min(rel_diffs(:,2)), find(rel_diffs(:,2)<1e-15,1)+1);
    fprintf('3. 推荐工作精度: %d位（实际误差<1e-%d）\n', ...
            precisions(find(rel_diffs(:,1)<1e-15,1)), round(-log10(min(rel_diffs(:,1)))));
end