function free_adjust_planes
    % 创建主窗口
    fig = figure('Name','Free Adjustment of Planes', 'NumberTitle','off');
    ax = axes('Parent', fig, 'Position', [0.1 0.4 0.8 0.5]);
    hold(ax, 'on');
    grid(ax, 'on');
    view(ax, 3);
    xlabel(ax, 'X');
    ylabel(ax, 'Y');
    zlabel(ax, 'Z');
    axis(ax, [-10 10 -10 10 -10 10]);
    
    % 初始化参数
    a = 0;
    b = 0;
    c = 0;
    
    % 创建滑动条
    slider_a = uicontrol('Parent', fig, 'Style', 'slider', ...
        'Min', -10, 'Max', 10, 'Value', a, ...
        'Position', [100 80 300 20], ...
        'Callback', @update_a);
    slider_b = uicontrol('Parent', fig, 'Style', 'slider', ...
        'Min', -10, 'Max', 10, 'Value', b, ...
        'Position', [100 50 300 20], ...
        'Callback', @update_b);
    slider_c = uicontrol('Parent', fig, 'Style', 'slider', ...
        'Min', -10, 'Max', 10, 'Value', c, ...
        'Position', [100 20 300 20], ...
        'Callback', @update_c);
    
    % 添加滑动条标签
    uicontrol('Parent', fig, 'Style', 'text', ...
        'Position', [20 80 60 20], ...
        'String', 'a');
    uicontrol('Parent', fig, 'Style', 'text', ...
        'Position', [20 50 60 20], ...
        'String', 'b');
    uicontrol('Parent', fig, 'Style', 'text', ...
        'Position', [20 20 60 20], ...
        'String', 'c');
    
    % 创建图形对象
    [y_grid, z_grid] = meshgrid(linspace(-10,10,20));
    x_plane = a * ones(size(y_grid));
    surf_x = surf(ax, x_plane, y_grid, z_grid, 'FaceColor', 'r', 'FaceAlpha', 0.3);
    
    [x_grid, z_grid] = meshgrid(linspace(-10,10,20));
    y_plane = b * ones(size(x_grid));
    surf_y = surf(ax, x_grid, y_plane, z_grid, 'FaceColor', 'g', 'FaceAlpha', 0.3);
    
    [x_grid, y_grid] = meshgrid(linspace(-10,10,20));
    z_plane = c * ones(size(x_grid));
    surf_z = surf(ax, x_grid, y_grid, z_plane, 'FaceColor', 'b', 'FaceAlpha', 0.3);
    
    scatter_A = scatter3(ax, a, 0, 0, 'ro', 'filled');
    scatter_B = scatter3(ax, 0, b, 0, 'go', 'filled');
    scatter_C = scatter3(ax, 0, 0, c, 'bo', 'filled');
    scatter_P = scatter3(ax, a, b, c, 'k', 'filled', 'SizeData', 100);
    
    line_PA = plot3(ax, [a, a], [0, b], [0, c], 'k--');
    line_PB = plot3(ax, [0, a], [b, b], [0, c], 'k--');
    line_PC = plot3(ax, [0, a], [0, b], [c, c], 'k--');
    
    text_dist = text(ax, a, b, c+1, '', 'FontSize', 10, 'Color', 'k');
    
    % 嵌套回调函数
    function update_a(~,~)
        a = get(slider_a, 'Value');
        update_plot();
    end
    function update_b(~,~)
        b = get(slider_b, 'Value');
        update_plot();
    end
    function update_c(~,~)
        c = get(slider_c, 'Value');
        update_plot();
    end
    
    function update_plot()
        % 更新点A、B、C和P的位置
        set(scatter_A, 'XData', a, 'YData', 0, 'ZData', 0);
        set(scatter_B, 'XData', 0, 'YData', b, 'ZData', 0);
        set(scatter_C, 'XData', 0, 'YData', 0, 'ZData', c);
        set(scatter_P, 'XData', a, 'YData', b, 'ZData', c);
        
        % 更新平面
        [y, z] = meshgrid(linspace(-10,10,20));
        x_new = a * ones(size(y));
        set(surf_x, 'XData', x_new, 'YData', y, 'ZData', z);
        
        [x, z] = meshgrid(linspace(-10,10,20));
        y_new = b * ones(size(x));
        set(surf_y, 'XData', x, 'YData', y_new, 'ZData', z);
        
        [x, y] = meshgrid(linspace(-10,10,20));
        z_new = c * ones(size(x));
        set(surf_z, 'XData', x, 'YData', y, 'ZData', z_new);
        
        % 更新连接线
        set(line_PA, 'XData', [a, a], 'YData', [0, b], 'ZData', [0, c]);
        set(line_PB, 'XData', [0, a], 'YData', [b, b], 'ZData', [0, c]);
        set(line_PC, 'XData', [0, a], 'YData', [0, b], 'ZData', [c, c]);
        
        % 计算距离
        PA_val = sqrt(b^2 + c^2);
        PB_val = sqrt(a^2 + c^2);
        PC_val = sqrt(a^2 + b^2);
        
        % 更新文本
        set(text_dist, 'Position', [a, b, c+1], ...
            'String', sprintf('PA=%.2f\nPB=%.2f\nPC=%.2f', PA_val, PB_val, PC_val));
        
        drawnow;
    end
end