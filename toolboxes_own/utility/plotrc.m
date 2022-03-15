function fig_ = plotrc(rc,dates)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

% vechRC = timetable(dates,vech3d(rc));
% vechRC = retime(vechRC, 'daily', 'fillwithmissing');
% dates = vechRC.dates;
% rc = ivech3d(vechRC{:,1});

volas = sqrt(diag3d(rc));

corrs = tril3d(cov2corr3d(rc),-1);

fig_ = figure;
% plot(dates, diag3d(rc), 'Color', [0 0.4470 0.7410, .2])

t = tiledlayout(3,1,'TileSpacing','Compact','Padding','Compact');

nexttile
plot(dates, logdet3d(rc))%,'Marker','.','MarkerSize',5, 'MarkerEdgeColor', [1 0 0])
axis_ = gca;
set(axis_, ...%'yscale','log', ...
        'Color', [0.85, 0.85, 0.85], 'XGrid', 'on', 'YGrid', 'on', ...
        'Box', 'off');
% ylabel('det(RC)')    
axis_.XAxis.Visible = 'off';
axis_.YAxis.Visible = 'off';
axis_.GridColor = [1 1 1];
axis_.MinorGridColor = [1 1 1];
axis_.GridAlpha = 0.5;
text([-10, -10, -10, -10, -10, -10, -10], ...
     [300, 400, 500, 600, 700, 800, 900], ...
     {'300', '400', '500', '600', '700', '800', '900'}, ...
     'Horiz', 'right', 'Vert', 'middle', 'FontSize', 8)
xt = axis_.XTick;
xtv = datestr(xt,'yyyy');
text(xt, axis_.YLim(1)*ones(size(xt)), xtv, 'Horiz','center', 'Vert','bottom', 'Interpreter','latex')
text(3500,800,'logdet($\mathbf{RC}$)','Interpreter','latex','HorizontalAlignment','center')

nexttile
my_stairs(dates, volas, 'Color', [0 0.4470 0.7410], 'LineWidth', 1.5);
axis_ = gca;
set(axis_, ...%'yscale','log', ...
        'Color', [0.85, 0.85, 0.85], 'XGrid', 'on', 'YGrid', 'on', ...
        'Box', 'off');
% ylabel('RC_{ii}');    
axis_.XAxis.Visible = 'off';
axis_.YAxis.Visible = 'off';
yticks([1e-1,1,1e1,1e2])
axis_.GridColor = [1 1 1];
axis_.MinorGridColor = [0.85, 0.85, 0.85];
axis_.GridAlpha = 0.5;
text([-10, -10, -10, -10], [1e-1,1,1e1,1e2], {'10^{-1}', '10^{0}', '10^{1}', '10^{2}'}, 'Horiz', 'right', 'Vert', 'middle', 'FontSize', 8)
xt = axis_.XTick;
xtv = datestr(xt,'yyyy');
text(xt, 10^(-2)*ones(size(xt)), xtv, 'Horiz','center', 'Vert','top','Interpreter','latex')
text(3670,375,'$\mathbf{RC}_{ii}$','Interpreter','latex','HorizontalAlignment','center')

nexttile
my_stairs(dates, corrs, 'Color', [0 0.4470 0.7410])
axis_ = gca;
set(axis_, 'Color', [0.85, 0.85, 0.85], 'XGrid', 'on', 'YGrid', 'on', ...
         'Box', 'off');
% ylabel('Corr_{ij}')     
axis_.XAxis.Visible = 'off';
axis_.YAxis.Visible = 'off';
axis_.GridColor = [1 1 1];
axis_.GridAlpha = 0.5;
text([-10, -10, -10, -10, -10], [-1,-0.5,0,0.5,1], {'-1', '-0.5', '0', '0.5', '1'}, 'Horiz', 'right' , 'Vert', 'middle', 'FontSize', 8)
xt = axis_.XTick;
xtv = datestr(xt,'yyyy');
text(xt, axis_.YLim(1)*ones(size(xt)), xtv, 'Horiz','center', 'Vert','top','Interpreter','latex')
text(3670,-.6,'$\mathbf{RCorr}_{ij}$','Interpreter','latex','HorizontalAlignment','center')

end
    
