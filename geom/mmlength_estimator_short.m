% author: ziyan (zoe) zhu
% email: zzhu1@g.harvard.edu
% Fig. 1b short paper  
clear all
f_size=22;
set(groot, 'DefaultTextInterpreter', 'Latex')
set(groot, 'DefaultLegendInterpreter', 'Latex')
set(groot, 'DefaultAxesTickLabelInterpreter', 'Latex')
set(0,'DefaultAxesFontSize',f_size)


q1_list = [1:0.01:5];
q2_list = q1_list;
% 
% q1_list = 2.6;
% q2_list = 3.3;

a0=2.4728;
A0 = a0*[1 1/2;
         0 sqrt(3)/2];

theta1=2.6;
theta2=2.8;

for q1_idx = 1:length(q1_list)
    for q2_idx = 1:length(q2_list)
        q12 = deg2rad(-q1_list(q1_idx));
        q23 = deg2rad(q2_list(q2_idx));

        Am_dom = calc_moire_dom_tri(q12, q23, A0);
        ml_dom(q1_idx, q2_idx) = norm(Am_dom(:,1));
    end 
end 

[q1grid, q2grid] = meshgrid(q1_list, q2_list);

clevels = [0:10:20,40:20:160];

%%
figure;
hold all

[C,h]=contourf(q1grid, q2grid, ml_dom/10, clevels);
scatter3(theta1,theta2,max(ml_dom(:)),500,'k','p','filled');
c = colorbar;
c.TickLabelInterpreter = 'latex';
c.Label.Interpreter = 'latex';
c.Label.String = 'Dominant moir\''e of moir\''e length (nm)';
c.Ticks = [10:20:160];
xlabel('$\theta_{12}$');
ylabel('$\theta_{23}$');
box on;
axis equal
yticks(xticks)
set(gca,'ColorScale','log')

cmap=colormap(brewermap([],'Spectral'));
cmap=flip(cmap,1);
colormap(cmap);
% colormap(jet);
c1 = 'k';

text(4.05, 4.6, '\textbf{(1,1)}', 'Color', c1,'FontSize',f_size*0.8,'Interpreter','latex');
text(2.1, 4.5, '\textbf{(2,1)}', 'Color', c1, 'FontSize',f_size*0.8,'Interpreter','latex');
text(1.45, 4.8, '\textbf{(3,1)}', 'Color', c1,'FontSize',f_size*0.8,'Interpreter','latex');
text(1.3, 3.7, '\textbf{(3,2)}', 'Color', c1,'FontSize',f_size*0.8,'Interpreter','latex');
ylim([min(abs(q1_list)),max(abs(q1_list))])

q_eq = linspace(min(q1grid(:)), max(q1grid(:)));
plot3(q_eq, q_eq, linspace(1,1)'.*max(ml_dom(:)), 'k--', 'LineWidth', 3);

set(gca, 'Layer', 'top');
