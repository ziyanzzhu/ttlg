clear all
%close all
f_size = 22;
set(groot, 'DefaultTextInterpreter', 'Latex')
set(groot, 'DefaultLegendInterpreter', 'Latex')
set(groot, 'DefaultAxesTickLabelInterpreter', 'Latex')
set(0,'DefaultAxesFontSize',f_size)


theta_list = [-2 0 2]; % twisting angles (global)
k_cutoff = 3;      % cutoff in momentum space
grid_search = 20;    % [G1,G2] are in [-grid_search,grid_search]^2 

proj = [1, 2, 3];  % sheet to project eigenvector weights onto
q_cut_type = 5; % what kind of line-cut we do in momentum space
savedata = 1; 
savepng = 0;
weight_cut = 0; % only plot parts of the bands, not very useful
wcut = 0.2;
E_field = 0.0; % electric field 
num_eigs = 120;

int_off = 0; 

nq = 10; % number of k points to sample

colors = [.7, .2, 0; 
          0, .7, .2;
          .2, 0, .7];
      
if length(proj) == 3
    figname = ['all_proj_q12_' num2str(-theta_list(1)) '_q23_' num2str(abs(theta_list(3)))...
                '_kcut_' num2str(k_cutoff) '_qtype_' num2str(q_cut_type)];

else 
    figname = ['q12_' num2str(-theta_list(1)) '_q23_' num2str(abs(theta_list(3)))...
                '_kcut_' num2str(k_cutoff) '_qtype_' num2str(q_cut_type) ''];
end 
        
if int_off
    figname = ['int12off' figname];
end 
      
tic 

% create layer data structures
for t = 1:3
   layers(t) = Layer(t,deg2rad(theta_list(t)));
end

G1 = layers(1).G; 
G2 = layers(2).G;
G3 = layers(3).G;
b12 = G1 - G2; 
b23 = G2 - G3; 
b_tri = b12-b23;

G11 = G2(:, 1);
G12 = G2(:, 2);
G13 = G11 + G12; 
K0 = 1/3 * (G11 + G13);

% get the K-point of each monolayer
for t = 1:3
   th = layers(t).theta;
   K(t,:) = [cos(th) -sin(th); sin(th) cos(th)]*K0;
end

% define the k-point sampling
samp = linspace(0,1,nq)';
samp = samp(1:end-1);

d = (1.0/2.0)*sqrt(sum(K(2,:) - K(1,:).^2));

K_1 = K(1,:);
K_2 = K(2,:);
K_3 = K(3,:);
M = (K(1, :) + K(2, :))/2;
M_23 = (K(3, :) + K(2, :))/2;

rot120 = [cos(2*pi/3), sin(2*pi/3); -sin(2*pi/3), cos(2*pi/3)];

q1_12 = K_1'-K_2';
q2_12 = rot120 * q1_12;
q3_12 = rot120 * q2_12; 

q1_23 = K_2'-K_3';
q2_23 = rot120 * q1_23;
q3_23 = rot120 * q2_23; 

if q_cut_type == 1  % K-Gamma-M of the L12 supercell
    k_sc = q3_12;
    m_sc = 0.5*(q3_12-q2_12);
    gamma_sc = [0, 0];
    
    pt(1,:) = k_sc;
    pt(2,:) = gamma_sc;
    pt(3,:) = m_sc;
    pt(4,:) = k_sc;
    type = 'L12 supercell';
    xt_labels = {'$K_{12}$', '$\Gamma_{12}$', '$M_{12}$', '$K_{12}$'};
    
elseif q_cut_type == 2  % All 3 cones
    pt(1,:) = -q1_12;
    pt(2,:) = [0, 0];
    pt(3,:) = q1_23;
    type = '3 cones';
    xt_labels = {'$K_{L1}$', '$K_{L2}$', '$K_{L3}$'};
    
elseif q_cut_type == 3  % K-Gamma-M of the L23 supercell
    
    k_sc = -q2_23;
    m_sc = 0.5*(q3_23-q2_23);
    gamma_sc = [0, 0];

    pt(1,:) = k_sc;
    pt(2,:) = gamma_sc;
    pt(3,:) = m_sc;
    pt(4,:) = k_sc;
    
    type = 'L23 supercell';
    
    xt_labels = {'$K_{23}$', '$\Gamma_{23}$', '$M_{23}$', '$K_{23}$'};
    
elseif q_cut_type == 4 % if the two twisting angles are the same, approximate the reciprocal lattice of the moire cell
    k_sc = 1/3*(2*b_tri(:, 1) + b_tri(:, 2));
    gamma_sc = zeros(size(k_sc)); 
    m_sc = 0.5*b_tri(:,1); 
    
    pt(1,:) = gamma_sc;
    pt(2,:) = m_sc;
    pt(3,:) = k_sc;
    pt(4,:) = gamma_sc;
    type = 'Trilayer supercell';
    
    xt_labels = {'$\bar{\Gamma}$', '$\bar{M}$', '$\bar{K}$', '$\bar{\Gamma}$'};
    
elseif q_cut_type == 5
    pt(1,:) = -q1_12;
    pt(2,:) = [0, 0];
    pt(3,:) = q3_12; 
    pt(4,:) = -q2_23;
    pt(5,:) = [0, 0];
    pt(6,:) = q1_23;
     
    type = '';
    xt_labels = {'$K_{L1}$', '$K_{L2}$', '$\Gamma_{12}$', '$\Gamma_{23}$', '$K_{L2}$', '$K_{L3}$'};
    
end

max_seg = size(pt,1)-1;

% making the line segments through the selected points
for seg = 1:max_seg
    if seg == 1
        q_list_x = pt(1,1)*(1-samp) + pt(2,1)*samp;
        q_list_y = pt(1,2)*(1-samp) + pt(2,2)*samp;
    else
        q_list_x = [q_list_x; pt(seg,1)*(1-samp) + pt(seg+1,1)*samp];
        q_list_y = [q_list_y; pt(seg,2)*(1-samp) + pt(seg+1,2)*samp];
    end  
end

q_list = [q_list_x, q_list_y]; % the list of k points 

q_list(end+1, :) = [pt(max_seg+1,1), pt(max_seg+1,2)];

% if q_cut_type ~= 2 && q_cut_type ~= 4
%     q_list = [q_list; -q_list];
% end

% q_list = (rot120*q_list')';

if q_cut_type == 1
    q_list = q_list+q3_12'; 
elseif q_cut_type == 3
    q_list = q_list-q2_23'; 
elseif q_cut_type == 4 
    q_list = q_list;%+b23(:,2)';
end 

ni(1) = 1; 
for i = 2:max_seg+1 
    ni(i) = (i-1)*size(samp,1);
end 
ni(max_seg+1) = ni(max_seg+1)+1;

for p_idx = 1:max_seg
    dis_here = norm(pt(p_idx+1,:)-pt(p_idx,:));
    if p_idx == 1
        karr = linspace(0,dis_here,ni(p_idx+1)-ni(p_idx)+1);
    else 
        karr(ni(p_idx)+1:ni(p_idx+1)) = linspace(karr(length(karr))+dis_here/(ni(p_idx+1)-ni(p_idx)), ...
            karr(length(karr))+dis_here,ni(p_idx+1)-ni(p_idx));
    end 
end 

for i = 1:length(ni)
    xt(i) = karr(ni(i));
end 


%% 
% setup kp model (output the list of scattered k's)
dof_list = kDoF_tri(layers,k_cutoff,grid_search);

% generate structure
dof_list.gen_dof()

% get the kpoints
k_list = dof_list.k_list();
ndof = size(k_list,1);
fprintf("%d total k points \n",ndof)


layer_list = k_list(:, 3);
k1 = k_list(layer_list == 1, :);
k2 = k_list(layer_list == 2, :);
k3 = k_list(layer_list == 3, :);
sz = 20;

figure
hold all;
scatter(k1(:, 1), k1(:, 2), sz, 'filled', 'r');
scatter(k2(:, 1), k2(:, 2), sz, 'r');
scatter(k3(:, 1), k3(:, 2), sz, 'filled', 'b');
scatter(q_list(:, 1), q_list(:,2), sz, 'filled');
axis equal
% xlim(xl)
% ylim(yl)
legend('L1', 'L2', 'L3');
title(type)
box on; 

H_inter = gen_interlayer_terms_mbd(k_list,layers);

for k_idx = 1:size(k_list,1)
    for tar_sheet = 1:3 
        if min( k_list(k_idx,3:end) == [tar_sheet 0 0 0 0 0 0])
           tar_dofs(tar_sheet, :) = 2*k_idx + [-1 0] ; 
        end
    end 
end

%% 
% get monolayer terms for each k point, then the band structure
for q_idx = 1:size(q_list,1)
    id = mod(q_idx,size(q_list,1));
    if id == 0 
        id = size(q_list,1);
    end 

    if q_idx>size(q_list,1)
        tar_q = q_list(id,:).*(-1);
    else
        tar_q = q_list(id,:);
    end 
    
    H_intra = gen_intralayer_terms_dirac(k_list,layers,tar_q,E_field);
    % H_intra = gen_intralayer_terms(k_list,layers,tar_q);
    H = H_intra+H_inter; 
    
    if q_cut_type < 4
        sigma = 0;
    elseif q_cut_type >= 4
        sigma = 4e-3;
    end 
    
    [raw_vecs, raw_vals] = eigs(H,length(H)); 
    % [raw_vecs, raw_vals] = eigs(H,num_eigs,0); 
    % raw_vecs: eigenvectors with d = (2*n_dof, 2*n_dof), 1st index: dof, 2nd index: band
    [vals(q_idx,:), order] = sort(real(diag(raw_vals))); % sort eigenvalues from smallest to largest 
    for tar_sheet = 1:3
        weights(tar_sheet,q_idx,:) = sum(abs(raw_vecs(tar_dofs(tar_sheet),order)).^2,1);
    end 
    fprintf("Diagonalization done with %d / %d \n",q_idx,size(q_list,1));
end


% karr = karr * 2.8205;
[offset1,a]=min(abs(vals(1, :)));
%[offset2,a]=min(abs(vals(end/2+1, :)));


%%
colors = [.7, .2, 0; 
          0, .7, .2;
          .2, 0, .7];

% vals_cut = vals;

% if q_cut_type == 4 
%     for i = 1:size(vals,1)
%         vals_plus = vals(i,vals(i,:)>=0);
%         vals_minus = vals(i,vals(i,:)<0);
%         
%         [~,a] = sort(abs(vals_plus));
%         [~,b] = sort(abs(vals_minus), 'descend');
%         figure(2222)
%         hold on;
%         plot(vals_plus(b))
%         
%         vals_tmp(i,:) = [vals_minus(b), vals_plus(a)];
%     end 
%     
% end 

% figure('Position', [100 100 450 700]);
figure('Position', [66 343 785 462]);

% %clf
% % figure(123)
% box on;
% hold all

% weights = weights / max(weights(:));
% 
% for tar_sheet = 1:length(proj) 
%     scatter(karr(1),(vals(1, 1))*1e3,50,'o','MarkerFaceColor',colors(proj(tar_sheet), :),'MarkerEdgeColor','none')
% end 

% select certain bands to plot 
% if weight_cut 
%     vals_cut = [];
%     k = 1; 
%     for i = 1:size(vals,2)
%         weights_tmp = weights(:,:,i);
%         if max(weights_tmp(:))>wcut
%             vals_cut(:,k) = vals(:,i);
%             k = k + 1;
%         end 
%     end 
% else 
%     vals_cut = vals;
% end 

% if q_cut_type ~= 4
%     ylim([-200 200]);
% else 
%     ylim([-20 20]);
% end 
yl = ylim;

try 
    plot(karr, (vals(1:end/2,:))*1e3, 'k.', 'LineWidth', 0.8);
    plot(karr, (vals(end/2+1:end,:))*1e3, 'k.', 'LineWidth', 0.8);
catch 
    plot(karr, (vals(:,min(abs(vals)<max(yl))))*1e3, 'k.', 'MarkerSize',5);
end 

xlim([0 max(karr)]) 
yl = ylim;
ylim([-200 200])

% 
% for q_idx = 1:size(q_list,1)
%     id = mod(q_idx-1, length(karr))+1;
%         
%     for d = 1:size(vals,2)
%         for i = 1:length(proj)
%             tar_sheet = proj(i);
%             if (weights(tar_sheet,q_idx,d) > 0.025) && abs((vals(q_idx,d))*1e3) < max(ylim)
%                 scatter(karr(id),(vals(q_idx,d))*1e3, 50, 'o', 'MarkerFaceAlpha', min(weights(tar_sheet, q_idx,d),1),...
%                     'MarkerFaceColor',colors(tar_sheet, :),'MarkerEdgeColor','none')
%             end
%         end 
%     end
%     fprintf("Done with plotting %d / %d \n",q_idx,size(q_list,1));
% end
% 
% 
ylabel('Energy (meV)');
xticks(xt);
xticklabels(xt_labels); 
title(['$\theta_{12} = ' num2str(q12) '^\circ, \theta_{23} = ' num2str(q23) '^\circ$']);

% 
% if length(proj)>1 && q_cut_type ~= 4
%     lg = {size(proj)};
%     for i = 1:length(proj)
%         lg{i} = ['L' num2str(proj(i))];
%     end 
%     legend(lg);
% end 
% 
% if int_off 
%     title(['$\theta_{12} = \theta_{23} = ' num2str(abs(theta_list(3))) '^\circ$, $T^{12} (\vec{r})$ turned off'])
% else 
%     title(['$\theta_{12} = ' num2str(-theta_list(1)) '^\circ$, $\theta_{23} = ' num2str(abs(theta_list(3))) '^\circ$'])
% end 
% 
% if savepng 
%     saveas(gcf, ['./figures_prelim/' figname '.png'])
% end
%%
if savedata 
    save(['./data_new/q12_' num2str(abs(theta_list(1))) '_q23_' num2str(abs(theta_list(3)))...
                '_kcut_' num2str(k_cutoff) '_qtype_' num2str(q_cut_type) '_23intoff.mat']);
end 

toc
