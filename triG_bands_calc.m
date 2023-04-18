
% authors: ziyan (zoe) zhu, stephen carr 
% email: zzhu1@g.harvard.edu
% Example calculation of the tTLG band structure
clear all
f_size = 22;
set(groot, 'DefaultTextInterpreter', 'Latex')
set(groot, 'DefaultLegendInterpreter', 'Latex')
set(groot, 'DefaultAxesTickLabelInterpreter', 'Latex')
set(0,'DefaultAxesFontSize',f_size)

% setup & geometry 
theta_list = [-1.1 0 1.6];  % twisting angles in degree (global)
k_cutoff = 4;           % cutoff
grid_search = 20;       % [G1,G2] are in [-grid_search,grid_search]^2 

proj = [1, 2, 3];       % sheet to project eigenvector weights onto
q_cut_type = 5;         % what kind of line-cut we do in momentum space
                        % 1: high symmetry line in the L12 bilayer moire Brillouin zone (single valley)
                        % 2: connecting the Dirac points of 3 layers 
                        % 3: high symmetry line in the L23 bilayer moire Brillouin zone (single valley)
                        % 4: high symmetry line in the trilayer moire of moire Brillouin zone (only works when twist angles are equal and not very useful) 
                        % 5: high symmetry line along K_L1 -> K_L2 -> Gamma_12 -> Gamma_23 -> K_L2 -> K_L3
savedata = 1;           % save useful variables to folder ./data/
E_field = 0.0;          % vertical displacement field in eV (total potential energy across the three layers)
num_eigs = 100;          % the number of eigenvalues to keep in the diagonalization (near 0 energy)
nq = 40;                % number of k points to sample on each high symmetry line segment 
color_on = 1;           % plot colors in the band structure (wavefunction weights) 
alpha = 1.43*sqrt(3);

% create layer data structures
for t = 1:3
   layers(t) = Layer(t,deg2rad(theta_list(t)),alpha);
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

% definition of q vectors (NN coupling separations in k space)
q1_12 = K_1'-K_2';
q2_12 = rot120 * q1_12;
q3_12 = rot120 * q2_12; 

q1_23 = K_2'-K_3';
q2_23 = rot120 * q1_23;
q3_23 = rot120 * q2_23; 

% defining high symmetry points
switch q_cut_type 
    case 1 % K-Gamma-M of the L12 supercell (note: this is shifted because the origin is not gamma point)
        k_sc = q3_12;
        m_sc = 0.5*(q3_12-q2_12);
        gamma_sc = [0, 0];

        pt(1,:) = k_sc;
        pt(2,:) = gamma_sc;
        pt(3,:) = m_sc;
        pt(4,:) = k_sc;
        type = 'L12 supercell';
        xt_labels = {'$K_{12}$', '$\Gamma_{12}$', '$M_{12}$', '$K_{12}$'};
        
    case 2 % All 3 cones 
        pt(1,:) = -q1_12;
        pt(2,:) = [0, 0];
        pt(3,:) = q1_23;
        pt(4,:) = -q1_12;
        type = '3 cones';
        xt_labels = {'$K_{L1}$', '$K_{L2}$', '$K_{L3}$', '$K_{L1}$'};
    
    case 3  % K-Gamma-M of the L23 supercell
    
        k_sc = -q2_23;
        m_sc = 0.5*(q3_23-q2_23);
        gamma_sc = [0, 0];

        pt(1,:) = k_sc;
        pt(2,:) = gamma_sc;
        pt(3,:) = m_sc;
        pt(4,:) = k_sc;
        type = 'L23 supercell';
        xt_labels = {'$K_{23}$', '$\Gamma_{23}$', '$M_{23}$', '$K_{23}$'};
    
    case 4 
        % if the two twisting angles are the same, approximate the reciprocal lattice of the moire cell
        % not very useful
        k_sc = 1/3*(2*b_tri(:, 1) + b_tri(:, 2));
        gamma_sc = zeros(size(k_sc)); 
        m_sc = 0.5*b_tri(:,1); 

        pt(1,:) = gamma_sc;
        pt(2,:) = m_sc;
        pt(3,:) = k_sc;
        pt(4,:) = gamma_sc;
        type = 'Trilayer supercell';
    
        xt_labels = {'$\bar{\Gamma}$', '$\bar{M}$', '$\bar{K}$', '$\bar{\Gamma}$'};
    
    case 5 % long line cut
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

q_list = [q_list_x, q_list_y]; % the list of center sites

q_list(end+1, :) = [pt(max_seg+1,1), pt(max_seg+1,2)];

% shift the origin of the line cut to be K_L2 
if q_cut_type == 1
    q_list = q_list+q3_12'; 
elseif q_cut_type == 3
    q_list = q_list-q2_23'; 
end 

% calculate the path length at q point 
ni(1) = 1; 
for i = 2:max_seg+1 
    ni(i) = (i-1)*size(samp,1);
end 
ni(max_seg+1) = ni(max_seg+1)+1;

for p_idx = 1:max_seg
    dis_here = norm(pt(p_idx+1,:)-pt(p_idx,:));
    if p_idx == 1
        qarr = linspace(0,dis_here,ni(p_idx+1)-ni(p_idx)+1);
    else 
        qarr(ni(p_idx)+1:ni(p_idx+1)) = linspace(qarr(length(qarr))+dis_here/(ni(p_idx+1)-ni(p_idx)), ...
            qarr(length(qarr))+dis_here,ni(p_idx+1)-ni(p_idx));
    end 
end 

for i = 1:length(ni)
    xt(i) = qarr(ni(i));
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

% plot the generated momentum degrees of freedom / basis
figure;
hold all;
scatter(k1(:, 1), k1(:, 2), sz, 'filled', 'r');
scatter(k2(:, 1), k2(:, 2), sz, 'r');
scatter(k3(:, 1), k3(:, 2), sz, 'filled', 'b');
scatter(q_list(:, 1), q_list(:,2), sz, 'filled');
axis equal
legend('L1', 'L2', 'L3');
title(type)
box on; 

% find the indices corresponding to the center site (A and B sublattices)
for k_idx = 1:size(k_list,1)
    for tar_sheet = 1:3 
        if min( k_list(k_idx,3:end) == [tar_sheet 0 0 0 0 0 0])
           tar_dofs(tar_sheet, :) = 2*k_idx + [-1 0] ; 
        end
    end 
end



%% 

tic 
% construction of interlayer Hamiltonian
H_inter = gen_interlayer_terms(k_list,layers);

% get monolayer terms for each k point, then the band structure
for q_idx = 1:size(q_list,1)
    id = mod(q_idx,size(q_list,1));
    if id == 0 
        id = size(q_list,1);
    end 
    
    tar_q = q_list(id,:);
    
    H_intra = gen_intralayer_terms_dirac(k_list,layers,tar_q,E_field);
    H = H_intra+H_inter; 
  
    [raw_vecs, raw_vals] = eigs(H,num_eigs,4e-3); 
    % raw_vecs: eigenvectors with d = (2*n_dof, 2*n_dof), 1st index: dof, 2nd index: band
    [vals(q_idx,:), order] = sort(real(diag(raw_vals))); % sort eigenvalues from smallest to largest 
    if color_on 
        for tar_sheet = 1:3
            weights(tar_sheet,q_idx,:) = sum(abs(raw_vecs(tar_dofs(tar_sheet),order)).^2,1);
        end 
    end 
    fprintf("Diagonalization done with %d / %d \n",q_idx,size(q_list,1));
end

%% plot the band structure
% eye pleasing colors 
b = [0, 40, 143]/255;
r = [224, 40, 0]/255;
g = [15, 255-80, 103]/255;


if color_on 
    weights = weights / max(weights(:));
end 
wcut = 1e-7; % do not plot is the weight is less than this value 
interval = 1; % whether to skip points when plotting the bands



figure(123)
clf
if q_cut_type ~= 5 && q_cut_type ~= 2
    set(gcf, 'Position', [66 343 390 462]);
else 
    set(gcf, 'Position', [66 343 785 462]);
end 
box on
hold all

if q_cut_type ~= 4
    ylim([-180 180]);
else 
    ylim([-50 50]);
end 

if color_on 
    for k_idx = 1:size(q_list,1)/interval

        q_idx = (k_idx-1)*interval+1; 

        for d = 1:size(vals,2)

            color_here = [0,0,0];
            wcond = max(weights(:,q_idx,d))>wcut;

            if wcond && abs((vals(q_idx,d))*1e3) < max(ylim)
                alpha = 1;
                weights_here = weights(:,q_idx,d)/max(weights(:,q_idx,d));
                color_tmp = log(weights_here+1)/max(log(weights_here+1));
                color_here = weights_here(1)*r + weights_here(2)*b + weights_here(3)*g;

                if max(weights(:,q_idx,d)) < 0.025
                    alpha = sum(color_here);
                end 

                if max(weights(:,q_idx,d)) < 0.005
                    alpha = min(color_here)*2;
                end 
                if alpha > 1 
                   alpha = 1;
                end 
                scatter(qarr(q_idx),(vals(q_idx,d))*1e3, 15, 'o', 'MarkerFaceColor', color_here, ...
                    'MarkerEdgeColor','none', 'MarkerFaceAlpha', alpha)
            end
        end
        fprintf("Done with plotting %d / %d \n",q_idx,size(q_list,1));
    end
    
else 
    for d = 1:size(vals,2)
       %  plot(qarr,(vals(:,d))*1e3, 'k', 'LineWidth', 2);
        scatter(qarr,(vals(:,d))*1e3, 15, 'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor','none')
    end 
    
end

ylabel('Energy (meV)');
xticks(xt);
xticklabels(xt_labels);
xlim([min(qarr) max(qarr)])
title(['$\theta_{12} = ' num2str(abs(theta_list(1))) '^\circ, \theta_{23} = ' ...
    num2str(abs(theta_list(3))) '^\circ$']);
exportgraphics(gcf, './figures/bands_tswg.pdf', 'ContentType','vector');

%%
if savedata 
    save(['./data/q12_' num2str(abs(theta_list(1))) '_q23_' num2str(abs(theta_list(3)))...
                '_kcut_' num2str(k_cutoff) '_qtype_' num2str(q_cut_type) '_bands.mat'], ...
                'theta_list','E_field','weights', 'q_list', 'vals', 'qarr', 'xt', ...
                'xt_labels', 'q_cut_type');
end 

toc
