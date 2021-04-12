% authors: ziyan (zoe) zhu, stephen carr 
% email: zzhu1@g.harvard.edu
% Example calculation of the TBG band structure
clear all
f_size = 12;
set(groot, 'DefaultTextInterpreter', 'Latex')
set(groot, 'DefaultLegendInterpreter', 'Latex')
set(groot, 'DefaultAxesTickLabelInterpreter', 'Latex')
set(0,'DefaultAxesFontSize',f_size)

% STILL TO DO:
% 1: symm breaking (~2 meV) by relaxations, check basis for configuration space...
% 2: add intralayer strain corrections (should be easy... ?)

% !! Settings

% geometry 
theta = 1.15;    % relative twist angle (in degrees)
relax_str = 1.0; % scaling of atomic relaxations (0.0 turns off relaxations)
E_field = 0.0;   % vertical displacement field in eV (total potential energy across the three layers)

% truncation of momentum basis
k_cutoff = 4;           % momentum cutoff, in units of norm(G1)
grid_search = 20;       % [G1,G2] are in [-grid_search,grid_search]^2 

% interlayer coupling sampling mesh
r_max = 8; % maximum radius, in Angstroms
dr = .25; % spacing of r mesh, in Angstroms
% interlayer fourier transform control
inter_q_cut_scale = 5; % maximum scattered momentum, in units of K0
inner_k_rad_scale = 5; % radius of each "island", in units of b12
dk_scale = 1/2; %  spacing of k mesh, in units of b12

% band structure line cut settings
q_cut_type = 1;         % what kind of line-cut we do in momentum space                        % 1: high symmetry line in the L12 bilayer moire Brillouin zone (single valley)
nq = 50;                % number of k points to sample on each high symmetry line segment 

% turn on/off saving
savedata = 0;           % save useful variables to folder ./data/

theta_list = theta*[-0.5, 0.5];  % twist of each layer
% create layer data structures
for t = 1:2
   layers(t) = Layer(t,deg2rad(theta_list(t)));
end

% generate orbital positions
nsheets = 2;
norbs = 2;
orb_pos = zeros(nsheets,norbs,2);
A1 = layers(1).A;
A2 = layers(2).A;
orb_pos(1,:,:) = layers(1).orbPos;
orb_pos(2,:,:) = layers(2).orbPos;

% reciprocal geometry
G1 = layers(1).G; 
G2 = layers(2).G;
b12 = G2 - G1; % moire reciprocal vectors

G11 = G1(:, 1);
G12 = G1(:, 2);
G13 = G11 + G12; 
K0 = 1/3 * (G11 + G13); % a K point of Layer 1

inter_q_cut = inter_q_cut_scale*norm(K0); % maximum scattered momentum
inner_k_rad = inner_k_rad_scale*norm(b12(:,1)); % size of scattering "island"
dk = dk_scale*norm(b12(:,1)); % spacing of k-mesh  


% get the K-point of each monolayer
for t = 1:2
   th = layers(t).theta - layers(1).theta;
   K(t,:) = [cos(th) -sin(th); sin(th) cos(th)]*K0;
end

% define the k-point sampling
samp = linspace(0,1,nq)';
samp = samp(1:end-1);

d = (1.0/2.0)*sqrt(sum(K(2,:) - K(1,:).^2)); 

K_1 = K(1,:);
K_2 = K(2,:);
M = (K(1, :) + K(2, :))/2; % M point of supercell

rot120 = [cos(2*pi/3), sin(2*pi/3); -sin(2*pi/3), cos(2*pi/3)];

% definition of q vectors (NN coupling separations in k space)
% these are used in the BMD model, not in the DFT model.
q1_12 = K_2'-K_1';
q2_12 = rot120 * q1_12;
q3_12 = rot120 * q2_12; 

% defining high symmetry points
switch q_cut_type 
    case 1 % K-Gamma-M of the L12 supercell
        k_sc = q3_12;
        m_sc = 0.5*(q3_12-q2_12);
        gamma_sc = [0, 0];

        pt(1,:) = k_sc;
        pt(2,:) = gamma_sc;
        pt(3,:) = m_sc;
        pt(4,:) = k_sc;
        type = 'L12 supercell';
        xt_labels = {'$K_{12}$', '$\Gamma_{12}$', '$M_{12}$', '$K_{12}$'};
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

q_list = q_list - k_sc';
% calculate the path length at q point 
ni(1) = 1; 
for i = 2:max_seg+1 
    ni(i) = (i-1)*size(samp,1)+1;
end 
%ni(max_seg+1) = ni(max_seg+1)+1;

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

% Generate momentum lattice
% setup kp model (output the list of scattered k's)
dof_list = kDoF_bi(layers,k_cutoff,grid_search);

% generate structure
dof_list.gen_dof()

% get the kpoints
k_list = dof_list.k_list();
ndof = size(k_list,1);
fprintf("%d total k points \n",ndof)

% find the indices corresponding to the center site (of K lattice, not H)
layer_list = k_list(:, 5);
k1 = k_list(layer_list == 1, :);
k2 = k_list(layer_list == 2, :);
for k_idx = 1:size(k_list,1)
    for tar_sheet = 1:2
        if min( k_list(k_idx,5:end) == [tar_sheet 0 0 0 0])
           tar_dofs(tar_sheet, :) = k_idx;
        end
    end 
end


% generate interlayercouplings object (used for interlayer Hamiltonian)

% create R and K meshes
intercoupling = InterCouplings(layers, 1, 2, r_max, dr, inter_q_cut, dk, inner_k_rad, K0, theta, orb_pos);

% apply relaxations
tic
fprintf("Applying atomic relaxations... \n")
intercoupling.relax_configs(relax_str);
fprintf("Relaxations done: ")
toc

% compute realspace couplings
intercoupling.gen_tR();

% apply fourier transform
tic
fprintf("FT of interlayer coupling starting... \n")
intercoupling.gen_tK();
fprintf("FT done: ")
toc

% check average AA vs AB coupling on first shell
H_inter_K0 = gen_interlayer_terms_dft(k_list(tar_dofs,:),layers,intercoupling,K0,inter_q_cut);
w0 = abs(H_inter_K0(3,1));
w1 = abs(H_inter_K0(4,1));
fprintf("w0 (inter_AA) = %s meV \n",num2str(1000*w0,'%.0f'));
fprintf("w1 (inter_AB) = %s meV \n",num2str(1000*w1,'%.0f'));
%%
% get H for each k point, and compute the band structure
for q_idx = 1:size(q_list,1)
    
    tar_q = q_list(q_idx,:)+K0';
    %tar_q_bmd= q_list(q_idx,:);

    H_inter = gen_interlayer_terms_dft(k_list,layers,intercoupling,tar_q,inter_q_cut);
    %H_inter = H_inter_bmd;

    H_intra = gen_intralayer_terms_dft(k_list,layers,tar_q,E_field);
    %H_intra = gen_intralayer_terms_dirac(k_list,layers,tar_q_bmd,E_field);

    H = H_intra+H_inter; 
  
    num_eigs = size(H,1);
    [raw_vecs, raw_vals] = eigs(H,num_eigs); 
    % raw_vecs: eigenvectors with d = (2*n_dof, 2*n_dof), 1st index: dof, 2nd index: band
    [vals(q_idx,:), order] = sort(real(diag(raw_vals))); % sort eigenvalues from smallest to largest 
    for tar_sheet = 1:2
        %weights(tar_sheet,q_idx,:) = sum(abs(raw_vecs(tar_dofs(tar_sheet),order)).^2,1);
    end 
    if (mod(q_idx,10) == 0 || q_idx == 1)
        fprintf("H Diag: %d / %d \n",q_idx,size(q_list,1));
    end
end

%% plot the band structure
% eye pleasing colors
b = [0, 40, 143]/255;
r = [224, 40, 0]/255;
g = [15, 255-80, 103]/255;

%weights = weights / max(weights(:));
%wcut = 1e-7; % do not plot is the weight is less than this value 
%interval = 1; % whether to skip points when plotting the bands

%figure()
clf
if q_cut_type ~= 5 && q_cut_type ~= 2
    set(gcf, 'Position', [66 343 390 462]);
else 
    set(gcf, 'Position', [66 343 785 462]);
end 
box on
hold all

plot(qarr,1e3*vals,'-k')
theta_scale = 200*theta/1.15;
axis([qarr(1) qarr(end) -theta_scale theta_scale])

ylabel('Energy (meV)');
xticks(xt);
xticklabels(xt_labels);
xlim([min(qarr) max(qarr)])
title(['$\theta = ' num2str(abs(theta_list(2)-theta_list(1))) '^\circ$']);

%%
if savedata 
    save(['./data/q12_' num2str(abs(theta_list(1))) '_q23_' num2str(abs(theta_list(3)))...
                '_kcut_' num2str(k_cutoff) '_qtype_' num2str(q_cut_type) '_bands.mat'], ...
                'theta_list','E_field','weights', 'q_list', 'vals', 'qarr', 'xt', ...
                'xt_labels', 'q_cut_type');
end 

%toc
