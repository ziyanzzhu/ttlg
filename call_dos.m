% author: ziyan (zoe) zhu
% email: zzhu1@g.harvard.edu
% example calculatiion of the DOS
clear all
q1_list = 1.8;                    % list of \theta_{12} in deg.
q2_list = [-2.8];              % list of \theta_{23} in deg.
k_cutoff = 3;                     % k space cutoff in the unit of reciprocal lattice constant   
E_list = linspace(-0.5,0.5,1e3);  % list of energies in eV
q_cut_type = 1;                   % type of Brillouin zone sampling. 
                                  % type 1: monolayer Brillouin zone
                                  % type 2: L12 moire Brillouin zone
                                  % type 3: L23 moire Brillouin zone
num_eigs = 80;                    % number of eigenvalues to keep in the diagonalization
nq = 31;                          % grid size
E_field = linspace(0.0,0.05,6);                      % vertical displacement field 

% parallelization; define the number of parallel processes
% if running on a cluster 
% np = str2num(getenv('SLURM_CPUS_PER_TASK')); % number of workers

% if running locally
np = 4;

k = 1;
for i = 1:length(q1_list)
    for j = 1:length(q2_list)
        twist(:, k) = [q1_list(i), q2_list(j)];
        k = k + 1;
    end 
end 

% total number of twist angles to calculate
tot_pt = length(q1_list)*length(q2_list)*length(E_field); 

figure(234)
set(gcf,'Position',[211 101 438 669])
clf

% if exist('np','var')
%     parpool('local',np)
% else 
%     parpool('local',1)
% end 
for i = 1:tot_pt
   
    % adapt the Gaussian width based on twist angles (parameters used in
    % the paper)
    if q_cut_type > 1
        if twist(q_cut_type-1,i) < 2 
            w_inv = 1e4/3;
        elseif twist(q_cut_type-1,i) <= 3.9
            w_inv = 1e3;
        else 
            w_inv = 500;
        end 
    else 
        w_inv = 1e3;
    end 
    
    w_inv = 100;
    if length(E_field) > 1
        lg{i} = ['$D = ' num2str(E_field(i)*1e3) '\, \mathrm{meV}$'];
    else
        lg{i} = ['$\theta_{12} = ' num2str(twist(1,i)) '^\circ$' 10 '$\theta_{23} = ' ...
            num2str(twist(2,i)) '^\circ$'];
    end 

    fprintf("Running twist angle %d/%d \n", i, tot_pt)
    tic
    if length(E_field) > 1
        fname = dos_calc_tri(twist(1),twist(2),E_field(i),k_cutoff,w_inv,nq,num_eigs,E_list,q_cut_type);
    else 
        fname = dos_calc_tri(twist(1,i),twist(2,i),E_field,k_cutoff,w_inv,nq,num_eigs,E_list,q_cut_type);
    end
    load(['./data/' fname])
    toc
    
    % check monolayer fit 
    subplot(1,2,1)
    box on
    hold all; 
    plot(dos_tot_mono, E_list-E_field(i)*2, 'LineWidth', 2)
%     plot(dos_tot_mono(cond), E_list(cond))
%     plot(polyval(dos_fit,E_list),E_list,'k--','LineWidth',2)
    ylim([-0.2 0.2])
%     xlim([0 max(dos_tot_mono(:))*1.1]);
    title('Monolayer')
    ylabel('Energy (eV')
    xlabel('DoS $\mathrm{(eV^{-1}\cdot\AA^{-2})}$');
    
    subplot(1,2,2)
    box on;
    hold all
    plot(dos_tot,E_list, 'LineWidth', 2)
    ylim([-0.2 0.2])
    yticklabels([])
    legend(lg)
    title('Full DOS')
    xlabel('DoS $\mathrm{(eV^{-1}\cdot\AA^{-2})}$');
    
    disp("=====================================")
end 

delete(gcp('nocreate'))


