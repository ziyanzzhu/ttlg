% author: ziyan (zoe) zhu
% email: zzhu1@g.harvard.edu
% calculate the density of states
% inputs: q1_list: a list of the 1st twist angles (theta_12) in deg.
%         q2_list: a list of the 2nd twist angles (theta_23) in deg.
%         E_field: vertical displacement field (energy in eV across 3
%         layers) 
%         k_cutoff: k space cutoff in the unit of reciprocal space lattice
%         constant 
%         param: Gaussian FWHM
%         nq: grid size. total grid size is (nq+1)x(nq+1)
%         num_eigs: number of eigenvalues to keep in the diagonalization
%         E_list: the list of energy to sample
%         q_cut_type: type of Brillouin zone sampling. 
%                     type 0: monolayer Brillouin zone
%                     type 1: L12 moire Brillouin zone
%                     type 2: L23 moire Brillouin zone
% return: dos_tot: total density of states (normalized) 
%         dos_tot_mono: total density of states of monolayer graphene
%         (normalized) 
%         dos_fit: polyfit to dos_tot_mono
%         prefac: normalization constant

function [dos_tot, dos_tot_mono, dos_fit, prefac]=dos_calc_tri(alpha,q1_list,q2_list,E_field,k_cutoff,param,nq,num_eigs,E_list,q_cut_type,save_data)

    linecut = 0; % if on, calculate the DOS for a line (useful for debugging)
    run_idx = 1; 
 
    for q1_idx = 1:length(q1_list) 
        for q2_idx = 1:length(q2_list)
                
            q12 = q1_list(q1_idx);
            q23 = q2_list(q2_idx);
            theta_list = [q12 0 q23];                         % twisting angles (global)
            grid_search = 15;                                 % [G1,G2] are in [-grid_search,grid_search]^2 
            disp("Global twist angles in deg.: ") 
            disp(theta_list)
   
            % create layer data structures
            for t = 1:3
               layers(t) = Layer(t,deg2rad(theta_list(t)),alpha);
            end

            G1 = layers(1).G; 
            G2 = layers(2).G;
            G3 = layers(3).G;
            b12 = G1 - G2; 
            b23 = G2 - G3; 

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
            
            if q_cut_type == 0     % monolayer Brillouin zone
                bz_base_1 = G2(:,1);
                bz_base_2 = G2(:,2);

                type = 'Monolayer';

            elseif q_cut_type == 1  % L12 moire Brillouin zone
                bz_base_1 = b12(:,1);
                bz_base_2 = b12(:,2);

                type = 'L12 supercell';

            elseif q_cut_type == 2  % L23 smoire Brillouin zone
                bz_base_1 = b23(:,1);
                bz_base_2 = b23(:,2);

                type = 'L23 supercell';

            end

            bz_n = nq;
            dk = 1/bz_n;
            grid_mesh = 0:dk:1;

            if linecut 
                [x_grid,y_grid] = meshgrid(grid_mesh,1);
            else 
                [x_grid,y_grid] = meshgrid(grid_mesh,grid_mesh);
            end 

            kx_grid = bz_base_1(1)*x_grid + bz_base_2(1)*y_grid; % shift the origin to be K
            ky_grid = bz_base_1(2)*x_grid + bz_base_2(2)*y_grid;

            q_list = [kx_grid(:) ky_grid(:)];

            % setup kp model (output the list of scattered k's)
            dof_list = kDoF_tri(layers,k_cutoff,grid_search);

            % generate structure
            dof_list.gen_dof()

            % get the kpoints
            k_list = dof_list.k_list();
            ndof = size(k_list,1);
            fprintf("%d total k points \n",ndof)

            % find index corresponding to the center site 
            for k_idx = 1:size(k_list,1)
                for tar_sheet = 1:3 
                    if min( k_list(k_idx,3:end) == [tar_sheet 0 0 0 0 0 0])
                       tar_dofs(tar_sheet, :) = 2*k_idx + [-1 0] ; 
                    end
                end 
            end

            % interlayer term
            H_inter = gen_interlayer_terms(k_list,layers);

            clear dos

            % get monolayer terms for each k point, then the band structure
            parfor  q_idx = 1:size(q_list,1)
                id = mod(q_idx,size(q_list,1));
                if id == 0 
                    id = size(q_list,1);
                end 

                tar_q = q_list(id,:);

                H_intra = gen_intralayer_terms_dirac(k_list,layers,tar_q,E_field);
                H = H_inter + H_intra;

                % monolayer DOS
                dos_mono(q_idx, :) = dos_gauss_smear(H_intra, tar_dofs, param, E_list, num_eigs, 1e-4);
                
                dos(q_idx, :) = dos_gauss_smear(H, tar_dofs, param, E_list, num_eigs, 1e-4);

                fprintf("Diagonalization done with %d / %d \n",q_idx,size(q_list,1));
            end
         
            dos_tot_mono = sum(dos_mono,1);

            % fit to expected slope 
            cond = E_list>=0 & E_list < 100e-3;
            dos_fit = polyfit(E_list(cond), dos_tot_mono(cond),1); 

            % normalize to graphene DOS, g(E) = 2*A_c*E/pi/vF^2
            hbar = 6.582119569e-16; % unit: eV.s
            vF = 0.8e15; % unit: nm/s 
            slope = 6/pi/vF^2/hbar^2; % expected slope for monolayer graphene DOS. factor of 3 from three layers 

            prefac(run_idx) = slope/dos_fit(1); 
            dos_tot_mono = prefac(run_idx) * dos_tot_mono;
                
            
            disp(prefac(run_idx))

            clear dos_tot dos_int
            
            dos_tot = sum(dos,1) * prefac(run_idx);
 
            fname = ['dos_q12_' num2str(abs(theta_list(1))) '_q23_' num2str(abs(theta_list(3)))...
                            '_kcut_' num2str(k_cutoff) '_qtype_' num2str(q_cut_type) '_nq_' num2str((nq+1)^2) '_zip.mat'];
            
            % save data
            if save_data
                save(['./data/' fname ], 'dos_tot', 'dos_tot_mono', 'param',...
                    'ndof', 'E_list', 'prefac', 'q12', 'q23','cond','dos_fit')
           
            end 

            run_idx = run_idx + 1; 
        end
    end 
