% authors: ziyan (zoe) zhu
% email: zzhu1@g.harvard.edu
% Example calculation of the tTLG DOS

clear all
f_size = 22;
set(groot, 'DefaultTextInterpreter', 'Latex')
set(groot, 'DefaultLegendInterpreter', 'Latex')
set(groot, 'DefaultAxesTickLabelInterpreter', 'Latex')
set(0,'DefaultAxesFontSize',f_size)

% the list of twist angles
q1_list = -1.5; 
q2_list = [2.5:0.3:3.3];

monoG = 1;  % if on, turn off the interlayer coupling
savedata = 0; 
savetxt = 0; 

k_cutoff = 2;        % cutoff in momentum space in the of reciprocal lattice constant
grid_search = 15;    % [G1,G2] are in [-grid_search,grid_search]^2 
num_eigs = 20;       % the number of eigenvalues to keep 
q_cut_type = 1;      % Brillouin zone sampling
E_field = 0;         % electric field (potential energy in eV across the 3 layers)
nq = 20;             % number of k points to sample

% the total number of twist angles to calculate
tot_pt = length(q1_list)*length(q2_list);
qidx_here = 0;

% loop over all twist angle combinations
for q1_idx = 1:length(q1_list) 
    for q2_idx = 1:length(q2_list)
        
        qidx_here = qidx_here + 1; 
        fprintf("Running twist angle %d/%d \n", qidx_here, tot_pt)
       
        theta_list = [q1_list(q1_idx) 0 q2_list(q2_idx)]; % twisting angles (global)

        % w_inv ~ 1/width of the Gaussian
        w_inv = 1e3/2;
        E_list = linspace(-1,1,1e3);

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

        if q_cut_type == 0  % monolayer 
            bz_base_1 = G2(:,1);
            bz_base_2 = G2(:,2);

        elseif q_cut_type == 1  % L12 supercell
            bz_base_1 = b12(:,1);
            bz_base_2 = b12(:,2);

            type = 'L12 supercell';

        elseif q_cut_type == 2  % L23 supercell

            bz_base_1 = b12(:,1);
            bz_base_2 = b12(:,2);

            type = 'L23 supercell';
        end

        for nq_idx = 1:length(nq) 
            bz_n = nq(nq_idx);
            dk = 1/bz_n;
            grid_mesh = 0:dk:1;
            [x_grid,y_grid] = meshgrid(grid_mesh,grid_mesh);
            kx_grid = bz_base_1(1)*x_grid + bz_base_2(1)*y_grid; % shift the origin to be K
            ky_grid = bz_base_1(2)*x_grid + bz_base_2(2)*y_grid;

            q_list = [kx_grid(:) ky_grid(:)];

        %     if q_cut_type == 4
        %         q_list = q_list+b23(:,2)';
        %     end 

            % figure;
            % scatter(kx_grid(:), ky_grid(:))
            % axis equal;

            %% 
            % 
            % figure(123)
            % set(gcf,'Position',[-15 521 884 284])
            % hold on;
            
            % loop over the list of 
            for c_idx = 1:length(k_cutoff)
                tic
                % setup kp model (output the list of scattered k's)
                dof_list = kDoF_tri(layers,k_cutoff(c_idx),grid_search);

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
                % 
                % figure
                % hold all;
                % scatter(k1(:, 1), k1(:, 2), sz, 'filled', 'r');
                % scatter(k2(:, 1), k2(:, 2), sz, 'r');
                % scatter(k3(:, 1), k3(:, 2), sz, 'filled', 'b');
                % scatter(q_list(:, 1), q_list(:,2), sz);
                % axis equal
                % % xlim(xl)
                % % ylim(yl)
                % legend('L1', 'L2', 'L3');
                % title(type)
                % box on; 

                H_inter = gen_interlayer_terms_mbd(k_list,layers);

                for k_idx = 1:size(k_list,1)
                    for tar_sheet = 1:3 
                        if min( k_list(k_idx,3:end) == [tar_sheet 0 0 0 0 0 0])
                           tar_dofs(tar_sheet, :) = 2*k_idx + [-1 0] ; 
                        end
                    end 
                end

                clear dos
                
                % get monolayer terms for each k point
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
                   
                    if monoG == 0 
                        H = H_intra+H_inter; 
                    else 
                        H = H_intra;
                    end 

                    dos(q_idx, :) = dos_gauss_smear(H, tar_dofs, w_inv, E_list, length(H), sigma);

                    fprintf("Diagonalization done with %d / %d \n",q_idx,size(q_list,1));
                end

                % normalization: making sure that the dos integrate to
                % unity 
                dos_tot = sum(dos,1);
                % integrated dos (unit: eV) 
                dos_int = sum(dos_tot)*mean(diff(E_list));
                int_norm = 1/dos_int; 
                dos_tot = int_norm*dos_tot; 
                              
                % fit the slope 
                dos_fit = polyfit(E_list(E_list>=0 & E_list < 350e-3), ...
                    dos_tot(E_list>=0 & E_list < 350e-3),1); 
                
                % normalize to graphene DOS, g(E) = 2*A_c*E/pi/vF^2
                hbar = 6.582119569e-16; % unit: eV.s
                vF = 0.8e16; % unit: A^2/s 
                slope = 6/pi/vF^2/hbar^2; % factor of 3 from three layers 
                
                prefac(nq_idx) = slope/dos_fit(1); 
                dos_tot = prefac * dos_tot;
                
                if monoG
                    disp(prefac(nq_idx))
                end
                
                a0 = 1.42*sqrt(3);
                A0 = sqrt(3)/2*a0^2;
                t = 2*6.582*0.8/3/a0*1e3;

                % dos_tot = dos_tot/length(q_list);
                % dos_tot = 2*dos_tot/max(dos_tot)/length(q_list); % normalize by the number of k sampling points     
%                 a0 = 1.42*sqrt(3);
%                 Asc = sqrt(3)/2*(a0*0.1)^2/(4*sin(deg2rad(theta_list(3))/2)^4); % trilayer supercell area in nm^2
%                 ne = 8/Asc; % density DOS should integrate to 
%                 prefac(nq_idx) = ne/dos_int; % unit: nm^(-2)*eV^(-1)
%                 disp(prefac(nq_idx))
                
%                 figure
%                 hold on;
%                 plot(dos_tot,E_list*1e3-4);
%                 plot(E_list*slope+dos_fit(2)*prefac, E_list*1e3-4)
%                 ylabel('Energy (meV)')
%                 xlabel('DoS $\mathrm{(eV^{-1}\cdot\AA^{-2})}$');
%                 box on; 
%                 xlim([0 max(dos_tot)*1.5]);
%                 ylim([-500 500]);
                
% 
%                 dos = reshape(dos, [nq(nq_idx)+1, nq(nq_idx)+1, length(E_list)]);
%                 rot120 = [cos(2*pi/3), -sin(2*pi/3); 
%                         sin(2*pi/3), cos(2*pi/3)];
%                 for i = 1:10
%                     E_idx = length(E_list)/2-50+10*i;
%                     figure 
%                     k_tmp = [kx_grid(:), ky_grid(:)]';
%                     hold on
%                     for j = 1:3
%                         k_tmp = rot120*k_tmp;
%                         kx_here = reshape(k_tmp(1, :), size(kx_grid));
%                         ky_here = reshape(k_tmp(2, :), size(ky_grid));
%                         surf(kx_here, ky_here, squeeze(log(1e-20+dos(:, :, E_idx))), 'EdgeColor', 'none');
%                     end 
%                     axis equal
%                     colormap copper
%                     view(2)
%                     caxis([-20 0])
%                     shading interp
%                     title(['E = ' num2str(E_list(E_idx)*1e3) ' meV'])
%                     yticks([])
%                     xticks([])
%                 end 

                
%               title(['$\theta_{12} = ' num2str(-theta_list(1)) '^\circ$, $\theta_{23} = ' num2str(abs(theta_list(3))) '^\circ$'])

% 
%                 if savepng 
%                     saveas(gcf, ['./figures_prelim/' figname '.png'])
%                 end
                %%
                if savedata
                    save(['./data_dos_gauss/dos_q12_' num2str(abs(theta_list(1))) '_q23_' num2str(abs(theta_list(3)))...
                                '_kcut_' num2str(k_cutoff) '_qtype_' num2str(q_cut_type) '_nq_' num2str((nq+1)^2) '.mat']);
                end 
                
                if savetxt 
                    fname_txt = ['./data_dos_gauss/dos_q12_' num2str(abs(theta_list(1))) '_q23_' num2str(abs(theta_list(3)))...
                                '_kcut_' num2str(k_cutoff) '_qtype_' num2str(q_cut_type) '_nq_' num2str((nq+1)^2) '.txt'];
                    % Create a table with the data and variable names
                    data = [E_list(:)'; dos_tot(:)'];
                    % Write data to text file
                    fileID = fopen(fname_txt,'w');
                    fprintf(fileID,'%5s %3s\n','Elist','DOS');
                    for d_idx = 1:length(data)
                        fprintf(fileID,'%12.8f %12.8f\n',data(:,d_idx));
                    end 
                    fclose(fileID);
                end 
                
                fprintf("k cutoff done with %d / %d \n",c_idx,size(k_cutoff,1));
                toc
            end 
        end
    
        
        
    end 
end 

disp(prefac)