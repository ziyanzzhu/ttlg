% author: ziyan (zoe) zhu
% email: zzhu1@g.harvard.edu
% generate the interlayer Hamiltonian
% inputs: k_list: momentum space basis; 
%         layers: the `layers` object created using Layer.m
function H_inter = gen_interlayer_terms_dft(k_list, layers, intercoupling, q_here, inter_q_cut)
    ndof = size(k_list,1);
    H_inter = zeros(2*ndof, 2*ndof); % factor of 2 for A,B sublattice
    G = layers(1).G;
    G1 = G(:,1);
    G2 = G(:,2);
    K0 = 1/3 * (2*G1+G2);
    
    % compute scattering pairs
    for i = 1%:2
        l_1 = i;
        l_2 = mod(l_1,2)+1;
        from_idx = find(k_list(:,5) == l_1);
        to_idx = find(k_list(:,5) == l_2);
        orbPosFrom = layers(l_1).orbPos;
        orbPosTo = layers(l_2).orbPos;
        
        for o1 = 1:2
            for o2 = 1:2
                tK_h = squeeze(intercoupling.tK(o1,o2,:));
                F{o1,o2} = scatteredInterpolant(intercoupling.K_meshx, intercoupling.K_meshy, tK_h);
            end
        end
        
        k1 = k_list(from_idx, 1:2);
        k2 = k_list(to_idx, 1:2); 
        q_diffx = k1(:,1)+k2(:,1)' + q_here(1);
        q_diffy = k1(:,2)+k2(:,2)' + q_here(2);
        from_mat = from_idx + zeros(size(to_idx))';
        to_mat = zeros(size(from_idx)) + to_idx';

        % find elements within cutoff radius
        id = (q_diffx.^2 + q_diffy.^2) < inter_q_cut^2;
        AA_idx_from = 2*(from_mat(id)-1)+1;
        AA_idx_to = 2*(to_mat(id)-1)+1;
        
        R120 = [cosd(120) -sind(120); sind(120) cosd(120)];
        sigma_z = [1 0; 0 -1];
        sigma_x = [0 1; 1 0];
        
        nR = 1; % 3
        nC2 = 1; % 2  
        nM = 1; % 2
            
        for R_sym = 1:nR % C3

            R120_h = R120^(R_sym-1);
            q_x_R = R120_h(1,1)*q_diffx(id) + R120_h(1,2)*q_diffy(id);
            q_y_R = R120_h(2,1)*q_diffx(id) + R120_h(2,2)*q_diffy(id);

            phase_mat_R = expm( 1j*2*pi/3*sigma_z)^(R_sym-1);

            for C2_sym = 1:nC2 % C2T
                phase_mat_C2 = sigma_x^(C2_sym-1);

                for M_sym = 1:nM % Mirror
                    phase_mat_M = sigma_x^(M_sym-1);
                    M_fac = (-1)^(M_sym-1);
                    % DOESNT WORK WITH CURRENT K-MESH! KEEP nM = 1
                    q_x = M_fac*q_x_R; % M flips x, but apply -M not M
                    q_y = q_y_R;
                    H_temp = zeros(2*ndof,2*ndof);
                    
                    for o1 = 1:2
                        for o2 = 1:2
                            d_h = orbPosTo(o1,:) - orbPosFrom(o2,:);
                            phase_h = exp( 1j*(q_x*d_h(1) + q_y*d_h(2)));
                            tk_h = F{o1,o2}(q_x,q_y);
                            H_temp = H_temp + sparse(AA_idx_to+o1-1, AA_idx_from+o2-1, tk_h.*phase_h, size(H_inter,1), size(H_inter,2));
                        end
                    end
                    
                    if (nR*nC2*nM > 1)
                        for idx_t = 1:length(AA_idx_to)
                            sm_t = AA_idx_to(idx_t)+[0,1];
                            sm_f = AA_idx_from(idx_t)+[0,1];
                            H_2b2 = H_temp(sm_t,sm_f); % get the 2x2 matrix sublocks
                            H_2b2 = phase_mat_R*H_2b2*conj(phase_mat_R)/nR; % C3
                            if (C2_sym == 2)
                                H_2b2 = conj(H_2b2);
                            end
                            H_2b2 = phase_mat_C2*H_2b2*phase_mat_C2/nC2; % C2T

                            if (M_sym == 2)
                               H_2b2 = H_2b2'; 
                            end
                            H_2b2 = phase_mat_M*H_2b2*phase_mat_M/nM; % Mirror

                            H_temp(sm_t,sm_f) = H_2b2; % reassign values
                        end
                    end
                    
                    H_inter = H_inter + H_temp;
                
                end
            end
        end
        
        %{
        for o1 = 1:2
            for o2 = 1:2
                d_h = orbPosTo(o1,:) - orbPosFrom(o2,:);
                tK_h = squeeze(intercoupling.tK(o1,o2,:));
                phase_h = exp(-1j*(q_diffx(id)*d_h(1) + q_diffy(id)*d_h(2)));
                F = scatteredInterpolant(intercoupling.K_meshx, intercoupling.K_meshy, tK_h);
                tk_h = F(q_diffx(id),q_diffy(id));
                H_temp = sparse(AA_idx_to+o2-1, AA_idx_from+o1-1, tk_h.*phase_h, size(H_inter,1), size(H_inter,2));
                H_inter = H_inter + H_temp;
            end
        end
        %}
            
                            
    % make the Hamiltonian Hermitian
    H_inter = H_inter + H_inter'; 
    end
    
end