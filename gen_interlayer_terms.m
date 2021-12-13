% author: ziyan (zoe) zhu
% email: zzhu1@g.harvard.edu
% generate the interlayer Hamiltonian
% inputs: k_list: momentum space basis; 
%         layers: the `layers` object created using Layer.m
function H_inter = gen_interlayer_terms(k_list, layers)
    ndof = size(k_list,1);
    H_inter = sparse(2*ndof, 2*ndof); % factor of 2 for A,B sublattice
    G = layers(2).G;
    G1 = G(:,1);
    G2 = G(:,2);
    K0 = 1/3 * (2*G1+G2);
    rot120 = [cos(2*pi/3), sin(2*pi/3); -sin(2*pi/3), cos(2*pi/3)];
    
    for i = 1:3 
        th = layers(i).theta;
        K(:,i) = [cos(th) -sin(th); sin(th) cos(th)]*K0;
    end 
    
    ph = [0, -2*pi/3, 2*pi/3];
    
    % interlayer coupling strengths of AA and AB-type
    w_aa = 0.07; 
    w_ab = 0.11; 
    
    % defining T matrices 
    for i = 1:length(ph)
        Tmat(:, :, i) = w_aa .* [1, 0; 0, 1] + w_ab .* [0, exp(1j*ph(i)); conj(exp(1j*(ph(i)))), 0];
    end 

    tol = 1e-5;
    
    % running over two pairs scattering pairs
    for i = 1:2
        l_1 = i;
        l_2 = mod(l_1,3)+1;
        l_3 = mod(l_1+1,3)+1;
        from_idx = find(k_list(:,3)==l_1);
        to_idx = find(k_list(:,3)==l_2);
        
        q1_12 = K(:,l_2)-K(:,l_3);
        q1_12 = K(:,l_1)-K(:,l_2);
        q2_12 = rot120 * q1_12;
        q3_12 = rot120 * q2_12;
        q_scatt = [q1_12, q2_12, q3_12];
%         disp(q_scatt)

        for x1 = 1:length(from_idx)
            idx_1 = from_idx(x1);
            k1 = k_list(idx_1, 1:2)';
            k2 = k_list(to_idx, 1:2)'; 
            
            % find coupled momentum states (NN coupling only)
            for j = 1:length(q_scatt)
                T_here = squeeze(Tmat(:, :, j));
                q_here = q_scatt(:, j); 
                q_diff = k1-k2;
                id = find(abs(q_diff(1, :)-q_here(1))<=tol & abs(q_diff(2, :)-q_here(2))<=tol);
%                 disp(q_diff(:,id))
%                 disp(isempty(id))
%                 disp(length(id))
                if ~isempty(id)
                    
                    AA_idx(1) = 2*(idx_1-1)+1;
                    AA_idx(2) = 2*(to_idx(id)-1)+1;
 
                    H_inter(AA_idx(1), AA_idx(2)) = T_here(1,1);
                    H_inter(AA_idx(1), AA_idx(2)+1) = T_here(1,2);
                    H_inter(AA_idx(1)+1, AA_idx(2)) = T_here(2,1);
                    H_inter(AA_idx(1)+1, AA_idx(2)+1) = T_here(2,2);
                end    
            end
            
                        
        end 
    end
     
    % making sure the Hamiltonian is Hermitian
    H_inter = H_inter + H_inter';    
    
end