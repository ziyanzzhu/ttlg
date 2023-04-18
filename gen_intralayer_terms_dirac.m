% authors: ziyan (zoe) zhu, stephen carr
% email: zzhu1@g.harvard.edu
% generating intralayer Hamiltonian (rotated Dirac Hamiltonian)
% input: k_list: list of k in the basis
%        layers: `layers` objected generated using `Layer.m`
%        tar_q: center site
%        E_field: non-zero value if E_field is turned on
function H = gen_intralayer_terms_dirac(k_list,layers,tar_q,E_field)
    
    vf = 6.582*0.8; % Fermi velocity, 6.582 is the conversion factor from hbar
    
    ndof = size(k_list,1); 
    H = sparse(2*ndof,2*ndof); 
    
    layer_list = k_list(:, 3);
    for j = 1:3
        G = layers(j).G;
        K(:,j) = 1/3 * (2*G(:,1) + G(:,2));
    end 
    
    ind = 0;
    
    if E_field 
        V = [-E_field/2, 0, E_field/2]; 
    else 
        V = [0, 0, 0];
    end
    
    for i = 1:3
        l_1 = i;
        l_2 = mod(l_1,3)+1;
        l_3 = mod(l_1+1,3)+1;

        k_here = k_list(layer_list==l_1, :);
        
        for n = 1:size(k_here, 1)
            theta_here = layers(l_1).theta;    
            
            k_val = k_here(n,1:2) + tar_q;    
            kx = k_val(1);
            ky = k_val(2);
            
            ind = ind+1;
            ibeg = 2*(ind-1)+1;
            
            % rotation
            prefac_on = 1; 
            prefac = exp(-1j*prefac_on*theta_here);
            
            H(ibeg:ibeg+1, ibeg:ibeg+1) = [V(i), vf*(kx+1j*ky)*prefac; 
                                             vf*(kx-1j*ky)*conj(prefac), V(i)];
                 
        end 
        
        
    end

end

