function H = gen_intralayer_terms_dirac(k_list,layers,tar_q,E_field)
    
    vf = 6.582*0.8;
    
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
        l1_index = 4 + (l_1 - 1)*2;
        l2_index = 4 + (l_2 - 1)*2;
        l3_index = 4 + (l_3 - 1)*2;
          
        q1=K(:,l_1)-K(:,l_2);
        q1=K(:,l_1)-K(:,2);
        rot120 = [cos(2*pi/3), sin(2*pi/3); -sin(2*pi/3), cos(2*pi/3)];
        q2=rot120*q1;
        q3=rot120*q2;

        k_here = k_list(layer_list==l_1, :);
        
        for n = 1:size(k_here, 1)
            theta_here = layers(l_1).theta;    
            
            k_val = k_here(n,1:2) + tar_q; 
                   
            kx = k_val(1);
            ky = k_val(2);
            
            ind = ind+1;
            k_all(ind,:) = k_val;
            ibeg = 2*(ind-1)+1;
            
            prefac_on = 1; 
            prefac = exp(-1j*prefac_on*theta_here);
            
            H(ibeg:ibeg+1, ibeg:ibeg+1) = -vf*[V(i), (kx+1j*ky)*prefac; 
                                               (kx-1j*ky)*conj(prefac), V(i)];
                 
        end 
        
        
    end

end

