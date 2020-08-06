% authors: stephen carr
% generating intralayer Hamiltonian (rotated Dirac Hamiltonian)
% input: k_list: list of k in the basis
%        layers: `layers` objected generated using `Layer.m`
%        tar_q: center site
%        E_field: non-zero value if E_field is turned on
function H = gen_intralayer_terms_dft(k_list,layers,q_here,E_field)
    
    dirac_shift = 0.729; % to make Dirac cones at 0 eV in monolayer
    
    ndof = size(k_list,1); 
    H = zeros(2*ndof,2*ndof); 
    
    layer_list = k_list(:, 5);
    
    % keeps track of the index of the Hamiltonian when unzipping
    ind = 1;
    
    if E_field 
        V = [-E_field/2, E_field/2]; 
    else 
        V = [0, 0];
    end
    
    for i = 1:2
        
        l_1 = i;
        %l_2 = mod(l_1,3)+1;
        A = layers(i).A;
        orb_pos = layers(i).orbPos;

        Q = k_list(layer_list==l_1, :);
        Q = Q(:,1:2) + q_here;
        
        % use vectorized form to quickly compute all H terms
        H_zip = zeros(size(Q,1),2,2);
        for j = -3:3
            for k = -3:3
                R = A*[j;k];
                B = squeeze(ml_graphene_TBH_dft(R(1),R(2),orb_pos));
                for f = 1:2 % from orb A/B
                    for t = 1:2 % to orb A/B
                        phase = exp(-1i*Q*R);
                        H_zip(:,t,f) = H_zip(:,t,f) + B(t,f)*phase;
                    end
                end

            end
        end
        
        % onsite additions
        H_zip(:,1,1) = H_zip(:,1,1) + dirac_shift + V(i);
        H_zip(:,2,2) = H_zip(:,2,2) + dirac_shift + V(i);
        
        % unzip the vectorized Hamiltonian
        for n = 1:size(Q, 1)
             
            ibeg = 2*(ind-1)+1;
            H(ibeg:ibeg+1, ibeg:ibeg+1) = H_zip(n,:,:);
            ind = ind+1;
                 
        end 
        
        %{
        for n = 1:size(k_here, 1)
            
            k_val = k_here(n,1:2) + q_here;    
            kx = k_val(1);
            ky = k_val(2);
            
            ind = ind+1;
            ibeg = 2*(ind-1)+1;
            
            % from orbital 2, to orbital 1: so (1,2) element of H_ml
            term_h = t1*(1 + exp(-1j*dot(A(:,1),k_val)) + exp(-1j*dot(A(:,2),k_val)));
            
            H(ibeg:ibeg+1, ibeg:ibeg+1) = [V(i), term_h; 
                                      conj(term_h), V(i)];
                 
        end 
        %}
        
    end

end

