classdef kDoF_tri < handle
    %kDoF creates and stores the k-space Degrees of Freedom for a trilayer
    %updated 11/22/2019, change how
    
    properties
        layer_list      % collection of Layer.m objects (input)
        max_k           % cutoff in momentum space (input)
        num_k           % total number of k-points (output)
        k_list          % list of k-points (output)
        grid_search;    % maximum size of the searching square (input)
    end
    
    methods
        function obj = kDoF_tri(Layers,max_k,grid_search)
            % set input variables (see properties for more info)
            obj.layer_list = Layers;
            obj.max_k = max_k;
            obj.grid_search = grid_search;
        end
        
        function [] = gen_dof(obj)
            % void function to create the k-point list
            
            clear obj.k_list;
            clear obj.num_k;
            grid_s = obj.grid_search;
            k_idx = 0;
            
            [~,ind] = min(abs([obj.layer_list(1).theta, obj.layer_list(3).theta]));
            kcut = norm(obj.layer_list(2*(ind-1)+1).G(:,1)-obj.layer_list(2).G(:,1));
            
            for i = 1:3 
                K(:,i) = 1/3 * (2*obj.layer_list(i).G(:,1) + obj.layer_list(i).G(:,2));
            end 
                     
            for l = 1:3 
                l_1 = l; % the layer we are in 
                l_2 = mod(l,3)+1;
                l_3 = mod(l+1,3)+1;
                
                kcut1 = obj.max_k*norm(K(:,2))*sqrt(3);
     
                b1 = obj.layer_list(l_1).G;
                b2 = obj.layer_list(l_2).G;
                b3 = obj.layer_list(l_3).G;
                
                % kcut2 = 4*norm(b1(:,2)); 
                kcut2 = kcut1; % constraining the norm of monolayer G vectors
                
                
                shift_vec=K(:,l_1)-K(:,2);
                grid_s1 = grid_s; 
                grid_s2 = grid_s;
%                 if l == 3 
%                     grid_s1 = 0;
%                 elseif l == 2
%                     grid_s2 = 0;
%                 else 
%                     grid_s1 = 0;
%                     grid_s2 = 0;
%                 end 
                for n2_1 = -grid_s1:grid_s1
                    for n2_2 = -grid_s1:grid_s1
                        k2 = b2(:,1)*n2_1 + b2(:,2)*n2_2; 
    
                        cond1 = norm(k2)<kcut2;
                        
                        if cond1
                            for n3_1 = -grid_s2:grid_s2
                                for n3_2 = -grid_s2:grid_s2

                                    k3 = b3(:,1)*n3_1 + b3(:,2)*n3_2;
                                    
                                    k1 = b1(:,1)*(n3_1+n2_1) + b1(:,2)*(n3_2+n2_2);
                                    k1p = (b2(:,1)*n3_1+b2(:,2)*n3_2)+(b3(:,1)*n2_1+b3(:,2)*n2_2);
                                    
                                    cond2 = norm(k2-k3)<kcut1;
                                    % cond2 = 1;
                                    
                                    cond3 = norm(k3)<kcut2;
                                    % cond3 = 1;
                                  
                                    if cond2 && cond3
                                        % mapping back to be near the origin of the given layer
                                        % and the origin is the K point of layer 2 (shifting by shift_vec) 
                                        k_here = k2+k3-k1+shift_vec;       
    
                                        k_in(1) = k_here(1);
                                        k_in(2) = k_here(2); 
                                        k_in(3) = l_1; 

                                        l1_index = 4 + (l_1 - 1)*2; % location of layer 1 in the data set
                                        l2_index = 4 + (l_2 - 1)*2; % location of layer 2 in the data set
                                        l3_index = 4 + (l_3 - 1)*2; % location of layer 3 in the data set  

                                        k_in(l1_index) = 0; 
                                        k_in(l1_index+1) = 0; 

                                        k_in(l2_index) = n2_1; % save grid index for this lattice
                                        k_in(l2_index+1) = n2_2;

                                        k_in(l3_index) = n3_1; % save grid index for this lattice
                                        k_in(l3_index+1) = n3_2;

                                        k_idx = k_idx+1;
                                        k_data(k_idx,:) = k_in; % push back the new DoF
                                    end 

                                end 
                            end 
                            
                        end 
                    end 
                end 
            end 
            
            
            obj.num_k = size(k_data, 1);
            
            % get rid of repeated indices
         
            
%             disp(tar_dofs)

            k_vec = k_data(:, 1:3);
            % disp(length(k_vec))
            
            idx = 1; 
            idx_list = ones([length(k_vec), 1]);
            tol = 1e-4;
            
            for a = 1:2
                for i = 1:size(k_vec,1)
                    k_here = k_vec(i,:);

                    if idx_list(i)
                        if a == 1
                            jbeg = 1;
                            jend = i-1;
                        else 
                            jbeg = i+1;
                            jend = size(k_vec,1);
                        end 
                            
                        for j = jbeg:jend
                            if abs(k_vec(j,1) - k_here(1)) < tol && ...
                                    abs(k_vec(j,2) - k_here(2)) < tol && abs(k_vec(j,3) - k_here(3)) < tol
                                if ~min( k_data(j,3:end) == [k_data(j,3) 0 0 0 0 0 0])
                                    idx_list(j) = 0;
                                    idx = idx + 1;
                                end 
                            end 
                        end 

                    end 
                end 
                
            end 
            
            idx_list = logical(idx_list);
            
            % disp(size(k_data))
            k_data = k_data(idx_list, :);
%             disp(size(k_data))
 
            obj.k_list = k_data;
        end
    end
end

