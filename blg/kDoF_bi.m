classdef kDoF_bi < handle
    %kDoF creates and stores the k-space Degrees of Freedom for a bilayer
    
    properties
        layer_list      % collection of Layer.m objects (input)
        max_k           % cutoff in momentum space (input)
        num_k           % total number of k-points (output)
        k_list          % list of k-points (output)
        grid_search;    % maximum size of the searching square (input)
    end
    
    methods
        function obj = kDoF_bi(Layers,max_k,grid_search)
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
            
            for i = 1:2 
                K(:,i) = 1/3 * (2*obj.layer_list(i).G(:,1) + obj.layer_list(i).G(:,2));
            end 
                     
            for l = 1:2
                l_1 = l; % the layer we are in 
                l_2 = mod(l,2)+1;
                
                kcut = obj.max_k*norm(K(:,l_1))*sqrt(3);
     
                b1 = obj.layer_list(l_1).G;
                b2 = obj.layer_list(l_2).G;
                
                shift_vec=K(:,1)-K(:,l_1);
                grid_s1 = grid_s; 
                for n2_1 = -grid_s1:grid_s1
                    for n2_2 = -grid_s1:grid_s1
                        k2 = b2(:,1)*n2_1 + b2(:,2)*n2_2; 
                        
                        if  norm(k2) < kcut 
                            k1 =  b1(:,1)*n2_1 + b1(:,2)*n2_2; 
                            % mapping back to be near the origin of the given layer
                            % and the origin is the K point of layer 2 (shifting by shift_vec) 
                            k_here = k2-k1+shift_vec;       

                            k_in(1) = k2(1);
                            k_in(2) = k2(2);
                            k_in(3) = k_here(1);
                            k_in(4) = k_here(2);
                            k_in(5) = l_1; 

                            l1_index = 6 + (l_1 - 1)*2; % location of layer 1 in the data set
                            l2_index = 6 + (l_2 - 1)*2; % location of layer 2 in the data set

                            k_in(l1_index) = 0; 
                            k_in(l1_index+1) = 0; 

                            k_in(l2_index) = n2_1; % save grid index for this lattice
                            k_in(l2_index+1) = n2_2;

                            k_idx = k_idx+1;
                            k_data(k_idx,:) = k_in; % push back the new DoF                            
                        end 
                    end 
                end 
            end 
            
            
            obj.num_k = size(k_data, 1);
            obj.k_list = k_data;
        
        end
    end
end

