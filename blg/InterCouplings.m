% inherit from "handle" to get OOP
classdef InterCouplings < handle
    % Basic datastructure to store geometric info for a 2D layer.
    
    properties
        layers
        layer_from        % index in the heterostructure
        layer_to
        tar_k
        
        R_meshx % follow moire pattern? inv(G1 - G2)?
        R_meshy
        
        B_meshx
        B_meshy
        B_meshz
        
        K_meshx % basis should follow moire recip. vectors, G1 - G2
        K_meshy
        
        tR
        tK
        
        theta % twist angle, in degrees
        orb_pos % positions of orbitals
        dA % area of volume element
        
        u_relax
    end
    
    methods
        function obj = InterCouplings(ls,l_from,l_to, r_max, r_l, k_max, k_l, k_inner, tar_k, theta_in, pos_in)
            %Layer Construct an instance of this class
            %   Detailed explanation goes here
            obj.layers = ls;
            obj.layer_from = l_from;
            obj.layer_to = l_to;
            obj.theta = theta_in;
            obj.orb_pos = pos_in;
            
            G1 = ls(l_from).G;
            G2 = ls(l_to).G;
            G0 = (G1 + G2)/2;
            obj.tar_k = tar_k;
            
            tot_theta = ls(2).theta - ls(1).theta;
            
            b12 = (G2 - G1); % reciprocal vectors.
            r12 = 2*pi*inv(b12'); % moire supercell
            
            nr = floor(norm(r12(:,1))/r_l);
            max_grid = floor(1.5*r_max/r_l);
            Rsc_base = r12/nr;
            R_unitarea = cross([Rsc_base(:,1);0],[Rsc_base(:,2);0]);
            obj.dA = R_unitarea(3);
            
            [mesh_x, mesh_y] = meshgrid(-max_grid:max_grid, -max_grid:max_grid);
            
            R_bigmeshx = mesh_x*Rsc_base(1,1) + mesh_y*Rsc_base(1,2); 
            R_bigmeshy = mesh_x*Rsc_base(2,1) + mesh_y*Rsc_base(2,2);
           
            R_keep_idx = (R_bigmeshx(:).^2 + R_bigmeshy(:).^2 < r_max^2);
            
            obj.R_meshx = R_bigmeshx(R_keep_idx);
            obj.R_meshy = R_bigmeshy(R_keep_idx);
            
            norbs = size(obj.orb_pos,2);
            for o1 = 1:norbs
                for o2 = 1:norbs
                    obj.B_meshx(o1,o2,:) = obj.R_meshx;
                    obj.B_meshy(o1,o2,:) = obj.R_meshy;
                    obj.B_meshz(o1,o2,:) = 3.35 + zeros(size(obj.R_meshx));
                    obj.u_relax(o1,o2,:,:) = zeros(size(obj.R_meshx));
                end
            end
            
            nk = floor(norm(norm(b12(:,1))/k_l));
            max_gridk = floor(1.5*k_inner/k_l);
            Gsc_base = b12/nk;
            
            max_gridk = max_gridk;
            
            [mesh_kx, mesh_ky] = meshgrid(-max_gridk:max_gridk, -max_gridk:max_gridk);
            
            K_bigmeshx = mesh_kx*Gsc_base(1,1) + mesh_ky*Gsc_base(1,2); 
            K_bigmeshy = mesh_kx*Gsc_base(2,1) + mesh_ky*Gsc_base(2,2);
           
            K_keep_idx = (K_bigmeshx(:).^2 + K_bigmeshy(:).^2 < k_inner^2);
            
            K_localmeshx = K_bigmeshx(K_keep_idx);
            K_localmeshy = K_bigmeshy(K_keep_idx);
            
            max_gridG = 2*floor(2.5*k_max/norm(G1(:,1)));
            [mesh_Gx, mesh_Gy] = meshgrid(-max_gridG:max_gridG, -max_gridG:max_gridG);

            dKx = mesh_Gx*G0(1,1) + mesh_Gy*G0(1,2);
            dKy = mesh_Gx*G0(2,1) + mesh_Gy*G0(2,2);
            
            K_finalmeshx = dKx(:)' + K_localmeshx + tar_k(1); 
            K_finalmeshy = dKy(:)' + K_localmeshy + tar_k(2); 

            %{
            for gidx1 = -max_gridG:max_gridG
                for gidx2 = -max_gridG:max_gridG
                    dKx = gidx1*G0(1,1) + gidx2*G0(1,2);
                    dKy = gidx1*G0(2,1) + gidx2*G0(2,2);
                    K_localmeshx + dKx;
                    K_localmeshy + dKy;
                end
            end
            %}
           
            K_final_keep_idx = (K_finalmeshx(:).^2 + K_finalmeshy(:).^2 < k_max^2);

            obj.K_meshx = K_finalmeshx(K_final_keep_idx);
            obj.K_meshy = K_finalmeshy(K_final_keep_idx);
   
        end
        
        function relax_configs(obj, relax_str)

            x = obj.R_meshx;
            y = obj.R_meshy;
            %z = obj.B_meshz;
            theta_rad = obj.theta*pi/180;
            s = false; %true: top layer, false: bottom layer
            
            L = obj.layers(1).A;
    
            R = [cos(theta_rad) -sin(theta_rad); sin(theta_rad) cos(theta_rad)];
            I = eye(2);
            T = 2*pi*inv(L');

            G = (I-R')*T;

            supercell = inv(I-R);
            
            norbs = size(obj.orb_pos,2);
            for o1 = 1:norbs
                for o2 = 1:norbs

                    % going FROM layer 1, going TO layer 2
                    pos_from = obj.orb_pos(1,o2,:);
                    pos_to = obj.orb_pos(2,o1,:);
                    dpos = squeeze(pos_to - pos_from);

                    D = supercell*([x(:)';y(:)'] - dpos);
                    D = D';

                    z = 3.35*ones(1,size(D,1));
                    atom_pos_list = [D(:,1)';D(:,2)';z];

                    u_relax_h = relax_str*relax_geom_from_file(theta_rad,(G*[1;0]),(G*[0;1]),s,atom_pos_list');

                    Bx(o1,o2,:) = x + 2*u_relax_h(:,1);
                    By(o1,o2,:) = y + 2*u_relax_h(:,2);
                    Bz(o1,o2,:) = z' - 2*u_relax_h(:,3);
                    u_rlx(o1,o2,:,:) = u_relax_h;
                end
            end
            
            obj.B_meshx = Bx;
            obj.B_meshy = By;
            obj.B_meshz = Bz;
            obj.u_relax = u_rlx;
            
        end
 
        function gen_tR(obj)
            
            norbs = size(obj.orb_pos,2);
            a = norm(obj.layers(1).A(:,1));
            x = obj.B_meshx;
            y = obj.B_meshy;
            z = obj.B_meshz;
            theta_h = obj.theta*pi/180;
            for o1 = 1:norbs
                for o2 = 1:norbs
                    xh = squeeze(x(o1,o2,:));
                    yh = squeeze(y(o1,o2,:));
                    zh = squeeze(z(o1,o2,:));
                    Pos = [xh, yh, zh];
                    tR_mat(o1,o2,:) = dft_interlayer_coupling(Pos,theta_h,0,a,o1,o2); 
                end
            end
            obj.tR = tR_mat;
            
        end
        
        function gen_tK(obj)
            
            A1 = obj.layers(1).A;
            a_uc = cross([A1(:,1);0],[A1(:,2);0]);
            a_uc = a_uc(3);
            dArea = obj.dA;

            Rx = obj.R_meshx;
            Ry = obj.R_meshy;
            Kx = obj.K_meshx;
            Ky = obj.K_meshy;  

            norbs = size(obj.orb_pos,2);


            for o1 = 1:norbs
                for o2 = 1:norbs
                    tR_h = squeeze(obj.tR(o1,o2,:));
                    phase_fac = exp(-1j*(Rx*Kx' + Ry*Ky'));
                    tK_h= dArea*transpose(sum(phase_fac.*tR_h,1))/a_uc;
                    tK_mat(o1,o2,:) = tK_h;
                end
            end

            obj.tK = tK_mat;

            
        end
        
    end
end

