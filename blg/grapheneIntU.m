function E = grapheneIntU(x,y,theta,s,a,o,m,orb_pos) % o:TO and m:FROM are orbitals
                                                        % s = 0 or 1, 1 is
                                                        % top

L = [1 cos(pi/3);0 sin(pi/3)]*a;

R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
I = eye(2);
T = 2*pi*(L')^(-1);

G = (I-R')*T;
%G = (R-I)*T;


zer = [0;0];
supercell = (I-R')^(-1);% * L;
%supercell = (R-I)^(-1);% * L;

% reciprocal lattice supercell = (I - R_(-theta)) * T

%[X,Y] = meshgrid(1:N);

%D = (1/N)*supercell*[X(:)';Y(:)'];
%D = D';

% going FROM layer 1, going TO layer 2
pos_from = orb_pos(1,m,:);
pos_to = orb_pos(2,o,:);
dpos = squeeze(pos_to - pos_from);

D = supercell*([x(:)';y(:)'] - dpos);
D = D';

z = 3.35*ones(1,size(D,1));
atom_pos_list = [D(:,1)';D(:,2)';z];

P = relax_geom_from_file(theta,(G*[1;0]),(G*[0;1]),zer',zer',s,atom_pos_list');

Pos = [x(:),y(:),z']+2*P;
%max(Pos(:,3))
%min(Pos(:,3))

% o = orbital TO, m = orbital FROM
E = dft_interlayer_coupling(Pos,theta,0,a,o,m);     % note that this only considers h^{12}
                                                % h^{21}(-r) = h^{12}(r)

%E = reshape(E,size(x,1),size(x,2));

%quiver(D(:,1),D(:,2),P(:,1),P(:,2))

%Relax_model_Stephen_ver3(rot_theta,sc_b1,sc_b2,sc_t1,sc_t2,top_layer,atom_pos_list)

end
