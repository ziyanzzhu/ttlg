% author: ziyan (zoe) zhu
% email: zzhu1@g.harvard.edu
% input: A0: lattice vecotrs
% theta: twist angle in rad. 
% delta: lattice mismatch
function Am = moireh_calc(A0, theta, delta, m, n) 
    
    rot = [cos(theta), -sin(theta);
           sin(theta), cos(theta)]; % ccw rotation
    A1 = A0; 
    A2 = (1+delta)*inv(rot)*A0; 

    g1 = 2*pi*transpose(inv(A1));
    g2 = 2*pi*transpose(inv(A2));

    G = m*g1 - n*g2;
    Am = transpose(inv(G)) * (2*pi);

    
    
