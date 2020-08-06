% author: stephen carr 
% generate the reciprocal lattice vectors 
% input: A: 
function G = getRecip(A)
% getRecip returns the recipriocal vector (identical to 2*pi*inv(A)' )

    R = [0 -1; 1 0];
    G(:,1) = 2*pi*(R*A(:,2))/(dot(A(:,1),R*A(:,2)));
    G(:,2) = 2*pi*(R*A(:,1))/(dot(A(:,2),R*A(:,1)));

end

