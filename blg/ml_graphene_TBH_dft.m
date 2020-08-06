function H = ml_graphene_TBH_dft(x,y,orb_pos)
    delta = .01;

    a = 1.42*sqrt(3);
    %L = [1 cos(pi/3);0 sin(pi/3)]*a;
    %orb_pos = zeros(2,2);
    %orb_pos(1,:) = (L*[1/3; 1/3])';
    %orb_pos(2,:) = (L*[2/3; 2/3])';

    H = zeros(size(x,1),2,2);

    t1 = -2.892;
    t2 =  0.243;
    t3 = -0.266;
    t4 =  0.024;

    r_list = [a/sqrt(3), a, a*sqrt(3)/1.5, a*sqrt(7/3),sqrt(3)*a];
    r_sq = r_list.^2;

    for t = 1:2
        for f = 1:2
            v = [x,y] + ones(size(x,1) ,1)*(orb_pos(t,:) - orb_pos(f,:)); % here are the positions for given orbital set
            rs = v(:,1).*v(:,1) + v(:,2).*v(:,2); % radius squared per site distance 

            H(:,t,f) =              t1 * (rs > r_sq(1) - delta & rs < r_sq(1) + delta) ...
                                +   t2 * (rs > r_sq(2) - delta & rs < r_sq(2) + delta) ...
                                +   t3 * (rs > r_sq(3) - delta & rs < r_sq(3) + delta) ...
                                +   t4 * (rs > r_sq(4) - delta & rs < r_sq(4) + delta);
        end
    end

end

