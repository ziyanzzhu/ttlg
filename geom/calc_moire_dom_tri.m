function Am_dom = calc_moire_dom_tri(q12, q23, A0)
    grid_search = 8; 

    A12 = moireh_calc(A0, q12, 0, 1, 1);
    A23 = moireh_calc(A0, q23, 0, 1, 1);
    G12 = 2*pi*transpose(inv(A12));
    G23 = 2*pi*transpose(inv(A23));

    a12 = -A12(:,1);
    a23 = A23(:,1);

    g12 = -G12(:,1);
    g23 = G23(:,1);

    % bilayer moire length
    ml12 = norm(a12);
    ml23 = norm(a23);

    % angle between the bilayer moires 
    theta = acos(dot(a12, a23)/(ml12*ml23));

    % find the dominant harmonic 
    idx = 1;
    for i1 = 1:grid_search 
        for i2 = 1:grid_search 
            g_vec = i1*g12 - i2*g23; 
            g_vec_norm(i1,i2) = norm(g_vec);
        end 
    end 

    gvecnorm_list = sort(g_vec_norm);

    [m_dom,n_dom] = find(g_vec_norm==min(g_vec_norm(:)));
    m_dom = m_dom;
    n_dom = n_dom;
%     disp('---------------')
%     disp([m_dom, n_dom]);

    delta = ml12/ml23*m_dom/n_dom-1;
    Am_dom = moireh_calc(-A12, theta, ml23/ml12-1, m_dom, n_dom);
end 
