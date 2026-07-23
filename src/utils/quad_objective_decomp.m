function [quad] = quad_objective_decomp(M_quad, ind_p, ind_q)
    %QUAD_OBJECTIVE untangle the quadratic objective into a
    %linearizable formulation

    %R = T' U^-1 T, R >0
   
    %use eigenvalue arguments here

    Qq = M_quad(ind_q, ind_q);
    Sq = M_quad(ind_q, ind_p);
    Rq = M_quad(ind_p, ind_p);


    [RqV, RqD] = eig(Rq);
    eRq = diag(RqD);
    ind_pos = find(abs(eRq) > 1e-12);

    Tq = RqV(:, ind_pos)';
    Uq = diag(1./eRq(ind_pos));

    quad = quad_param(Qq, Sq, Uq, Tq);
end