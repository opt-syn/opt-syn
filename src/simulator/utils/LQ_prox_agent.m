function prox_out = LQ_prox_agent(c, v_self, v_other, Q_self, Q_other, b_self)
    n = length(v_self);
    vec_t = (1/c)*v_self - Q_other*v_other- b_self;
    prox_out = ((1/c)*eye(n) + Q_self) \ vec_t;
end