function bl_rate = bl_rate_calc(C,R,m_idx,M,mu)
if (m_idx == 0)
    % bl_rate = C*2*R/3;
    bl_rate = C*mean(R);
else
    a_list = C*R;
    bl_rate_term = 0;
    for k = m_idx:M
        bl_rate_term_nr = prod(a_list(m_idx:k));
        bl_rate_term_dr = prod(a_list(1:k) + mu);
        bl_rate_term = bl_rate_term + bl_rate_term_nr/bl_rate_term_dr;
    end
    bl_rate = mu^m_idx*bl_rate_term;
end
end