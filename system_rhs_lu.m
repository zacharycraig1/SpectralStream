function sys = system_rhs_lu(t,w,dummy,A,B,C,L,U,P)
    nu = .001;
    n = 64;
    y = L \ (P*w);
    psi = U \ y;
    psi_x = B*psi;
    psi_y = C* psi;
    w_x = B * w;
    w_y = C * w;
    sys = -psi_x .* w_y + psi_y .* w_x + nu* (A*w); 
end
