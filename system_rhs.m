function sys = system_rhs(t,w,dummy,A,B,C)
    nu = .001;
    
    
    psi = A \ w;
    psi_x = B*psi;
    psi_y = C* psi;
    w_x = B * w;
    w_y = C * w;
    sys = -psi_x .* w_y + psi_y .* w_x + nu* (A*w); 
end