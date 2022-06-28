function sys = sys_rhs_gmres(t,w,dummy,A,B,C)
    nu = .001;
    
    tol = 10^-6;
    
    [psi,flag,relres,iter] = gmres(A,w,10,tol,500);
    
    
    psi_x = B*psi;
    psi_y = C* psi;
    w_x = B * w;
    w_y = C * w;
    sys = -psi_x .* w_y + psi_y .* w_x + nu* (A*w); 
end