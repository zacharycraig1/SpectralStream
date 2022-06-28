function sys = sys_rhs_bic(t,w,dummy,A,B,C)
    tol = 10^-1;
    nu=.001;
    [psi,flag,relres,iter] = bicgstab(A,w,tol);
    
    psi_x = B*psi;
    psi_y = C* psi;
    w_x = B * w;
    w_y = C * w;
    sys = -psi_x .* w_y + psi_y .* w_x + nu* (A*w); 
end