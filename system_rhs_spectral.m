function sys = system_rhs_spectral(t,w,dummy,A,B,C)
    nu = .001;
    L = 20;
    n = 64;
    kx = (2*pi / L) * [0:(n/2-1) (-n / 2):-1];
    ky = (2*pi / L) * [0:(n/2-1) (-n / 2):-1];
    kx(1) = 1e-06;
    ky(1) = 1e-06;
    [Kx,Ky] = meshgrid(kx,ky);
    fomega = fft2(reshape(w,[n,n]));
    
    fpsi = -fomega./(Kx.^2 + Ky.^2);
    inversefpsi = ifft2(fpsi);
    
    psi = reshape(ifft2(fpsi),[1,4096]).';
    
    psi_x = B*psi;
    psi_y = C*psi;
    w_x = B*w;
    w_y = C*w;
    sys = -(psi_x .* w_y) + (psi_y .* w_x) + nu * (A * w);
    
    
    
end