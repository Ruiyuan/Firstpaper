function [inc_scat_inv] = incoherent_scat(ncx,ncxdot,dA,dK,kx,ky,kz,a0,vmode,sigmaSi);

    k = sqrt(kx^2 + ky^2 + kz^2)/sigmaSi;
    R = 0.6204*ncxdot*a0; % (3/4/pi)^(1/3)*ncxdot*a0;
    chi = k*R;
    G = pi*R^2;
    alpha = 1/sqrt(2);
    delta = a0/4*sqrt(3); % delta is nearest neighbor distance. 
    h1 = (alpha^2/4*dA^2 + 3*alpha^8*dK^2*(sin(alpha*k*delta/2))^4/((alpha*k*delta/2)^4))*(pi*(cos(4*chi)-1+4*chi*sin(4*chi)+32*chi^4-8*chi^2)/16/chi^6); 
    sig_far = G*chi^4*h1;
    dq = sqrt(1+dA)/sqrt(1+dK) - 1;
    h2 = (1 - sin(2*chi*dq)/(chi*dq) + (sin(chi*dq))^2/(chi*dq)^2);
    sig_near = 2*G*h2; 
    sig_inc = 1/(1/sig_far + 1/sig_near);
    eta = 1/(ncx*a0)^3;
    inc_scat_inv = eta*sig_inc*abs(vmode);
end
