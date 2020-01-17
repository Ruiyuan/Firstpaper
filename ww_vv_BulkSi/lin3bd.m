%**********************************************************************
% Purpose:
%  Linearize three-body potential
%
% Arguments:
%  xi, xj, xk -- atomic coordinates
%  atypei, atypej, atypek -- atomic types
%  D(i1,xyz1,i2,xyz2) -- second derivative matrix
%    i1 -- atomic index
%    xyz1 -- coordinate index
%    i2 -- atomic index
%    xyz2 -- coordinate index
%
% Note:
%  Atom i should be the one in the center
%   SUBROUTINE lin3bd(xi, xj, xk, atypei, atypej, atypek, D)
function [D]=lin3bd(xi, xj, xk, atypei, atypej, atypek, boxlx, isigma, leps, gamma)

% TO TEST lin3bd function
%  xi=r(:,j), xj=r(:,i), xk=r(:,k), atypei=atype(j), atypej=atype(i), atypek=atype(k), boxlx= boxL,

%     REAL(MYREAL),DIMENSION(3)	:: xi, xj, xk
%     INTEGER			:: atypei, atypej, atypek
%     REAL(MYREAL),DIMENSION(3,3,3,3) :: D
% 
%     REAL(MYREAL),DIMENSION(3)	:: xij, xik
%     REAL(MYREAL)		:: rij, rik
    % costh -- cosine of angle jik
    % dcosthi -- d(costh)/d(xi)
    % dcosthj -- d(costh)/d(xj)
    % dcosthk -- d(costh)/d(xk)
    % dcosth(1,:) = dcosthi
    % dcosth(2,:) = dcosthj
    % dcosth(3,:) = dcosth3
%     REAL(MYREAL)		:: costh
%     REAL(MYREAL),DIMENSION(3)	:: dcosthi, dcosthj, dcosthk
%     REAL(MYREAL),DIMENSION(3,3)	:: dcosth
%     INTEGER		:: i, j, i0, j0
    % Initilizaton
    D = zeros(3,3,3,3);
%     global eps isigma leps epssig gamma; 
    xij = (xj - xi) - round((xj - xi)./boxlx).*boxlx;
%     xij = xij*isigma(atypei,atypej); % ?? Do I need use isigma or not
%     here ??
    rij = sqrt(dot(xij,xij));

    xik = (xk - xi) - round((xk - xi)./boxlx).*boxlx;
%     xik = xik*isigma(atypei,atypek);  % ?? Do I need use isigma or not
%     here ??
    rik = sqrt(dot(xik,xik));

    costh = dot(xij, xik)/(rij*rik);

%     dcosthj = xik/(rij*rik) - xij/(rij*rij)*costh;
%     dcosthk = xij/(rij*rik) - xik/(rik*rik)*costh;
%     dcosthj = xik/(rij*rik) - diag(diag(xij*xij')*xik')/(rij*rij)/(rij*rik);
    dcosthj = (xik - rik/rij*costh*xij)/(rij*rik);
    dcosthk = (xij - rij/rik*costh*xik)/(rij*rik); % xik/(rik*rik)*costh;
    dcosthi = -dcosthj - dcosthk; %%?? Is this correct? I don't think so. 

    dcosth(1,:) = dcosthi;
    dcosth(2,:) = dcosthj;
    dcosth(3,:) = dcosthk;
    

    
    for i = 1: 3;
    for j = 1: 3;
      for i0 = 1: 3;
      for j0 = 1: 3;
        D1(i,i0,j,j0) = dcosth(i,i0)*dcosth(j,j0);  % ?? which one is correct? I don't know. 
        D(i0,j0,i,j) = dcosth(i,i0)*dcosth(j,j0);
        D2(i,j,i0,j0) = dcosth(i,i0)*dcosth(j,j0);
      end % DO % j0
      end % DO % i0
    end % DO % j
    end % DO % i
    al = 1.8; 
    D1 = D1*2.0*leps(atypei,atypej,atypek)*(isigma(atypei,atypej)^2)*exp(gamma/(rij - al) + gamma/(rik -al));
    D = D*2.0*leps(atypei,atypej,atypek)*(isigma(atypei,atypej)^2)*exp(gamma/(rij - al) + gamma/(rik -al));
    

  end % SUBROUTINE lin3bd
