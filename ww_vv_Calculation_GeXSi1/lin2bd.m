%**********************************************************************
% Purpose:
%  Linearize two-body potential
% Arguments:
%  xi, xj -- atomic coordinates
%  atypei, atypej -- atomic types
%  D(i1,xyz1,i2,xyz2) -- second derivative matrix
%    i1 -- atomic index
%    xyz1 -- coordinate index
%    i2 -- atomic index
%    xyz2 -- coordinate index
function[D]=lin2bd(xi, xj, atypei, atypej, boxlx, eps, isigma)
%     REAL(MYREAL),DIMENSION(3)	:: xi, xj
%     INTEGER			:: atypei, atypej
%     REAL(MYREAL),DIMENSION(2,3,2,3) :: D

    % xij -- separation
    % rij -- xij/sigma
    % K -- spring constant of two-body interaction
    %      without epsilon and sigma
%     REAL(MYREAL),DIMENSION(3)	:: xij
%     REAL(MYREAL)		:: rij
%     REAL(MYREAL)		:: K
%     INTEGER			:: n, m
    D = zeros(3,3,2,2);
%     global eps isigma leps epssig gamma; 
    xij = (xi - xj) - round((xi - xj)./boxlx).*boxlx;
%     yij = (yi - yj) - round((yi - yj)/boxly)*boxly;
%     zij = (zi - zj) - round((zi - zj)/boxlz)*boxlz;
    
   % xij = xij*isigma(atypei,atypej); % ?? Do I need use isigma here or not? ??
    rij = sqrt(dot(xij, xij));

    for n = 1: 3;
    for m = 1: 3;
      D(n,m,1,1) = xij(n)*xij(m)/(rij^2);
    end % DO % m
    end % DO % n

    D(:,:,1,2) = -D(:,:,1,1); 
    D(:,:,2,1) = -D(:,:,1,1); 
    D(:,:,2,2) =  D(:,:,1,1); 

    D = D*eps(atypei,atypej)*(isigma(atypei,atypej)^2);
    K =20.8791000; % 0999900000; % 9; % 1; % 99999; % 8239242; % 48239242; %  20.879148239242; %   
    D = D*K;


end %SUBROUTINE lin2bd
