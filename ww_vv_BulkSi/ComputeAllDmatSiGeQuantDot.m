 % ! Computate dynamical matrix for all "special" atoms (all atoms) 
 function [dmatallD,dmatallD2,fmat,D2,D3] = ComputeAllDmatSiGeQuantDot(xs,ys,zs,ncx,ncy,ncz,A,atype, eps, isigma,leps, gamma)

 % Global variables 
 
%     lambdaSi = 21.0; lambdaGe = 31.0; 
%     epsSi = 3.473928e-19; epsGe = 3.085e-19; 
%     % epsSi = 3.4723e-19; epsGe = 3.085e-19;
%     
%     sigmaSi = 2.0951e-10; sigmaGe = 2.0951e-10; 
%    % massSi=28.085*1.6605e-27; massGe=72.64*1.6605e-27;
%     massSi = 4.6637e-26; massGe = 1.2057e-25; 
%     [eps,isigma,leps,epssig, imass] = SiGein(lambdaSi, lambdaGe, epsSi, epsGe, sigmaSi, sigmaGe, massSi, massGe);
%     gamma=1.20d0;
%     global eps isigma leps epssig gamma; 
 
 natm = length(xs);
 r = [xs; ys; zs];
 boxL=[ncx;ncy;ncz]*A;
 
 nround2=-12;
 nround3= -12; 
 
 % Initialization 
 dmatallD=zeros(3,3,2,natm);
 dmatallD2=zeros(3,3,natm,natm); % % % 
 % ModNlist
 [nbr] = InitNlist(r, boxL, A);
 
%  for i = 1:1; % natm; % 
nic=0;
nib=0;
for i = 1:2;
    %-----------------------------
     % 2bd interaction
     for n1 = 1:4; 
         j = nbr(i,n1);
         D2 = lin2bd(r(:,i), r(:,j), atype(i), atype(j), boxL, eps, isigma);
%          D2 = roundn(D2, nround2);
         list2(1) = i;
         list2(2) = j;
         dmatallD = Add2bdDmat(dmatallD, list2, D2, nbr);
     end % n1
%-------------------------------         
     % 3bd interaction: i is the center point. 
     for n1 =1: 3; % 1:4; % 
         for n2 =n1+1:4; % 1:4; % 
             j = nbr(i,n1);
             k = nbr(i,n2);
             if(j==i || k==i || k==j)
                 continue;
             end
             
             D3 = lin3bd(r(:,i), r(:,j), r(:,k), atype(i), atype(j), atype(k), boxL, isigma, leps, gamma);
%              D3 = roundn(D3, nround3);
             nic=nic+1;
             list3(1) = i;
             list3(2) = j;
             list3(3) = k;
             dmatallD = Add3bdDmati(dmatallD, list3, D3, nbr, i);
         end % n2
     end % n1
     
     %--------------------------------------------------------
     % 3bd interaction: i is not the center point. 
     % 3bd interaction
     for n1 = 1: 4;
         j = nbr(i,n1);
         for n2 = 1:4;
             
             k = nbr(j,n2);
             if(k==i)        % ! This means just two atoms. wrong. 
                 continue;
             end
             
%              for n3 = 1:4;
                 if(min(abs(nbr(i,1:4) - k)) == 0);     % k is not i's nearest neighbor. we already calcuated it in above. 
                     continue;
                 end
                   
             
             D3 = lin3bd(r(:,j), r(:,i), r(:,k), atype(j), atype(i), atype(k), boxL, isigma, leps, gamma);
%              D3 = roundn(D3, nround3);
             nib=nib+1;
             list3(1) = j;
             list3(2) = i;
             list3(3) = k;
             dmatallD = Add3bdDmatj(dmatallD, list3, D3, nbr, i);
         end % n2
     end % n1
%     
end

  nm=0;
  na=natm;
  ncell=2;
  fmat=zeros(6,6,3,3,3);
  x=xs; y=ys; z=zs;
  boxlx=boxL(1); boxly=boxL(2); boxlz=boxL(3);

       
%--------------------------------------
% (1,1,1) as center to calculat dmatallD. 
%--------------------------------------
       for i = 1: natm; % %  % Calculate All atoms in SuperCell % % %  primcell(i0,j0):primcell(i0,j0)+1;
    %-----------------------------
     % 2bd interaction
     for n1 = 1:4; 
         j = nbr(i,n1);
         D2 = lin2bd(r(:,i), r(:,j), atype(i), atype(j), boxL, eps, isigma);
%          D2 = roundn(D2, nround2);
         list2(1) = i;
         list2(2) = j;
         dmatallD2 = Add2bdDmat(dmatallD2, list2, D2, nbr);
     end % n1
%-------------------------------         
     % 3bd interaction: i is the center point. 
     for n1 =1: 3; % 1:4; % 
         for n2 =n1+1:4; % 1:4; % 
             j = nbr(i,n1);
             k = nbr(i,n2);
             if(j==i || k==i || k==j)
                 continue;
             end
             
             D3 = lin3bd(r(:,i), r(:,j), r(:,k), atype(i), atype(j), atype(k), boxL, isigma, leps, gamma);
%              D3 = roundn(D3, nround3);
             nic=nic+1;
             list3(1) = i;
             list3(2) = j;
             list3(3) = k;
             dmatallD2 = Add3bdDmati(dmatallD2, list3, D3, nbr, i);
         end % n2
     end % n1
     
     %--------------------------------------------------------
     % 3bd interaction: i is not the center point. 
     % 3bd interaction
     for n1 = 1: 4;
         j = nbr(i,n1);
         for n2 = 1:4;
             
             k = nbr(j,n2);
             if(k==i)        % ! This means just two atoms. wrong. 
                 continue;
             end
             
%              for n3 = 1:4;
                 if(min(abs(nbr(i,1:4) - k)) == 0);     % k is not i's nearest neighbor. we already calcuated it in above. 
                     continue;
                 end
                   
             
             D3 = lin3bd(r(:,j), r(:,i), r(:,k), atype(j), atype(i), atype(k), boxL, isigma, leps, gamma);
%              D3 = roundn(D3, nround3);
             nib=nib+1;
             list3(1) = j;
             list3(2) = i;
             list3(3) = k;
             dmatallD2 = Add3bdDmatj(dmatallD2, list3, D3, nbr, i);
         end % n2
     end % n1
%     
       end
%-END (1,1,1) as center. 
%-----------------------------------------
      
      
 end
