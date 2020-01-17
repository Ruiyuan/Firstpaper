function [Phiall, rxall, ryall, rzall]=dispersioncurveSuperCell_Modify(A,ncell,rc,Ksp,naxis,mm,nbasis,forceflag,fmatpure,xs,ys,zs,ncx,ncy,ncz,ncxsuper,nsupercell)
% kv = kvI; mm = mI ; fmatpure = fmatSi; 

x=xs;
y=ys;
z=zs;
rs=[xs;ys;zs];
% ncxsuper = 1; 
% ncx= 2; %  1; % 5;
% ncy= 2; % 5;
% ncz= 2; % 5;
% % % A1 = A*ncxsuper; % ncx;
boxlx=ncx*A;
boxly=ncy*A;
boxlz=ncz*A;

%    kx=kv(1);
%    ky=kv(2);
%    kz=kv(3);
 
 %----------positions of bulk material atoms. ---------------%
na=length(xs);
% % % ncell1 = ncell*nbasis*(ncxsuper^3); % ncx*ncy*ncz; % ncell*ncx*ncy*ncz; % ncell*nbasis ; % 
% % % nsupercell = ncell1; 
nprimcell = ncell;
nconvcell = ncell*nbasis;

%--over-----positions of bulk material atoms. -----over----%
 %%%%%%%%%%%%%%%%%%%Primitive Cell: N = 1*conventional cell %%%%%%%%%%%

        nm=0;
       for j=1:na/ncell;
           for i=1:ncell;
               nm=nm+1;
               primcell(i,j)=nm;
               if(rs(1:3,primcell(i,j))/A==[1;1;1])   % here is to find the center for furture calculation of fmatSi, fmatGe, fmatSiSiGe, and fmatSiGeGe. 
                   i0=i;
                   j0=j;
               end
           end
       end
 %%%%%%%%%%%%%%%%%%%SuperCell: N*conventional cell %%%%%%%%%%%
               nm=0;
       for j=1:na/nsupercell;
           for i=1:nsupercell;
               nm=nm+1;
               supercell(i,j)=nm;
%                if(rs(1:3,primcell(i,j))/A==[1;1;1])   % here is to find the center for furture calculation of fmatSi, fmatGe, fmatSiSiGe, and fmatSiGeGe. 
%                    i0=i;
%                    j0=j;
%                end
           end
       end
 %%%%%%%%%%%%%%%%%%%Conventional Cell: 1*conventional cell %%%%%%%%%%%
               nm=0;
       for j=1:na/nconvcell;
           for i=1:nconvcell;
               nm=nm+1;
               convcell(i,j)=nm;
%                if(rs(1:3,primcell(i,j))/A==[1;1;1])   % here is to find the center for furture calculation of fmatSi, fmatGe, fmatSiSiGe, and fmatSiGeGe. 
%                    i0=i;
%                    j0=j;
%                end
           end
       end
       
   D=zeros(naxis*ncell,naxis*ncell);
   Dsuper = zeros(naxis*nsupercell,naxis*nsupercell);
   Phiall = zeros(naxis*nsupercell,naxis*na);
   rxall = zeros(naxis*nsupercell,naxis*na);
   ryall = zeros(naxis*nsupercell,naxis*na);
   rzall = zeros(naxis*nsupercell,naxis*na);
   
 nbasis=4;  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%PRIMITIVE CELL CALCULATION FOR Dynamics
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%MATRIX.%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   for i=1:nsupercell; 
       
    for jsuper=1:nsupercell;
        DABIJ=zeros(naxis,naxis);
        for ksuper=1:na/nsupercell;
            i1super=supercell(i,1);
            j1super=supercell(jsuper,ksuper);
            
            i1 = i1super; %  primcell(i1primy,i1primx);
            j1 = j1super; % primcell(j1primy,j1primx);
            
            i1primx = round(i1super/nprimcell); % box number. 
            i1primy = mod(i1super+1,nprimcell)+1; % atom number. 
            j1primx = round(j1super/nprimcell); % box number; 
            j1primy = mod(j1super+1,nprimcell)+1; % atom number. 
            
% % %             i1convx = round(i1super/nconvcell); % box number. 
% % %             i1convy = mod(i1super+1,nconvcell)+1; % atom number. 
% % %             j1convx = round(j1super/nconvcell); % box number; 
% % %             j1convy = mod(j1super+1,nconvcell)+1; % atom number. 
            
% % %             i1voidcheck = mod(i1convx,36);
% % %             j1voidcheck = mod(j1convx,36);
% % %             voidconv = [8,11,26,29];
% % %             voidcheck = 0; 
% % %             for icheck = 1:length(voidconv)
% % %                 if or(i1voidcheck == voidconv(icheck), j1voidcheck == voidconv(icheck))
% % %                     voidcheck = 1; 
% % %                 end
% % %             end
% % %             if (voidcheck == 1)
% % %                 continue;
% % %             end
            
            ii=primcell(1,i1primx);
            jj=primcell(1,j1primx);
                    dxij=(x(jj)-x(ii)); %/(A/2);
                    dyij=(y(jj)-y(ii)); %/(A/2);
                    dzij=(z(jj)-z(ii)); %/(A/2);
                    dxij=dxij-boxlx*round(dxij/boxlx);
                    dyij=dyij-boxly*round(dyij/boxly);
                    dzij=dzij-boxlz*round(dzij/boxlz);
                    % end periodic boundary condition
                    dxij=round(dxij/(A/2));
                    dyij=round(dyij/(A/2));
                    dzij=round(dzij/(A/2));
                   
                   Phi=zeros(3,3); 
                   if(forceflag==1)
                      Phi=PhiABij(i,j1,xs,ys,zs,ncx,ncy,ncz,A,ncell,rc,nbasis,forceflag);
                   else
                     if(-1<=dxij & dxij<=1)
                       if(-1<=dyij & dyij<=1)
                           if(-1<=dzij & dzij<=1)
                               
%                                    Phi=fmatpure((i-1)*3+1:i*3, (jprime-1)*3+1:jprime*3, dxij+2, dyij+2, dzij+2);
                                   Phi=fmatpure((i1primy-1)*3+1:i1primy*3, (j1primy-1)*3+1:j1primy*3, dxij+2, dyij+2, dzij+2);
                               
                           end
                       end
                     end
                   end
           
            dx=xs(i1)-xs(j1);
            dy=ys(i1)-ys(j1);
            dz=zs(i1)-zs(j1);
            dx=dx-round(dx/boxlx)*boxlx;
            dy=dy-round(dy/boxly)*boxly;
            dz=dz-round(dz/boxlz)*boxlz;
            rirj=-[dx dy dz];
%             dabc=exp(1i*dot(kv,rirj));
%             dDABIJ=Phi*dabc*(1/sqrt(mm*mm));
%             DABIJ=DABIJ+dDABIJ        ;    
            
            % % % ### Construct Phiall for ncsuperll and all atoms.                 
                indxall = naxis*(i-1)+(1:naxis); %
                indyall = naxis*(j1super-1)+(1:naxis); %
                Phiall(indxall, indyall) = Phi(1:naxis, 1:naxis)*(1/sqrt(mm(i1)*mm(j1))); % 
                rxall(indxall, indyall) = -dx*ones(3,3);
                ryall(indxall, indyall) = -dy*ones(3,3);
                rzall(indxall, indyall) = -dz*ones(3,3);                
              
        end
        
% 
% %             for al=1:naxis;
% %                 for be=1:naxis;
      
% %                 indxsuper = naxis*(i-1)+1:naxis; % al;
% %                 indysuper = naxis*(jsuper-1)+1:naxis; % be;
% %                 Dsuper(indxsuper, indysuper) = DABIJ(1:naxis,1:naxis); 
% %                 end
% %             end
                       
    end
   end
% %  [UU,SS,VV]=svdecon(Dsuper);
% %  eVout=zeros(naxis*nsupercell,naxis*nsupercell);
% %  eVout=UU ; % 
% %  for i=1:naxis*nsupercell;
% %     om(i,1)=SS(i,i);
% %  end
% % om=(sqrt(om));   
%    om = sqrt(eig(Dsuper));
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%END CALCULATION OF DYNAMICAL MATRIX%%%%%%%

end
