clear;
clc;
clear  lambdaSi lambdaGe epsSi epsGe sigmaSi sigmaGe massSi massGe eps isigma leps epssig gamma;




% Atom parameters; %
% lambdaSi = 21.0; lambdaGe = 31.0; epsSi = 3.473928e-19; epsGe = 3.085e-19; sigmaSi = 2.0951e-10; sigmaGe = 2.0951e-10; 
%     massSi=28.085*1.6605e-27; massGe=72.64*1.6605e-27;
%     mass=[massSi,massGe];
%----------------------------------------------
lambdaSi = 21.0; lambdaGe = 31.0; 
     epsSi = 3.473928e-19; epsGe = 3.085e-19; 
    % epsSi = 3.47392e-19; epsGe = 3.085e-19; 
    % epsSi = 3.4723e-19; epsGe = 3.085e-19;
    % epsSi = 3.473768059e-19; epsGe = 3.085e-19;  % make up values to make sure 3bd value is same to Hong Zhao's 3bd value. 
    
    sigmaSi = 2.0951e-10; sigmaGe = 2.0951e-10; 
   % massSi=28.085*1.6605e-27; massGe=72.64*1.6605e-27;
    massSi = 4.6637e-26; massGe = 1.2057e-25; 
    [eps,isigma,leps,epssig, imass] = SiGein(lambdaSi, lambdaGe, epsSi, epsGe, sigmaSi, sigmaGe, massSi, massGe);
    gamma=1.20d0;
    a0=5.43095e-10;
   % global lambdaSi lambdaGe epsSi epsGe sigmaSi sigmaGe massSi massGe eps isigma leps epssig gamma; 
%END Atom parameters. 

A= a0/sigmaSi; %  15168727030; 2.592215; % 2.592215168; %
ncell=2;
nbasis=4;
naxis=3;
forceflag=4;  % forceflag=1 : spring; forceflag=2: Silicon; forceflag=3: Germanium.  means choose s-w potential.
% ncx=2;
% ncy=2;
% ncz=2;
% na=ncx*ncy*ncz*ncell*nbasis;
% boxlx=(ncx)*A;
% boxly=(ncy)*A;
% boxlz=(ncz)*A;
rc=1.8;

KA=1;
KB=1;
KC=1;
mI=1;
% mT=0.8;
nT=0;
mTmatrix=[massGe/massSi]; % because of PhiGe has included the mass of massGe, so I just need reduced here as mT=1. 
colormatrix=['b','y','g'];

%---------start--contstruct "force matrix" of pure Si or heavy Si-------
ncx=2;
ncy=2;
ncz=2;

boxlx=ncx*A;
boxly=ncy*A;
boxlz=ncz*A;
na=ncx*ncy*ncz*ncell*nbasis;
[xs ys zs]=creatdiamond(A,ncx,ncy,ncz,ncell,nbasis);   
% Phipure=zeros(naxis,naxis,ncell,na);
% PhiSi=zeros(naxis,naxis,ncell,na);
% PhiGe=zeros(naxis,naxis,ncell,na);
% atypepure=zeros(2,na);
% atypepure(1,:)=ones(1,na);
% atypepure(2,:)=2*ones(1,na);  % here "2" means Ge. 
% for i=1:ncell; 
%        
%     for jprime=1:ncell;
%         DABIJ=zeros(naxis,naxis);
%         for kprime=1:na/ncell;
%             j1=(kprime-1)*ncell+jprime;
%             
%            for ii=1:length(atypepure(:,1))
%              Phi=PhiABij(i,j1,xs,ys,zs,ncx,ncy,ncz,A,ncell,rc,nbasis,forceflag,atypepure(ii,:));
% %            Phipure(:,:,i,j1)=Phi;
%                 if ii==1;
%                     PhiSi(:,:,i,j1)=Phi;
%                 else
%                     PhiGe(:,:,i,j1)=Phi;
%                 end
%            end
% 
%         end
% 
%     end
% end

 %continue to constructure fmatpure[-1:1;-1:1;-1:1;6,6] matrix for
       %convinence
      x=xs; y=ys; z=zs;
      rs=[x;y;z];
      % to for different atom types at different interfaces. 
      atype=zeros(4,na);
      ainterface=[1,2,3,4,5,6]; % different interface types: 1=Si, 2=Ge, 3=SiSiGe, 4=SiGeGe;
      for i=1:length(ainterface);
          for j=1:na;
              switch ainterface(i)
                  case 1;
                      atype(i,j)=1;   % 1=Si; 2=Ge;
                  case 2;   % fmatGe
                      atype(i,j)=2;   % 2=Ge;
                  case 3;   % fmatSiSiGe
                      if (x(j)<=5*A/4)
                          atype(i,j)=1;
                      else 
                          atype(i,j)=2;  % 2=Ge;
                      end
                  case 4;   % fmatSiGeGe
                      if (x(j)<A)
                          atype(i,j)=1;
                      else 
                          atype(i,j)=2;  % 2=Ge;
                      end
                  case 5;  % fmatGeGeSi
                      if (x(j)<=5*A/4)
                          atype(i,j)=2;
                      else 
                          atype(i,j)=1;  % 2=Ge;
                      end
                  case 6;  % fmatGeSiSi
                      if (x(j)<A)
                          atype(i,j)=2;
                      else 
                          atype(i,j)=1;  % 2=Ge;
                      end
              end
          end
      end
      % END different atom types at different interfaces.     
      
      %----------------------------------------
      %Calculate Phi & fmat;
      for i = 1:length(ainterface);
          switch ainterface(i)
              case 1;
                  [PhiSi, PhiSi111, fmatSi,D2,D3] = ComputeAllDmatMA(xs,ys,zs,ncx,ncy,ncz,A,atype(i,:), eps, isigma,leps, gamma);
              case 2;
                  [PhiGe, PhiGe111, fmatGe] = ComputeAllDmatMA(xs,ys,zs,ncx,ncy,ncz,A,atype(i,:), eps, isigma,leps, gamma);
              case 3;
                  [PhiSiSiGe, PhiSiSiGe111, fmatSiSiGe] = ComputeAllDmatMA(xs,ys,zs,ncx,ncy,ncz,A,atype(i,:), eps, isigma,leps, gamma);
              case 4;
                  [PhiSiGeGe, PhiGeGe111, fmatSiGeGe] = ComputeAllDmatMA(xs,ys,zs,ncx,ncy,ncz,A,atype(i,:), eps, isigma,leps, gamma);
              case 5;
                  [PhiGeGeSi, PhiGeGeSi111, fmatGeGeSi] = ComputeAllDmatMA(xs,ys,zs,ncx,ncy,ncz,A,atype(i,:), eps, isigma,leps, gamma);
              case 6;
                  [PhiGeSiSi, PhiGeSiSi111, fmatGeSiSi] = ComputeAllDmatMA(xs,ys,zs,ncx,ncy,ncz,A,atype(i,:), eps, isigma,leps, gamma);
          end
      end
% % % % % % %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
% % % % % % %------COMPARE MY FMAT WITH HONG ZHAO'S FMAT
% % % % % %       ASi=dlmread('fmat.Si');
% % % % % % AGe=dlmread('fmat.Ge');
% % % % % % ASiSiGe=dlmread('fmat.SiSiGe');
% % % % % % ASiGeGe=dlmread('fmat.SiGeGe');
% % % % % % 
% % % % % % 
% % % % % % n=0;
% % % % % % for i=-1:1
% % % % % % for j=-1:1;
% % % % % % for k= -1:1;
% % % % % % if(mod(i+j+k,2)~=0)
% % % % % % continue;
% % % % % % end
% % % % % % n=n+1;
% % % % % % fmatSiZHAO(:,:,i+2,j+2,k+2)=-ASi((n-1)*7+2:n*7,:);
% % % % % % fmatGeZHAO(:,:,i+2,j+2,k+2)=-AGe((n-1)*7+2:n*7,:);
% % % % % % fmatSiSiGeZHAO(:,:,i+2,j+2,k+2)=-ASiSiGe((n-1)*7+2:n*7,:);
% % % % % % fmatSiGeGeZHAO(:,:,i+2,j+2,k+2)=-ASiGeGe((n-1)*7+2:n*7,:);
% % % % % % %fmatGeGeSi(:,:,i+2,j+2,k+2)=-AGeGeSi((n-1)*7+2:n*7,:);
% % % % % % %fmatGeSiSi(:,:,i+2,j+2,k+2)=-AGeSiSi((n-1)*7+2:n*7,:);
% % % % % %  
% % % % % % end
% % % % % % end
% % % % % % end
% % % % % % 
% % % % % % dSi = fmatSi+fmatSiZHAO;
% % % % % % dGe = fmatGe+fmatGeZHAO;
% % % % % % dSiSiGe = fmatSiSiGe + fmatSiSiGeZHAO;
% % % % % % dSiGeGe = fmatSiGeGe + fmatSiGeGeZHAO; 
% % % % % % 
% % % % % % % %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% % % % % % % %PRINT OUT MY FMAT AND READ IN.  
% % % % % % % %ROUND MY FMAT TO %0.9E AS HONG ZHAO's format
% % % % % % % 
% % % % % % % 
% % % % % % % 
% % % % % % % fid1 = fopen('fmatSi_MA.txt', 'wt');
% % % % % % % fid2 = fopen('fmatGe_MA.txt', 'wt');
% % % % % % % fid3 = fopen('fmatSiSiGe_MA.txt', 'wt');
% % % % % % % fid4 = fopen('fmatSiGeGe_MA.txt', 'wt');
% % % % % % % [m,n,k,j,i]=size(fmatSi);
% % % % % % % for ii=1:1:i;
% % % % % % % for jj=1:1:j;
% % % % % % % for kk=1:1:k;
% % % % % % % for mm=1:1:m;
% % % % % % % for nn=1:1:n;
% % % % % % % if nn==n
% % % % % % % fprintf(fid1,'%0.9e\n', fmatSi(mm,nn,kk,jj,ii))
% % % % % % % fprintf(fid2,'%0.9e\n', fmatGe(mm,nn,kk,jj,ii))
% % % % % % % fprintf(fid3,'%0.9e\n', fmatSiSiGe(mm,nn,kk,jj,ii))
% % % % % % % fprintf(fid4,'%0.9e\n', fmatSiGeGe(mm,nn,kk,jj,ii))
% % % % % % % else
% % % % % % % fprintf(fid1,'%0.9e\t', fmatSi(mm,nn,kk,jj,ii))
% % % % % % % fprintf(fid2,'%0.9e\t', fmatGe(mm,nn,kk,jj,ii))
% % % % % % % fprintf(fid3,'%0.9e\t', fmatSiSiGe(mm,nn,kk,jj,ii))
% % % % % % % fprintf(fid4,'%0.9e\t', fmatSiGeGe(mm,nn,kk,jj,ii))
% % % % % % % end
% % % % % % % end
% % % % % % % end
% % % % % % % end
% % % % % % % end
% % % % % % % end
% % % % % % % fclose(fid1)
% % % % % % % fclose(fid2)
% % % % % % % fclose(fid3)
% % % % % % % fclose(fid4)
% % % % % % % 
% % % % % % % load('fmatSi_MA.txt');
% % % % % % % load('fmatGe_MA.txt');
% % % % % % % load('fmatSiSiGe_MA.txt');
% % % % % % % load('fmatSiGeGe_MA.txt');
% % % % % % % 
% % % % % % % clear fmatSi fmatGe fmatSiSiGe fmatSiGeGe
% % % % % % % 
% % % % % % % fmatSi=zeros(6,6,3,3,3);
% % % % % % % fmatGe=zeros(6,6,3,3,3);
% % % % % % % fmatSiSiGe=zeros(6,6,3,3,3);
% % % % % % % fmatSiGeGe=zeros(6,6,3,3,3);
% % % % % % % n=0;
% % % % % % % for i=1:3
% % % % % % % for j=1:3;
% % % % % % % for k= 1:3;
% % % % % % % %if(mod(i+j+k,2)~=0)
% % % % % % % %continue;
% % % % % % % %end
% % % % % % % n=n+1;
% % % % % % % fmatSi(:,:,k,j,i)=fmatSi_MA((n-1)*6+1:n*6,:);
% % % % % % % fmatGe(:,:,k,j,i)=fmatGe_MA((n-1)*6+1:n*6,:);
% % % % % % % fmatSiSiGe(:,:,k,j,i)=fmatSiSiGe_MA((n-1)*6+1:n*6,:);
% % % % % % % fmatSiGeGe(:,:,k,j,i)=fmatSiGeGe_MA((n-1)*6+1:n*6,:);
% % % % % % % 
% % % % % % % 
% % % % % % % end
% % % % % % % end
% % % % % % % end
% % % % % % % 
% % % % % % % dSi2 = fmatSi+fmatSiZHAO;
% % % % % % % dGe2 = fmatGe+fmatGeZHAO;
% % % % % % % dSiSiGe2 = fmatSiSiGe + fmatSiSiGeZHAO;
% % % % % % % dSiGeGe2 = fmatSiGeGe + fmatSiGeGeZHAO; 
% % % % % % % 
% % % % % % % %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% % % % % % % %END---PRINT OUT MY FMAT AND READ IN.  
% % % % % % 
% % % % % % save('fmatSi.mat','fmatSi')
% % % % % % save('fmatGe.mat','fmatGe')
% % % % % % save('fmatSiSiGe.mat','fmatSiSiGe')
% % % % % % save('fmatSiGeGe.mat','fmatSiGeGe')
% % % % % % save('fmatGeGeSi.mat','fmatGeGeSi')
% % % % % % save('fmatGeSiSi.mat','fmatGeSiSi')
% % % % % % %END SUBROUNTE 
% % % % % % 
% % % % % % 
