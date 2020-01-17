clear;
clc;
tic
%----------------------------------------------
    lambdaSi = 21.0; lambdaGe = 31.0; 
    epsSi = 3.473928e-19; epsGe = 3.085e-19; 
    sigmaSi = 2.0951e-10; sigmaGe = 2.0951e-10; 
   % massSi=28.085*1.6605e-27; massGe=72.64*1.6605e-27;
    massSi = 4.6637e-26; massGe = 1.2057e-25; 
    [eps,isigma,leps,epssig, imass] = SiGein(lambdaSi, lambdaGe, epsSi, epsGe, sigmaSi, sigmaGe, massSi, massGe);
    gamma=1.20d0;
    a0=5.43095e-10;
   % global lambdaSi lambdaGe epsSi epsGe sigmaSi sigmaGe massSi massGe eps isigma leps epssig gamma; 
%END Atom parameters. 

A = a0/sigmaSi; %  15168727030; 2.592215; % 2.592215168; %
ncell=2;
nbasis=4;
naxis=3;
forceflag=4;  % forceflag=1 : spring; forceflag=2: Silicon; forceflag=3: Germanium.  means choose s-w potential.

rc=1.8;

mI=1;
mTmatrix= massGe/massSi; % because of PhiGe has included the mass of massGe, so I just need reduced here as mT=1. 

%---------start--contstruct "force matrix" of pure Si or heavy Si-------
ncx = 3; % 5;
ncy = ncx; %5;
ncz = ncx; % 5;
ncxsuper = ncx;
ncxvoid = 1;
flagsupercell = 1;

N = 8;

boxlx=ncx*A;
boxly=ncy*A;
boxlz=ncz*A;
% [xs ys zs]=creatdiamond(A,ncx,ncy,ncz,ncell,nbasis);   
% [xs ys zs]=creatSuperCell(A,ncx,ncy,ncz,ncell,nbasis,ncxsuper);   
if(ncxvoid>1)
    [xs,ys,zs,nsupercell]=creatSuperCellwMultipleVoid(A,ncx,ncy,ncz,ncell,nbasis,ncxsuper,ncxvoid,flagsupercell);   
else
    [xs,ys,zs,nsupercell]=creatSuperCellwVoid(A,ncx,ncy,ncz,ncell,nbasis,ncxsuper,ncxvoid,flagsupercell);   
end
na=length(xs); % ncx*ncy*ncz*ncell*nbasis;

A1 = A*ncxsuper; % *ncx; % ; % 

% Phipure=zeros(naxis,naxis,ncell,na);
PhiSi=zeros(naxis,naxis,ncell,na);
PhiGe=zeros(naxis,naxis,ncell,na);
atypepure=zeros(2,na);
atypepure(1,:)=ones(1,na);
atypepure(2,:)=2*ones(1,na);  % here "2" means Ge. 

% -----------READ MA DATA-------
% % % load('C:\matlabcodes\Dispersion Curve\fmatSi.mat');
% % % load('C:\matlabcodes\Dispersion Curve\fmatGe.mat');
% % % load('C:\matlabcodes\Dispersion Curve\fmatSiSiGe.mat');
% % % load('C:\matlabcodes\Dispersion Curve\fmatSiGeGe.mat');
load('fmatSi.mat');
load('fmatGe.mat');
fmatSi = -fmatSi;
fmatGe = -fmatGe;
% % fmatSiSiGe = -fmatSiSiGe;
% % fmatSiGeGe = -fmatSiGeGe;

% for nmT=1:length(mTmatrix);%-0.3:0.4;
nmT=1; 
    clear transs reff checkk; 
    mT=mTmatrix(nmT); 
    nT=1;
    
Ksp=1.0;

ncell1 = ncell*nbasis;
nmode=naxis*ncell1;
wwindex=0;
tic
%%%%%%%%%%%%%%%%%%%%%%Dispersion Curve Calculation%%%%%%%%%%%%%%%%%%%%%
clear kx;


XX = ((0:1/N:1-1/N)+1/N/2)*1*pi/A1; %  (0.001:0.0625:1)*1*pi/(A1); 
YY = ((0:1/N:1-1/N)+1/N/2)*1*pi/A1; %(-1:1/32:1)*1*pi/(A1); %(-1:0.0625:1)*1*pi/(A1); % 0.0*kxx;% 
ZZ = ((0:1/N:1-1/N)+1/N/2)*1*pi/A1; %(-1:1/32:1)*1*pi/(A1); % (-1:0.0625:1)*1*pi/(A1); % 0.0*kxx;% (-1:0.0667:1)*1*pi/(A1); % 0.0*kxx; 
% dk = 0.001*2*pi/A/sigmaSi; 
dkx = (XX(2) - XX(1))/sigmaSi; 
dky = (YY(2) - YY(1))/sigmaSi; 
dkz = (ZZ(2) - ZZ(1))/sigmaSi; 


% % % dkk = [0.001; 0; 0]*1*pi/(A1);
 wwindex=0;
%  matlabpool local 2;
nxmax=length(XX);
nymax=length(YY);
nzmax=length(ZZ);

figure;
nk = 1;
for nx= 1:nxmax; % 1:25;  % nxmax;
    for ny= 1:nymax;
        for nz = 1:nzmax;
            if (ZZ(nz)<=YY(ny) && YY(ny)<=XX(nx) )
                kxx(nk) = XX(nx);
                kyy(nk) = YY(ny);
                kzz(nk) = ZZ(nz);
                if (YY(ny) == XX(nx) && ZZ(nz) == YY(ny))
%                     figure(11)
                    w(nk) = 1;
                    wx(nk) = 1; wy(nk) = 0; wz(nk) = 0;
%                      scatter3(kxx(nk),kyy(nk),kzz(nk),'r*')
%                      hold on;
                 else if(YY(ny) == XX(nx) && ZZ(nz) ~= YY(ny))
%                         figure(12)
                        w(nk) = 3;
                        wx(nk) = 2; wy(nk) = 0; wz(nk) = 1;
%                         scatter3(kxx(nk),kyy(nk),kzz(nk),'go')
%                         hold on;
                     else if(ZZ(nz) == XX(nx) && ZZ(nz) ~= YY(ny))
%                             figure(13)
                            w(nk) = 3;
                            wx(nk) = 1; wy(nk) = 1; wz(nk) = 1;
%                             scatter3(kxx(nk),kyy(nk),kzz(nk),'bs')
%                             hold on;
                          else if(ZZ(nz) == YY(ny) && ZZ(nz) ~= XX(nx))
%                                   figure(14)
                                  w(nk) = 3;
                                   wx(nk) = 1; wy(nk) = 2; wz(nk) = 0;
%                                    scatter3(kxx(nk),kyy(nk),kzz(nk),'k.')
%                                     hold on;
                              else
%                                   figure(15)
                                    w(nk) = 2*3;
%                                     wx(nk) = 3; wy(nk) = 3; wz(nk) = 0;
                                    wx(nk) = 2*1; wy(nk) = 2*1; wz(nk) = 2*1;
%                                     scatter3(kxx(nk),kyy(nk),kzz(nk),'m+')
%                                     hold on;
                              end
                          end
                     end
                end
               nk = nk+1;

            else
                continue;
                
            end

        
        end
    end
end

% % % % figure;
% % % % scatter3(kxx,kyy,kzz)
% % % % xlabel('kx');ylabel('ky');zlabel('kz')
% % % % xlim([0 pi/A1]); ylim([0 pi/A1]); zlim([0 pi/A1]); 

% % % dkk = [0.001; 0; 0]*1*pi/(A1);
 wwindex=0;
%  matlabpool local 2;
nxmax=length(kxx);
nymax=length(kyy);
nzmax=length(kzz);

kindex = nxmax;
nbranch = ncell*nbasis*(ncxsuper^3 - ncxvoid)*3;
% nbranch = 3*ncell*nbasis*(ncxsuper^3 - ncxvoid); 
% for kx=ks:dk:kend;
%     kx=kxx(1); ky=kyy(1); kz=kzz(1);
%     kvI=[kx;ky;kz];  % on the [111] plane or direction. 
%     [wwI1,DI1,eV1]=dispersioncurveSuperCell(kvI,A,ncell,rc,Ksp,naxis,mI,nbasis,forceflag,fmatSi,xs,ys,zs,ncx,ncy,ncz,ncxsuper,nsupercell);  
% 
% parpool(4)
% kindex = nxmax;
% nbranch = 3*ncell*nbasis*(ncxsuper^3 - ncxvoid); 
[Phiall, rxall, ryall, rzall]=dispersioncurveSuperCell_Modify(A,ncell,rc,Ksp,naxis,mI,nbasis,forceflag,fmatSi,xs,ys,zs,ncx,ncy,ncz,ncxsuper,nsupercell);  

ninc = 0
Cpsum = 0;
ksum = 0;

h = 1.054571726e-34; % reduced planck constant. 
kB = 1.3806488e-23;
% % % Davis 2011
B = 2.1e-19; 
C = 180;
D = 1.32e-45;

% % % % % % Hopkins 2010
% % % B = 3.73e-19; % unit: sK^-1;
% % % C = 157.3; % unit: K.
% % % D = 9.32e-45; % unit: s^3;
% % % % E = 2.3e-3; % unit: m. 

n=0;
% %  matlabmail('ruiyuanma@gmail.com','Alarm','M62 Code Start','townem62@gmail.com', '19580101')

    T= 300%  5:5:400 % :4150;  
for nk= 1:length(kxx); % 1:25;  % nxmax;
% % %     n = n+1;
% % %     for ny= 1:nymax;
    if(mod(nk,10)==0)
      disp(nk/length(kxx));
    end
        data = zeros(nbranch);
        vmodex = zeros(nbranch);
        vmodey = zeros(nbranch);
        vmodez = zeros(nbranch);

% % %         for nz = 1:nzmax;
% % for nx = 1:nxmax;
% %     kx = kxx(nx);
% %     ky = 0;
% %     kz = 0;
    kx=kxx(nk);
    ky=kyy(nk);
    kz=kzz(nk);
    Dsuper = zeros(naxis*nsupercell, naxis*nsupercell);
    Dsuperf = zeros(naxis*nsupercell, naxis*nsupercell);
    Dsfx = zeros(naxis*nsupercell, naxis*nsupercell);
    Dsbx = zeros(naxis*nsupercell, naxis*nsupercell);
    Dsfy = zeros(naxis*nsupercell, naxis*nsupercell);
    Dsby = zeros(naxis*nsupercell, naxis*nsupercell);
    Dsfz = zeros(naxis*nsupercell, naxis*nsupercell);
    Dsbz = zeros(naxis*nsupercell, naxis*nsupercell);
    
kxall = kx*ones(naxis*nsupercell,naxis*na);
kyall = ky*ones(naxis*nsupercell,naxis*na);
kzall = kz*ones(naxis*nsupercell,naxis*na);
Dall = Phiall.*exp(1i*(kxall.*rxall+kyall.*ryall+kzall.*rzall));
% % % Dall = Phiall.*exp(1i*(kx*rxall+ky*ryall+kz*rzall));

for i = 1:na/nsupercell;
Dsuper = Dsuper + Dall(:,(i-1)*naxis*nsupercell+(1:naxis*nsupercell));
end

% % %  [UU,SS,VV]=svdecon(Dsuper);
 [UU,SS]=eig(Dsuper);
% % %  eV=zeros(naxis*nsupercell,naxis*nsupercell);
 eV=UU ; % 
 wwI = sqrt(diag(SS));


% % % ### Velocity #3: 
        dk = 0.0001*1*pi/A1;
    Dfx = (Dall.*exp(1i*(dk)*rxall)); 
    Dbx = (Dall.*exp(1i*(-dk)*rxall)); 
    Dfy = (Dall.*exp(1i*(dk)*ryall)); 
    Dby = (Dall.*exp(1i*(-dk)*ryall)); 
    Dfz = (Dall.*exp(1i*(dk)*rzall)); 
    Dbz = (Dall.*exp(1i*(-dk)*rzall)); 
%     check = max(Dallf-Df2);
    for i = 1:na/nsupercell;
        Dsfx = Dsfx + Dfx(:,(i-1)*naxis*nsupercell+(1:naxis*nsupercell));
        Dsbx = Dsbx + Dbx(:,(i-1)*naxis*nsupercell+(1:naxis*nsupercell));
        Dsfy = Dsfy + Dfy(:,(i-1)*naxis*nsupercell+(1:naxis*nsupercell));
        Dsby = Dsby + Dby(:,(i-1)*naxis*nsupercell+(1:naxis*nsupercell));   
        Dsfz = Dsfz + Dfz(:,(i-1)*naxis*nsupercell+(1:naxis*nsupercell));
        Dsbz = Dsbz + Dbz(:,(i-1)*naxis*nsupercell+(1:naxis*nsupercell)); 
    end
% % %     vmode(nz,:) = diag(eV'*((Dsfx-Dsuper)*(1/massSi)/(dkx/sigmaSi))*eV)./(2*wwI*sqrt(1/massSi));
    vmodex = diag(eV'*((Dsfx-Dsbx)*(1/massSi)/(2*dk/sigmaSi))*eV)./(2*wwI*sqrt(1/massSi));
    vmodey = diag(eV'*((Dsfy-Dsby)*(1/massSi)/(2*dk/sigmaSi))*eV)./(2*wwI*sqrt(1/massSi));
    vmodez = diag(eV'*((Dsfz-Dsbz)*(1/massSi)/(2*dk/sigmaSi))*eV)./(2*wwI*sqrt(1/massSi));


% % % ### Velocity #3: END

%         end
        
        
    
  
    
    
     data = wwI*sqrt(1/massSi); 
%         end
        Cpmode = kB*(h*data/kB/T).^2.*(exp(h*data/kB/T)./((exp(h*data/kB/T)-1).^2));
        Cpsum = Cpsum + sum(sum(w(nk)*Cpmode)); % kB*(h*data/kB/T).^2.*(exp(h*data/kB/T)./((exp(h*data/kB/T)-1).^2))));
        
        
        wwmode = data;
            umk_scat_inv = B*(wwmode.*wwmode)*T*exp(-C/T);
    % impurity scattering. 
    imp_scat_inv = D*wwmode.^4;
    % boundary scattering.
    bound_scat_inv = 0;% (vmode)/E; % (1-p)/(1+p)*(1/dx + 1/dy)*abs(vmode)/2;
    
    tot_scat_inv = umk_scat_inv + imp_scat_inv  + bound_scat_inv; % ;%
    tot_scat = 1./tot_scat_inv;
    kmode_holland = wx(nk)*vmodex.*vmodex.*tot_scat.*Cpmode + wy(nk)*vmodey.*vmodey.*tot_scat.*Cpmode + wz(nk)*vmodez.*vmodez.*tot_scat.*Cpmode; %
    ksum = ksum + sum(sum(kmode_holland));
    ksum2 = sum(sum(kmode_holland));
%     end
    ksumx(nk) = ksum;
    ksumx2(nk) = ksum2;
    
    
end
% % % Subroutine: Calculate thermal conductivity.
% Method #1! Be careful about the Volume for Primitive and Super Cell! 
% % % Alos be careful are 8 time difference for Primitive and Conventional
% Cell. 
k = ksum*1/(N*N*N*(A1*sigmaSi)^3) %  Super Cell/Conventional Cell Only. 
% % % k2 = ksum*dkx*dky*dkz/(2*pi)^3 %  % Primitive Cell Only  
k2 = 8*ksum*dkx*dky*dkz/(2*pi)^3 % Super Cell/Conventional Cell Only.
Cp = Cpsum/(N*N*N*(A1*sigmaSi)^3)/2329 %  Super Cell/Conventional Cell Only. 
% % % Cp = Cpsum*8/(N*N*N*(A1*sigmaSi)^3)/2329 % Primitive Cell Only
% % % Following Equations are capable for Both Super Cell and Primitive
% Cell. 
Lx = 2*pi/(dkx); Ly = 2*pi/(dky); Lz = 2*pi/(dkz);
V = Lx*Ly*Lz;
k3 = (8*ksum)/V 
Cp3 = Cpsum*8/V/2329 % 

% % % Subroutine: End

timevalue = toc

figure
plot(kxx,ksumx*dkx*dky*dkz/(2*pi)^3,'o-')

xlswrite('kxx-N10-ncx6-ncxvoid2.xlsx',kxx)
xlswrite('ksumx-N10-ncx6-ncxvoid2.xlsx',ksumx)
xlswrite('ksumx2-N10-ncx6-ncxvoid2.xlsx',ksumx2)


% End 

% % % figure
% % % plot(N1,K1,'ko-','LineWidth',2,'MarkerSize',10)
% % % hold on
% % % plot(N2,K2,'ro-','LineWidth',2,'MarkerSize',10)
% % % hold on
% % % plot(N3,K3,'bo-','LineWidth',2,'MarkerSize',10)
% % % hold on
% % % plot(N4,K4,'go-','LineWidth',2,'MarkerSize',10)
% % % xlabel('N (number of k steps in one direction)')
% % % ylabel('k (W/mK)')
% % % legend('N=1','N=2','N=3','N=4')
% % % set(gca,'Linewidth',3.0)
% % % 
% % % figure
% % % plot(N1,T1,'ko-','LineWidth',2,'MarkerSize',10)
% % % hold on
% % % plot(N2,T2,'ro-','LineWidth',2,'MarkerSize',10)
% % % hold on
% % % plot(N3,T3,'bo-','LineWidth',2,'MarkerSize',10)
% % % hold on
% % % plot(N4,T4,'go-','LineWidth',2,'MarkerSize',10)
% % % xlabel('N (number of k steps in one direction)')
% % % ylabel('Time (s)')
% % % legend('N=1','N=2','N=3','N=4')
% % % set(gca,'Linewidth',3.0)

% % % % figure
% % % % plot(kxx,ksumx2*dkx*dky*dkz/(2*pi)^3,'o-')
% % % % ksumx2*dkx*dky*dkz/(2*pi)^3

% % % % % % % matlabmail('ruiyuanma@gmail.com','k','M62 Code Finish','townem62@gmail.com', '19580101')
