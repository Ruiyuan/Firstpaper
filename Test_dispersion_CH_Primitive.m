% % % # # # DISPERSION CURVE CODE CAN SWITCH BETWEEN PRIMITIVE CELL AND
% SUPERCELL # # # % % %
% % % Add a "flagsuper" switch in this main code, and
% "creatSuperCellwVoid.m" code. 
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
mTmatrix=[massGe/massSi]; % because of PhiGe has included the mass of massGe, so I just need reduced here as mT=1. 

%---------start--contstruct "force matrix" of pure Si or heavy Si-------
ncx = 2; % 5;
ncy = ncx; %5;
ncz = 2; % 5;
% % % ncxh=3;
% % % ncyh=3;
% % % nczh=3;

ncxsuper = ncx;
% % % ncysuper = ncy;
% % % nczsuper = ncz+nczh;
ncxvoid = 0;

% % % SWITCH BETWEEN PRIMITIVE CELL AND SUPERCELL!!!
flagsupercell = 0;
if(flagsupercell == 1)
    A1 = A*ncxsuper; % *ncx; % ; %
else 
    A1 = A*1;
end
% % % SWITCH BETWEEN PRIMITIVE CELL AND SUPERCELL!!!

boxlx=ncx*A;
boxly=ncy*A;
boxlz=ncz*A;
% [xs ys zs]=creatdiamond(A,ncx,ncy,ncz,ncell,nbasis);   
% [xs ys zs]=creatSuperCell(A,ncx,ncy,ncz,ncell,nbasis,ncxsuper);   

[xs,ys,zs,nsupercell]=creatSuperCellwVoid(A,ncx,ncy,ncz,ncell,nbasis,ncxsuper,ncxvoid,flagsupercell);   
na=length(xs); % ncx*ncy*ncz*ncell*nbasis;
figure;
scatter3(xs,ys,zs);

% Phipure=zeros(naxis,naxis,ncell,na);
% % % PhiSi=zeros(naxis,naxis,ncell,na);
% % % PhiGe=zeros(naxis,naxis,ncell,na);
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
kxx =  (0.001:0.01:1)*2*pi/(A1); %
kyy = 0; % 0.0*kxx; % (0.001:0.1:1-0.001)*1*pi/(A1); % 0.5*1*pi/A1; %  
kzz = 0; %  0.0*kxx; % 0.1*1*pi/A1; % %  
% dk = 0.001*2*pi/A/sigmaSi; 
% dkx = (kxx(2) - kxx(1))/sigmaSi; 
% dky = (kyy(2) - kyy(1))/sigmaSi; 
% dkz = (kzz(2) - kzz(1))/sigmaSi; 

% % % dkk = [0.001; 0; 0]*1*pi/(A1);
 wwindex=0;
%  matlabpool local 2;
nxmax=length(kxx);
nymax=length(kyy);
nzmax=length(kzz);

kindex = nxmax;
nbranch = length(xs)*3; % ncell*nbasis*(ncxsuper^3*3 - ncxvoid);
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
B = 2.1e-19; 
C = 180;
D = 1.32e-45;
n=0;
%  matlabmail('ruiyuanma@gmail.com','Alarm','M62 Code Start','townem62@gmail.com', '19580101')

    T= 300%  5:5:400 % :4150;  
% % % for nz= 1:nzmax; % 1:25;  % nxmax;
% % %     n = n+1
% % %     for nz= 1:nzmax
% % %         data = zeros(nxmax,nbranch);
% % %         for nx = 1:nxmax;
    ky = kyy;
%     ky = kyy;
    kz = kzz;
for nx = 1:nxmax;
% for nz = 1:nzmax;
    kx = kxx(nx);
    ky = 1*kx; % kyy;
    kz = 1*kx; % kzz;
%     kz=kzz(nz);
%     kz=kzz(nz);
    ninc = ninc + 1;
    Dsuper = zeros(naxis*nsupercell, naxis*nsupercell);
    Dsuperf = zeros(naxis*nsupercell, naxis*nsupercell);
    Dsf2 = zeros(naxis*nsupercell, naxis*nsupercell);
    

kxall = kx*ones(naxis*nsupercell,naxis*na);
kyall = ky*ones(naxis*nsupercell,naxis*na);
kzall = kz*ones(naxis*nsupercell,naxis*na);
Dall = Phiall.*exp(1i*(kxall.*rxall+kyall.*ryall+kzall.*rzall));
% % % Dall = Phiall.*exp(1i*(kx*rxall+ky*ryall+kz*rzall));

% % % Dsuper = Dall;
for i = 1:na/nsupercell;
Dsuper = Dsuper + Dall(:,(i-1)*naxis*nsupercell+(1:naxis*nsupercell));
end

 [UU,SS,VV]=svdecon(Dsuper);
% % %  eV=zeros(naxis*nsupercell,naxis*nsupercell);
 eV=UU ; % 
 wwI = sqrt(diag(SS));
% % %  for i=1:naxis*nsupercell;
% % %     om(i,1)=SS(i,i);
% % %  end
% % %  wwI=(sqrt(om));   
% % %  wwxx2(nx,nz,nz,:)=wwI*sqrt(1/massSi); %
 
 % % %### Dispersion Values at a smaller forward. 
 kxallf = (kx +0.001*1*pi/(A1))*ones(naxis*nsupercell,naxis*na);
 Dallf = Phiall.*exp(1i*(kxallf.*rxall+kyall.*ryall+kzall.*rzall));
% % % %  Dsuperf = Dallf;
    for i = 1:na/nsupercell;
        Dsuperf = Dsuperf + Dallf(:,(i-1)*naxis*nsupercell+(1:naxis*nsupercell));
    end
 [UU,SS,VV]=svdecon(Dsuperf);
% % % %  eVf=zeros(naxis*nsupercell,naxis*nsupercell);
 eVf=UU ; % 
  wwf=(sqrt(diag(SS)));  
% % % %  for i=1:naxis*nsupercell;
% % % %     omf(i,1)=SS(i,i);
% % % %  end
% % % %  wwf=(sqrt(omf));  
% % % % % % % % % ### GROUP VELOCITY ######## 
% % % %     checkeVf = abs((eV'*eVf)./(norm(eV).*norm(eVf')));
% % % %     [Mf, ipf] = max(checkeVf');
% % % % % %     wwxxf(nx,nz,nz,:)=wwf(ipf)*sqrt(1/massSi);
% % % % % %     vmode(nx,nz,nz,:) = (wwI - wwf(ipf))*sqrt(1/massSi)/(0.001*1*pi/A1/sigmaSi);
% % % %     vmode(nx,:) = (wwI - wwf(ipf))*sqrt(1/massSi)/(0.001*1*pi/A1/sigmaSi);
% % % % % % % ### Array-wise Operation: The END.
% % % % 
% % % %     data(nx,:) = wwI*sqrt(1/massSi); 
% % % %     
% % % %     % % % Potential Code for gruop velocity calculation. 
% % % %         dkx = 0.001*1*pi/A1;
% % % %     Df2 = (Dall.*exp(1i*(dkx)*rxall)); 
% % % %     check = max(Dallf-Df2);
% % % %     for i = 1:na/nsupercell;
% % % %         Dsf2 = Dsf2 + Df2(:,(i-1)*naxis*nsupercell+(1:naxis*nsupercell));
% % % %     end
% % % %     
% % % % % % % % %     Dcross = (Dsuper'*Dsf2);
% % % % % % % % %    
% % % % % % % % %             [UU,SS,VV]=svdecon(Dsuper);
% % % % % % % % %         [Uf,Sf,Vf]=svdecon(Dsuperf);
% % % % % % % % %         Dcross = Dsuper'*Dsuperf;
% % % % % % % % %         Scross =SS'*Sf;
% % % % % % % % %             [Mc, ipc] = max(Dcross);
% % % % % % % % %             [Ms, ips] = max(Scross);
% % % % % % % % %          [Uc,Sc,Vc] = svdecon(Dcross);
% % % % % % % % %     for i =1:nbranch;
% % % % % % % % %         Dcross(:,i) = abs(Dcross(:,i)./(SS(i,i)^2));
% % % % % % % % %     end
% % % %         for i=1:nbranch;
% % % % %             vg(nx,i) = eV(:,i)'*((Dall*(1/massSi).*exp(1i*(dkx)*rxall)-Dall*(1/massSi))/(dkx))*eV(:,i)/(2*wwI(i)*sqrt(1/massSi));
% % % %               vg(nx,i) = eV(:,i)'*((Dsf2-Dsuper)*(1/massSi)/(dkx/sigmaSi))*eV(:,i)/(2*wwI(i)*sqrt(1/massSi));
% % % %         end
% % % %             vg2(nx,:) = diag(eV'*((Dsf2-Dsuper)*(1/massSi)/(dkx/sigmaSi))*eV)./(2*wwI*sqrt(1/massSi));
% % % % % % % % % ### GROUP VELOCITY ######## 
           
wwxx(nx,:) = wwI*sqrt(1/massSi);  % wwxx(nz,nx,:) = wwI; 

        
% end

end

figure
plot(kxx,wwxx,'k-');
% ylim([0 5])

% % % % % for i = nbranch:-1:1
% % % % %     figure(11)
% % % % % surf(kzz,kyy,wwxx(:,:,i)*sqrt(1/massSi))
% % % % % % pause;
% % % % % hold on;
% % % % % % % % Check Symmetry:
% % % % % d_sym(:,:,i) = (wwxx(:,:,i) - wwxx(:,:,i)')/wwxx(:,:,i);
% % % % %     figure(12)
% % % % % surf(kzz,kyy,d_sym(:,:,i))
% % % % % hold on;
% % % % % end
% % % % % figure(11)
% % % % % xlabel('kz'); ylabel('ky'); zlabel('Frequency (rad/ps)')
% % % % % figure(12)
% % % % % xlabel('kz'); ylabel('ky'); zlabel('w-transpose(w)/w')
% % figure
% % plot(kxx,data,'k-')
% % figure
% % plot(kxx,vmode,'b.')
% % hold on;
% % plot(kxx,vg,'ro')
% % plot(kxx,vg2,'go')

