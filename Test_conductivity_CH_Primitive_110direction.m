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
ncx =  6;
ncy = ncx; %5;
ncz = ncy; % 5;
% % % ncxh=3;
% % % ncyh=3;
% % % nczh=3;

ncxsuper =  ncx; % ncx; % ncx; %  ncx; % ncxsuper must >= 2 for Primitive Cell! 
% % % ncysuper = ncy;
% % % nczsuper = ncz+nczh;
ncxvoid = 0; % 1;

% % % SWITCH BETWEEN PRIMITIVE CELL AND SUPERCELL!!!
flagsupercell = 1; % 1;
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
% % % figure;
% % % scatter3(xs,ys,zs);

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
N = 40; 
if(flagsupercell == 1)
    XX = ((0:1/N:1-1/N)+1/N/2)*1*pi/A1; % sqrt(2); %2; %   (0.001:0.0625:1)*1*pi/(A1); 
else 
    XX = ((0:1/N:1-1/N)+1/N/2)*2*pi/A1; %
end
% % YY = ((0:1/N:1-1/N)+1/N/2)*2*pi/A1; %(-1:1/32:1)*1*pi/(A1); %(-1:0.0625:1)*1*pi/(A1); % 0.0*kxx;% 
% % ZZ = ((0:1/N:1-1/N)+1/N/2)*2*pi/A1; %(-1:1/32:1)*1*pi/(A1); % (-1:0.0625:1)*1*pi/(A1); % 0.0*kxx;% (-1:0.0667:1)*1*pi/(A1); % 0.0*kxx; 
% dk = 0.001*2*pi/A/sigmaSi; 
dkx = (XX(2) - XX(1))/sigmaSi; 
% % dky = (YY(2) - YY(1))/sigmaSi; 
% % dkz = (ZZ(2) - ZZ(1))/sigmaSi; 

mI = ones(1,na);
% % % dkk = [0.001; 0; 0]*1*pi/(A1);
 wwindex=0;
%  matlabpool local 2;
nxmax=length(XX);
% % nymax=length(YY);
% % nzmax=length(ZZ);

% % % figure;
% % % nk = 1;
% % % for nx= 1:nxmax; % 1:25;  % nxmax;
% % %     for ny= 1:nymax;
% % %         for nz = 1:nzmax;
% % %             if ((XX(nx) + ZZ(nz) + YY(ny))> (3/2*2*pi/A1))
% % %                continue;
% % %             else
% % %                kxx(nk) = XX(nx);
% % %                kyy(nk) = YY(ny);
% % %                kzz(nk) = ZZ(nz);
% % %                nk = nk+1;
% % %             end 
% % %         end
% % %     end
% % % end

 wwindex=0;
%  matlabpool local 2;
kxx = XX;
nxmax=length(kxx);
% % % nymax=length(kyy);
% % % nzmax=length(kzz);

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


h = 1.054571726e-34; % reduced planck constant. 
kB = 1.3806488e-23;
% % % % B = 2.1e-19; 
% % % % C = 180;
% % % % D = 1.32*10^(-45);
% % % % E = 1e20; % L;
% % % % n=0;

B = 1.73e-19; 
C = 137.39;
D = 1.32*10^(-45);
E = 1e20; % L;
n=0;

% % % % % % % Hopkins-2009
% % % % B = 3.73e-19; % unit: sK^-1;
% % % % C = 157.3; % unit: K.
% % % % D = 9.32e-45; % unit: s^3;
% % % % E = 2.3e-3; % unit: m. 

%  matlabmail('ruiyuanma@gmail.com','Alarm','M62 Code Start','townem62@gmail.com', '19580101')

    
% % % for nz= 1:nzmax; % 1:25;  % nxmax;
% % %     n = n+1
% % %     for nz= 1:nzmax
% % %         data = zeros(nxmax,nbranch);
% % %         for nx = 1:nxmax;
% % %     ky = kyy;
%     ky = kyy;
% % %     kz = kzz;
if(ncxvoid == 0)
            L = 1e20;
        else
            L = (ncx/ncxvoid - ncxvoid)*a0;
end
TT = 300; % 40:10:400; %
for ii = 1:length(TT); %   5:5:400 % :4150;  
    T = TT(ii);
    ninc = 0
Cpsum = 0;
ksum = 0;

 ksumx_NBS = zeros(length(XX),1);
ksumx_3Part = zeros(length(XX),1);
ww = zeros(length(XX),nbranch);       


for nk = 1:length(XX);
        data = zeros(nbranch);
        vmodex = zeros(nbranch);
        vmodey = zeros(nbranch);
        vmodez = zeros(nbranch);
    if(mod(nk,10) == 0)
        disp(nk/length(XX));
    end
    kx=kxx(nk);
    ky= kx; %  kyy(nk);
    kz= 0; % kzz(nk);
    kall = sqrt(kx^2 + ky^2 + kz^2);
    dk = dkx*(kall/kx);
    ninc = ninc + 1;
    
% % %     Dsuper = zeros(naxis*nsupercell, naxis*nsupercell);
% % %     Dsuperf = zeros(naxis*nsupercell, naxis*nsupercell);
% % %     Dsuperb = zeros(naxis*nsupercell, naxis*nsupercell);

% % % %     Dsfx = zeros(naxis*nsupercell, naxis*nsupercell);
% % % %     Dsbx = zeros(naxis*nsupercell, naxis*nsupercell);
% % % %     Dsfy = zeros(naxis*nsupercell, naxis*nsupercell);
% % % %     Dsby = zeros(naxis*nsupercell, naxis*nsupercell);
% % % %     Dsfz = zeros(naxis*nsupercell, naxis*nsupercell);
% % % %     Dsbz = zeros(naxis*nsupercell, naxis*nsupercell);

% % % kxall = kx*ones(naxis*nsupercell,naxis*na);
% % % kyall = ky*ones(naxis*nsupercell,naxis*na);
% % % kzall = kz*ones(naxis*nsupercell,naxis*na);
% Dall = Phiall.*exp(1i*(kxall.*rxall+kyall.*ryall+kzall.*rzall));
% % % Dall = Phiall.*exp(1i*(kx*rxall+ky*ryall+kz*rzall));

    Dsuper = zeros(naxis*nsupercell, naxis*nsupercell);

Dall = Phiall.*exp(1i*(kx*rxall+ky*ryall+kz*rzall));

% % % Dsuper = Dall;
for i = 1:na/nsupercell;
Dsuper = Dsuper + Dall(:,(i-1)*naxis*nsupercell+(1:naxis*nsupercell));
end

  [eV,SS]=svd(Dsuper);
  wwI = sqrt(diag(SS));
%   ww(nk,:) = wwI;
%   continue;
% % %  for i=1:naxis*nsupercell;
% % %     om(i,1)=SS(i,i);
% % %  end
% % %  wwI=(sqrt(om));   
% % %  wwxx2(nx,nz,nz,:)=wwI*sqrt(1/massSi); %
  
% % % % % % %  % % %### Dispersion Values at a smaller forward. 
% % % % % % %  kxallf = (kx +0.001*1*pi/(A1))*ones(naxis*nsupercell,naxis*na);
% % % % % % %  kxallb = (kx -0.001*1*pi/(A1))*ones(naxis*nsupercell,naxis*na);
% % % % % % %  kyallf = ky; % (ky +0.001*1*pi/(A1))*ones(naxis*nsupercell,naxis*na);
% % % % % % %  kyallb = ky; % (ky -0.001*1*pi/(A1))*ones(naxis*nsupercell,naxis*na);
% % % % % % %  
% % % % % % %  Dallf = Phiall.*exp(1i*(kxallf.*rxall+kyallf.*ryall+kzall.*rzall));
% % % % % % %  Dallb = Phiall.*exp(1i*(kxallb.*rxall+kyallb.*ryall+kzall.*rzall));
% % % % % % % 
% % % % % % % % % % %  Dsuperf = Dallf;
% % % % % % %     for i = 1:na/nsupercell;
% % % % % % %         Dsuperf = Dsuperf + Dallf(:,(i-1)*naxis*nsupercell+(1:naxis*nsupercell));
% % % % % % %         Dsuperb = Dsuperb + Dallb(:,(i-1)*naxis*nsupercell+(1:naxis*nsupercell));
% % % % % % %     end
% % % % % % %  [UU,SS,VV]=svd(Dsuperf);
% % % % % % % % % % %  eVf=zeros(naxis*nsupercell,naxis*nsupercell);
% % % % % % %  eVf=UU ; % 
% % % % % % %   wwf=(sqrt(diag(SS)));  
% % % % % % %   [UU,SS,VV]=svd(Dsuperb);
% % % % % % % % % % %  eVf=zeros(naxis*nsupercell,naxis*nsupercell);
% % % % % % %  eVb=UU ; % 
% % % % % % %   wwb=(sqrt(diag(SS))); 
   dkxf = 0.001*1*pi/A1; % 
 % % %### Dispersion Values at a smaller forward. 
    
    Dsuper = zeros(naxis*nsupercell, naxis*nsupercell);
    kxf = (kx+dkxf); % (kx + 0.025*1*pi/(A1)); %(kx+0.01); %  
    kyf = ky +dkxf; % kxf; 
    kzf = kz; % +dkxf; 
    Dall = Phiall.*exp(1i*(kxf*rxall+kyf*ryall+kzf*rzall)); 
  
    for i = 1:na/nsupercell;
        Dsuper = Dsuper + Dall(:,(i-1)*naxis*nsupercell+(1:naxis*nsupercell));
    end
  [eVf,SS]=svd(Dsuper);
    wwf=(sqrt(diag(SS))); 
    % % % ***Backward Step. **********************    
    Dsuper = zeros(naxis*nsupercell, naxis*nsupercell);    
    kxb = (kx - dkxf);
    kyb = ky  - dkxf; % kxb; 
    
    Dall = Phiall.*exp(1i*(kxb*rxall+kyb*ryall+kz*rzall));
    for i = 1:na/nsupercell;
        Dsuper = Dsuper + Dall(:,(i-1)*naxis*nsupercell+(1:naxis*nsupercell));
    end 
    [eVb,SS]=svd(Dsuper);
    wwb=(sqrt(diag(SS))); 
% % % ***Backward Step. END********************** 

% % % % % ### GROUP VELOCITY ***Center Diff***######## 
    checkeV = abs((eVb'*eVf));% ./(norm(eVb).*norm(eVf')));
    [M, ipf] = max(checkeV');
%             clear checkeV M;
%     checkeV = []; M = []; 
    dkxvg = sqrt((kxf - kxb)^2 + (kyf - kyb)^2 );
    vmodef = (wwb - wwf(ipf))*sqrt(1/massSi)/(1*dkxvg/sigmaSi); % (1*0.001*kx/sigmaSi); % (1*0.025*1*pi/A1/sigmaSi); % (1*0.01/sigmaSi); % 
    vmode = vmodef;
    % % % % % ### GROUP VELOCITY ***RORWARD***######## 

% % % % % % % % ### GROUP VELOCITY ***FORWARD***######## 
% % %     checkeV = abs((eV'*eVf));% ./(norm(eVb).*norm(eVf')));
% % %     [M, ipf] = max(checkeV');
% % % %             clear checkeV M;
% % % %     checkeV = []; M = []; 
% % %     dkxvg = sqrt((kxf - kx)^2 + (kyf - ky)^2 );
% % %     vmodef = (wwI - wwf(ipf))*sqrt(1/massSi)/(1*dkxvg/sigmaSi); % (1*0.001*kx/sigmaSi); % (1*0.025*1*pi/A1/sigmaSi); % (1*0.01/sigmaSi); % 
% % %     vmode = vmodef;
% % %     % % % % % ### GROUP VELOCITY ***FORWARD***######## 

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
           
% % wwxx(nx,:) = wwI*sqrt(1/massSi);  % wwxx(nz,nx,:) = wwI; 

      wwmode = wwI*sqrt(1/massSi); 
%         end
        Cpmode = kB*(h*wwmode/kB/T).^2.*(exp(h*wwmode/kB/T)./((exp(h*wwmode/kB/T)-1).^2));
        Cpsum = Cpsum + sum(sum((Cpmode*(kall/sigmaSi)^2))); % kB*(h*data/kB/T).^2.*(exp(h*data/kB/T)./((exp(h*data/kB/T)-1).^2))));
        
        
            umk_scat_inv = B*(wwmode.*wwmode)*T*exp(-C/T);
    % impurity scattering. 
    imp_scat_inv = D*wwmode.^4;
    % boundary scattering.
    bound_scat_inv = (vmode)/E; % L; % E; %  (vmode)/E; %0; %  (1-p)/(1+p)*(1/dx + 1/dy)*abs(vmode)/2;
 % partical scattering. from Majumdar-1993
    rvoid = ncxvoid*a0/2;
    chi = rvoid*sqrt(kx^2+ky^2+kz^2)/sigmaSi;
    sig_void = pi*rvoid^2*(chi^4/(chi^4+1));
    rho_void = (ncxvoid/ncx)^3;
    part_scat_inv = sig_void*rho_void*abs(vmode);
    
    tot_scat_inv = umk_scat_inv + imp_scat_inv  + bound_scat_inv; % ;%
    tot_scat = 1./tot_scat_inv;
    kmode_holland = vmode.*vmode.*tot_scat.*Cpmode; %  + wy(nk)*vmodey.*vmodey.*tot_scat.*Cpmode + wz(nk)*vmodez.*vmodez.*tot_scat.*Cpmode; %
    ksum = ksum + sum(sum(kmode_holland.*(kall/sigmaSi)^2));
%     ksum2 = sum(sum(kmode_holland));
%     end
% % %     ksumx(nk) = ksum;
% % %     ksumx2(nk) = ksum2;
% % %     wwxx(nk,:) = wwI;
 % % % *** Don't consider Boundary Scattering. ***
    tot_scat_inv2 = umk_scat_inv + imp_scat_inv;
    tot_scat2 = 1./tot_scat_inv2;
     kmode_holland2_NBS = vmode.*vmode.*tot_scat2.*Cpmode; %
     ksumx_NBS(nk) = sum(sum(kmode_holland2_NBS*(kall/sigmaSi)^2));
    % % % *** END don't consider Boundary scattering. 
    
    % % % *** Consider particle scattering, not Boundary Scattering. 
    tot_scat_inv3 = umk_scat_inv + imp_scat_inv + part_scat_inv;
    tot_scat3 = 1./tot_scat_inv3;
     kmode_holland3_Part = vmode.*vmode.*tot_scat3.*Cpmode; %
     ksumx_3Part(nk) = sum(sum(kmode_holland3_Part*(kall/sigmaSi)^2));
    % % % *** END Particle Scattering, not Boundary scattering. 

    
end
% % % % Cp = Cpsum*8/(N*N*N*(A1*sigmaSi)^3)/2329 %% % Why need to multiply 8??
% % % % % % % Subroutine: Two different Methods calculate thermal conductivity. 
% % % % % % % k = ksum*8/(N*N*N*(A1*sigmaSi)^3)
% % % % % % % k2 = ksum*dkx*dky*dkz/(2*pi)^3
% % % % % % % Subroutine: End 
% % % % k_density = ksum/(nxmax*nymax*nzmax*(A1*sigmaSi)^3)/2329
% % % % 
k = ksum*(dk)/(6*pi^2) % 8/(N*N*N*(A1*sigmaSi)^3); 
% k_110 = ksum*dk*sqrt(2)/(3*pi^2)
% Cp_110 = Cpsum*dk*sqrt(2)/(3*pi^2)/2329
Cp = Cpsum*dk/(3*pi^2)/2329; 
% k1_110 = ksum*dk*(4*pi/3)/(pi/sqrt(2))^3/8  %  %  % sqrt(2)/(3*pi^2)
% % % Cp1_110 = Cpsum*dk*sqrt(2)/(3*pi^2)/2329
k_NBS = sum(ksumx_NBS)*dk/(6*pi^2) 
k_3Part = sum(ksumx_3Part)*dk/(6*pi^2)
timevalue = toc
Nresult = [Cp, k,k_NBS,k_3Part, timevalue]
Tresult(ii,:) = Nresult; 
end
save('Nresult_110ncx6.mat','Nresult','-ascii')

% % figure
% % plot(kxx,wwxx,'k-');
% ylim([0 5])
% % % % %  Cpsum*(2*pi)^3/((A1*sigmaSi)^3)/3/2329
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

