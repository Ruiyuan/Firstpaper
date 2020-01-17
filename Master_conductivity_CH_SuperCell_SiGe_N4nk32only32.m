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

A = a0/sigmaSi; %  15168727030; 2.592215; % 2.592215168; %l
ncell=2;
nbasis=4;
naxis=3;
forceflag=4;  % fororceflag=1 : spring; forceflag=2: Silicon; forceflag=3: Germanium.  means choose s-w potential.

rc=1.8;

mI=1;
mTmatrix=[massGe/massSi]; % because of PhiGe has included the mass of massGe, so I just need reduced here as mT=1. 

%---------start--contstruct "force matrix" of pure Si or heavy Si-------
ncx = 4; % 2; % 5;
ncy = ncx; %5;
ncz = ncx; % 5;
% % % ncxh=3;
% % % ncyh=3;
% % % nczh=3;

ncxlarge = ncx; % 8; % Thi is for Multiple Voids. 
ncxsuper = ncx;
% % % ncysuper = ncy;
% % % nczsuper = ncz+nczh;
ncxvoid = 0; % 0; % 2;
% N = 50; 

% % % SWITCH BETWEEN PRIMITIVE CELL AND SUPERCELL!!!
        flagsupercell = 1;

% % % boxlx=ncx*A;
% % % boxly=ncy*A;
% % % boxlz=ncz*A;
% [xs ys zs]=creatdiamond(A,ncx,ncy,ncz,ncell,nbasis);   
% [xs ys zs]=creatSuperCell(A,ncx,ncy,ncz,ncell,nbasis,ncxsuper);   

% [xs,ys,zs,nsupercell]=creatSuperCellwVoid(A,ncx,ncy,ncz,ncell,nbasis,ncxsuper,ncxvoid,flagsupercell);   
% na=length(xs); % ncx*ncy*ncz*ncell*nbasis;
% figure;
% scatter3(xs,ys,zs);

if(ncxvoid>1)
    [xs,ys,zs,nsupercell]=creatSuperCellwMultipleVoid(A,ncx,ncy,ncz,ncell,nbasis,ncxsuper,ncxvoid,flagsupercell);   
else
    [xs,ys,zs,nsupercell]=creatSuperCellwVoid(A,ncx,ncy,ncz,ncell,nbasis,ncxsuper,ncxvoid,flagsupercell);   
end

na=length(xs); % ncx*ncy*ncz*ncell*nbasis;
% % % **** Add Quantum Dots here ****% 
ncxdot = 1; 
mI = ones(1,na); % % ones(1,nsupercell);
atype = ones(1,na); % % ones(1,nsupercell);
% % % % % % % ***** Si3x3x3 with Ge1x1x1 Super Cell. 
% % % % for i = (ncx*ncy + ncx + (ncx - ncxdot)/2)*8+1:1:(ncx*ncy + ncx + (ncx - ncxdot)/2 + 1)*8;
% % % %     mI(i) = massGe/massSi; % 1; % 
% % % %     atype(i) = 2; 
% % % % end
    mI1 = mI; 
    % % % *** To create SuperCell 8*(Si3x3x3wGe1x1x1)
        Nsc1 = 0;
        
        for i=0:ncxlarge/ncx-1; 
            for j=0:ncxlarge/ncy-1;
                for k=0:ncxlarge/ncz-1;            
                    xs((1:nsupercell)+Nsc1*nsupercell) = xs(1:nsupercell)+i*A*ncxsuper;
                    ys((1:nsupercell)+Nsc1*nsupercell) = ys(1:nsupercell)+j*A*ncxsuper;
                    zs((1:nsupercell)+Nsc1*nsupercell) = zs(1:nsupercell)+k*A*ncxsuper;
                    mI((1:nsupercell)+Nsc1*nsupercell) = mI(1:nsupercell);
                    atype((1:nsupercell)+Nsc1*nsupercell) = atype(1:nsupercell);  
                    
                    Nsc1 = Nsc1 + 1;
                  
                end
            end
        end
        na= length(xs);
        nsupercell = na; 
        if(ncxvoid == 0)
            L = 1e20;
        else
            L = (ncx/ncxvoid - ncxvoid)*a0;
        end
        
        ncxsuper = ncxlarge;
        ncx = ncxlarge; 
        ncy = ncxlarge;
        ncz = ncxlarge;
        % % % SWITCH BETWEEN PRIMITIVE CELL AND SUPERCELL!!!
        if(flagsupercell == 1)
            A1 = A*ncxsuper; % *ncx; % ; %
        else 
            A1 = A*1;
        end
        
% % % ***** Si2x2x2 with 2 Ge atoms in one CC (1st and 2nd atoms) *****
% % % % % % % for i = 1: 1: ncx*ncy*ncz; % (ncx*ncy + ncx + (ncx - ncxdot)/2)*8+1:1:(ncx*ncy + ncx + (ncx - ncxdot)/2 + 1)*8;
% % % % % % %     ni = (i-1)*8 + (1:2);
% % % % % % %     mI(ni) = massGe/massSi; % 1; % 
% % % % % % %     atype(ni) = 2; 
% % % % % % % end
% % % **** End of Add Quantum Dots

% Phipure=zeros(naxis,naxis,ncell,na);
% % % PhiSi=zeros(naxis,naxis,ncell,na);
% % % PhiGe=zeros(naxis,naxis,ncell,na);
atypepure=zeros(2,na);
atypepure(1,:)=ones(1,na);
atypepure(2,:)=2*ones(1,na);  % here "2" means Ge. 

% -----------READ MA wwmode-------
% % % load('C:\matlabcodes\Dispersion Curve\fmatSi.mat');
% % % load('C:\matlabcodes\Dispersion Curve\fmatGe.mat');
% % % load('C:\matlabcodes\Dispersion Curve\fmatSiSiGe.mat');
% % % load('C:\matlabcodes\Dispersion Curve\fmatSiGeGe.mat');
load('fmatSi.mat');
% % load('fmatGe.mat');
fmatSi = -fmatSi;
% % fmatGe = -fmatGe;
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

% nbranch = 3*ncell*nbasis*(ncxsuper^3 - ncxvoid); 
% for kx=ks:dk:kend;
%     kx=kxx(1); ky=kyy(1); kz=kzz(1);
%     kvI=[kx;ky;kz];  % on the [111] plane or direction. 
%     [wwI1,DI1,eV1]=dispersioncurveSuperCell(kvI,A,ncell,rc,Ksp,naxis,mI,nbasis,forceflag,fmatSi,xs,ys,zs,ncx,ncy,ncz,ncxsuper,nsupercell);  
% 
% parpool(4)
% kindex = nxmax;
% nbranch = 3*ncell*nbasis*(ncxsuper^3 - ncxvoid); 
if (sum(mI) - 1*length(mI) > 0)
    [dmatallD,dmatallD2,fmat0,D2,D3] = ComputeAllDmatSiGeQuantDot(xs,ys,zs,ncx,ncy,ncz,A,atype, eps, isigma,leps, gamma);
    fmatSiGeQD = -dmatallD2; 
    [Phiall, rxall, ryall, rzall]=dispersioncurveSuperCell_SiGe(A,ncell,rc,Ksp,naxis,mI,nbasis,forceflag,fmatSiGeQD,xs,ys,zs,ncx,ncy,ncz,ncxsuper,nsupercell);  

else 
%     fmatSiGeQD = fmatSi; 
[Phiall, rxall, ryall, rzall]=dispersioncurveSuperCell_Modify(A,ncell,rc,Ksp,naxis,mI,nbasis,forceflag,fmatSi,xs,ys,zs,ncx,ncy,ncz,ncxsuper,nsupercell);  

end

ninc = 0;
Cpsum = 0;
ksum = 0;

h = 1.054571726e-34; % reduced planck constant. 
kB = 1.3806488e-23;
% % % *** Davis 2011 Data
B = 2.1e-19; 
C = 180;
D = 1.32e-45;
n=0;


% % % % % % Hopkins 2010
% % % B = 3.73e-19; % unit: sK^-1;
% % % C = 157.3; % unit: K.
% % % D = 9.32e-45; % unit: s^3;
% % % E = 2.3e-3; % unit: m. 

% % % % % %  *** Mingo 2003, Calculation Silicon Nanowire:
% % % B = 1.73e-19; % unit: s/K.
% % % C = 137.3; % unit: K. 
% % % D = 1.32e-45; % s^3. % % % N = 30; 
Nk = [32]; %48; %  60; % [ 48]; % 80]; % 72]; % 36, 48 , 60]; % ]; % ,36, 48 , 60,  80, 96, 120, 160];    
Nresult = zeros(length(Nk),6);


for ink = 1:length(Nk)
    N = Nk(ink);

if(flagsupercell == 1)
    XX = ((0:1/N:1-1/N)+1/N/2)*1*pi/A1; %  (0.001:0.0625:1)*1*pi/(A1); 
    YY = ((0:1/N:1-1/N)+1/N/2)*1*pi/A1; % 0; %  %(-1:1/32:1)*1*pi/(A1); %(-1:0.0625:1)*1*pi/(A1); % 0.0*kxx;% 
    ZZ = ((0:1/N:1-1/N)+1/N/2)*1*pi/A1; %  0; % (-1:1/32:1)*1*pi/(A1); % (-1:0.0625:1)*1*pi/(A1); % 0.0*kxx;% (-1:0.0667:1)*1*pi/(A1); % 0.0*kxx; 
nbranch = nsupercell*3; % length(xs)*3; % ncell*nbasis*(ncxsuper^3*3 - ncxvoid);

else 
    XX = ((0:1/N:1-1/N)+1/N/2)*2*pi/A1; %
    YY = ((0:1/N:1-1/N)+1/N/2)*2*pi/A1; % 0; %  (-1:1/32:1)*1*pi/(A1); %(-1:0.0625:1)*1*pi/(A1); % 0.0*kxx;% 
    ZZ = ((0:1/N:1-1/N)+1/N/2)*2*pi/A1; % 0; %  (-1:1/32:1)*1*pi/(A1); % (-1:0.0625:1)*1*pi/(A1); % 0.0*kxx;% (-1:0.0667:1)*1*pi/(A1); % 0.0*kxx; 
nbranch = ncell*3; % ncell*nbasis*(ncxsuper^3*3 - ncxvoid);

end
% XX = ((0:1/N:1-1/N)+1/N/2)*2*pi/A1; %  (0.001:0.0625:1)*1*pi/(A1); 
% dk = 0.001*2*pi/A/sigmaSi; 
dkx = (XX(2) - XX(1))/sigmaSi; 
dky = (YY(2) - YY(1))/sigmaSi;  % 0; % 
dkz = (ZZ(2) - ZZ(1))/sigmaSi;  % 0; % 

Lx = 2*pi/(dkx); Ly = 2*pi/(dky); Lz = 2*pi/(dkz);
V = Lx*Ly*Lz;

% % % dkk = [0.001; 0; 0]*1*pi/(A1);
 wwindex=0;
%  matlabpool local 2;
nxmax=length(XX);
nymax=length(YY);
nzmax=length(ZZ);

% figure;
nk = 1;
for nx= 1:nxmax; % 1:25;  % nxmax;
    for ny= 1:nymax;
        for nz = 1:nzmax;
            if ((XX(nx) + ZZ(nz) + YY(ny))> (3/2*2*pi/A1))
               continue;
            else
               kxx(nk) = XX(nx);
               kyy(nk) = YY(ny);
               kzz(nk) = ZZ(nz);
               nk = nk+1;
            end 
        end
    end
end

 wwindex=0;
%  matlabpool local 2;
nxmax=length(kxx);
nymax=length(kyy);
nzmax=length(kzz);

kindex = nxmax;

% * * * END of k space grid

%  matlabmail('ruiyuanma@gmail.com','Alarm','M62 Code Start','townem62@gmail.com', '19580101')

% * * * k space grid




    T= 300%  5:5:400 % :4150;  
% % % for nz= 1:nzmax; % 1:25;  % nxmax;
% % %     n = n+1
% % %     for nz= 1:nzmax
% % %         wwmode = zeros(nxmax,nbranch);
% % %         for nx = 1:nxmax;
% % %     ky = kyy;
%     ky = kyy;
% % %     kz = kzz;
% matlabpool('open',2);
ww = zeros(length(kxx),nbranch);
vv = zeros(length(kxx),nbranch);
ksumx2 = zeros(nbranch,1);
ksumx_NBS = zeros(nbranch,1);
ksumx_3Part = zeros(nbranch,1);

fileID = fopen('Fig1_BulkSi_N4_nk_32_ksum_P4_32768.csv','w'); % 'ncx','N');

% parfor nk = 1:length(kxx);
for nk = 1+length(kxx)/4*3:length(kxx)/4*4;


%	printf ("Processed %d of '%d'.", nk, N^3);
%	fflush(stdout)
    if(mod(nk,10)==0)
    	disp(nk/length(kxx)) 
    end
% % %         wwmode = zeros(nbranch);
% % %         vmodex = zeros(nbranch);
% % %         vmodey = zeros(nbranch);
% % %         vmodez = zeros(nbranch);


    kx=kxx(nk);
    ky=kyy(nk);
    kz=kzz(nk);
    
    ninc = ninc + 1;
    Dsuper = zeros(naxis*nsupercell, naxis*nsupercell);
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
% % % Dall = Phiall.*exp(1i*(kxall.*rxall+kyall.*ryall+kzall.*rzall));
% clear kxall kyall kzall Dall;
Dall = Phiall.*exp(1i*(kx*rxall+ky*ryall+kz*rzall));

% % % Dsuper = Dall;
for i = 1:na/nsupercell;
Dsuper = Dsuper + Dall(:,(i-1)*naxis*nsupercell+(1:naxis*nsupercell));
end
% % %  [UU,SS,VV]=svdecon(Dsuper);
% % %   [UU, SS]=eig(Dsuper);
% % %   eV=UU ; % 
% % %   wwI = sqrt(diag(SS));
% % %   [eV,SS]=eig(Dsuper);
  [eV,SS]=svd(Dsuper);
  wwI = sqrt(diag(SS));

% % %  eV=zeros(naxis*nsupercell,naxis*nsupercell);

% % %  for i=1:naxis*nsupercell;
% % %     om(i,1)=SS(i,i);
% % %  end
% % %  wwI=(sqrt(om));   
% % %  wwxx2(nx,nz,nz,:)=wwI*sqrt(1/massSi); %
 
 % % %### Dispersion Values at a smaller forward. 
% % %  kxallf = (kx +0.001*1*pi/(A1))*ones(naxis*nsupercell,naxis*na);
% % %  kxallb = (kx -0.001*1*pi/(A1))*ones(naxis*nsupercell,naxis*na);
% %     clear Dall UU SS VV; 
    Dsuper = zeros(naxis*nsupercell, naxis*nsupercell);
    kxf = (kx+0.05*kx); % (kx + 0.025*1*pi/(A1)); %(kx+0.01); %  
    Dall = Phiall.*exp(1i*(kxf*rxall+ky*ryall+kz*rzall)); 
  
    for i = 1:na/nsupercell;
        Dsuper = Dsuper + Dall(:,(i-1)*naxis*nsupercell+(1:naxis*nsupercell));
% % %         Dsuperb = Dsuperb + Dallb(:,(i-1)*naxis*nsupercell+(1:naxis*nsupercell));
    end
% % %     [UU,SS,VV]=svdecon(Dsuper);
  [eVf,SS]=svd(Dsuper);
%   [eVf,SS]=eig(Dsuper);

  % % % %  eVf=zeros(naxis*nsupercell,naxis*nsupercell);
% %      eVf=UU ; % 
    wwf=(sqrt(diag(SS))); 
% % % % % % % % % ***Backward Step. **********************    
% % % % % % % %     clear Dall UU SS VV; 
% % % % % %     Dsuper = zeros(naxis*nsupercell, naxis*nsupercell);    
% % % % % %     kxb = (kx -0.001*1*pi/(A1));
% % % % % % 
% % % % % %     Dall = Phiall.*exp(1i*(kxb*rxall+ky*ryall+kz*rzall));
% % % % % % % % %  Dallf = Phiall.*exp(1i*(kxallf.*rxall+kyall.*ryall+kzall.*rzall));
% % % % % % % % %  Dallb = Phiall.*exp(1i*(kxallb.*rxall+kyall.*ryall+kzall.*rzall));
% % % % % % % % % %  Dsuperf = Dallf;
% % % % % %     for i = 1:na/nsupercell;
% % % % % % %         Dsuperf = Dsuperf + Dallf(:,(i-1)*naxis*nsupercell+(1:naxis*nsupercell));
% % % % % %         Dsuper = Dsuper + Dall(:,(i-1)*naxis*nsupercell+(1:naxis*nsupercell));
% % % % % %     end 
% % % % % % % % %   [UU,SS,VV]=svdecon(Dsuper);
% % % % % %   [UU,SS]=eig(Dsuper);
% % % % % % % % % %  eVf=zeros(naxis*nsupercell,naxis*nsupercell);
% % % % % %  eVb=UU ; % 
% % % % % %   wwb=(sqrt(diag(SS))); 
% % % % % % % % % ***Backward Step. END**********************    

% % %     clear Dall Dsuper  UU SS VV;  
    Dall = []; Dsuper = []; UU = []; SS = []; VV = [];

% % % % % % % % % % ### GROUP VELOCITY ***BACKWARD***######## 
% % % % %     checkeV = abs((eVb'*eVf));% ./(norm(eVb).*norm(eVf')));
% % % % %     [M, ip] = max(checkeV');
% % % % % %             clear checkeV M;
% % % % % %     checkeV = []; M = []; 
% % % % %     vmode = (wwb - wwf(ip))*sqrt(1/massSi)/(2*0.001*1*pi/A1/sigmaSi);
% % % % % % % % % % ### GROUP VELOCITY ***BACKWARD***######## 

% % % % % ### GROUP VELOCITY ***FORWARD***######## 
    checkeV = abs((eV'*eVf));% ./(norm(eVb).*norm(eVf')));
    [M, ip] = max(checkeV');
%             clear checkeV M;
%     checkeV = []; M = []; 
    vmode = (wwI - wwf(ip))*sqrt(1/massSi)/(1*0.05*kx/sigmaSi); % (1*0.025*1*pi/A1/sigmaSi); % (1*0.01/sigmaSi); % 
% % % % % ### GROUP VELOCITY ***RORWARD***######## 

% % % ### Array-wise Operation: The END.
% % % % 
% % % %     wwmode(nx,:) = wwI*sqrt(1/massSi); 
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
        Cpsum = Cpsum + sum(sum((Cpmode))); % kB*(h*wwmode/kB/T).^2.*(exp(h*wwmode/kB/T)./((exp(h*wwmode/kB/T)-1).^2))));
        
        
            umk_scat_inv = B*(wwmode.*wwmode)*T*exp(-C/T);
    % impurity scattering. 
    imp_scat_inv = D*wwmode.^4;
    % boundary scattering.
    bound_scat_inv = abs(vmode)/L; % E; %  (vmode)/E; %0; %  (1-p)/(1+p)*(1/dx + 1/dy)*abs(vmode)/2;
    % partical scattering. from Majumdar-1993
    rvoid = ncxvoid*a0/2;
    chi = rvoid*sqrt(kx^2+ky^2+kz^2)/sigmaSi;
    sig_void = pi*rvoid^2*(chi^4/(chi^4+1));
    rho_void = (ncxvoid/ncx)^3;
    part_scat_inv = sig_void*rho_void*abs(vmode);
    
    tot_scat_inv = umk_scat_inv + imp_scat_inv  + bound_scat_inv; % ;%
    tot_scat = 1./tot_scat_inv;
    kmode_holland = vmode.*vmode.*tot_scat.*Cpmode; %  + wy(nk)*vmodey.*vmodey.*tot_scat.*Cpmode + wz(nk)*vmodez.*vmodez.*tot_scat.*Cpmode; %
    

% % %     ksum = ksum + sum(sum(kmode_holland));
% % %     ksumx(nk) = ksum;
% % %     ksum2 = sum(sum(kmode_holland));
% % %     ksumx2(nk) = ksum2;
    ksumx2(nk) = sum(sum(kmode_holland))*8/V;
    % % % *** Don't consider Boundary Scattering. ***
    tot_scat_inv2 = umk_scat_inv + imp_scat_inv;
    tot_scat2 = 1./tot_scat_inv2;
     kmode_holland2_NBS = vmode.*vmode.*tot_scat2.*Cpmode; %
     ksumx_NBS(nk) = sum(sum(kmode_holland2_NBS))*8/V;
    % % % *** END don't consider Boundary scattering. 
    
    % % % *** Consider particle scattering, not Boundary Scattering. 
    tot_scat_inv3 = umk_scat_inv + imp_scat_inv + part_scat_inv;
    tot_scat3 = 1./tot_scat_inv3;
     kmode_holland3_Part = vmode.*vmode.*tot_scat3.*Cpmode; %
     ksumx_3Part(nk) = sum(sum(kmode_holland3_Part))*8/V;
    % % % *** END Particle Scattering, not Boundary scattering. 
    
    ww(nk,:) = wwmode;
    vv(nk,:) = vmode;
% % %     clear vmode wwmode imp_scat_inv bound_scat_inv 
%     vmode = []; wwmode = []; imp_scat_inv = []; bound_scat_inv = []; 
%     fprintf(fileID,'%6s %12s\n','x','exp(x)');
% if(mod(nk,10)==0)
    fprintf(fileID,'%6d %12.8f %12.8f %12.8f\n', nk, ksumx2(nk), ksumx_NBS(nk), ksumx_3Part(nk));
% end

end
% matlabpool('close');
    fclose(fileID);

ksum = sum(sum(ksumx2));
% % % Subroutine: Calculate thermal conductivity.
% Method #1! Be careful about the Volume for Primitive and Super Cell! 
% % % Alos be careful are 8 time difference for Primitive and Conventional
% Cell. 
k = ksum*1/(N*N*N*(A1*sigmaSi)^3) %  Super Cell/Conventional Cell Only.
k2 = 8*ksum*dkx*dky*dkz/(2*pi)^3 % Super Cell/Conventional Cell Only.
Cp = Cpsum/(N*N*N*(A1*sigmaSi)^3)/2329 %  Super Cell/Conventional Cell Only. 

% % % k2 = ksum*dkx*dky*dkz/(2*pi)^3 %  % Primitive Cell Only  
% % % Cp = Cpsum*8/(N*N*N*(A1*sigmaSi)^3)/2329 % Primitive Cell Only
% % % Following Equations are capable for Both Super Cell and Primitive
% Cell. 
Lx = 2*pi/(dkx); Ly = 2*pi/(dky); Lz = 2*pi/(dkz);
V = Lx*Ly*Lz;
k3 = sum(sum(ksumx2)); %(8*ksum)/V 
Cp3 = Cpsum*8/V/2329 % 

% %% 
k_NBS = sum(sum(ksumx_NBS)); % sum(sum(ksumx_NBS))*8/V
k_Particle = sum(sum(ksumx_3Part)); % sum(sum(ksumx_3Part))*8/V

% % % Subroutine: End
% % % save('k_N8.mat','k','k2','k3','Cp');
timevalue = toc

Nresult(ink,:) = [Nk(ink),Cp3,k_NBS,k3,k_Particle,timevalue];
Cp3 = 0; Cpsum = 0; tic

end

xlswrite('Fig1_BulkSi_N4_nk_32_P4_32768.xlsx',Nresult)
% xlswrite('Fig1_BulkSi_N4_nk_24.xlsx',Nresult)
% 
% % % xlswrite('kxx-ncx6-ncxvoid2-Methok3d1-N5.xlsx',kxx)
% % % xlswrite('ksumx-ncx6-ncxvoid2-Method1-N5.xlsx',ksumx)

% % % % csvwrite('ksumx2-Si6x6x6wGe1x1x1-8-Method1-N15.xlsx',ksumx2)
% % % % csvwrite('ww-Si6x6x6wGe1x1x1-8-Method1-N15.xlsx',ww)
% % % % csvwrite('vv-Si6x6x6wGe1x1x1-8-Method1-N15.xlsx',vv)

% figure
% plot(kxx,wwxx,'k);
% ylim([0 5])

% % % %  Cpsum*(2*pi)^3/((A1*sigmaSi)^3)/3/2329
% % % %  figure;
% % % %  plot(kxx,ww,'bo-');
% % % %     xlabel('kx')
% % % %     ylabel('Angular Frequency (rad/s)')
% % % %     legend('Ge1x1x1 in Si3x3x3')
% % % %     set(gca,'Linewidth',3.0)
% % % %  figure
% % % %  plot(kxx,vv,'bo');
% % % %     xlabel('kx')
% % % %     ylabel('Group velocity (m/s)')
% % % %     legend('Ge1x1x1 in Si3x3x3')
% % % %     set(gca,'Linewidth',3.0)
    
    
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
% % plot(kxx,wwmode,'k-')
% % figure
% % plot(kxx,vmode,'b.')
% % hold on;
% % plot(kxx,vg,'ro')
% % plot(kxx,vg2,'go')

plot(N2,K2,'--ks','LineWidth',2,'MarkerSize',8)
