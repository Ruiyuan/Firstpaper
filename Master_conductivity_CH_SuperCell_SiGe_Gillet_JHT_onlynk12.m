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
ncx = 3; % 2; % 5;
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
ncydot = ncxdot;
nczdot = ncxdot;

mI = ones(1,na); % % ones(1,nsupercell);
atype = ones(1,na); % % ones(1,nsupercell);

% % % ***** Si5x5x5 with embedded Ge3x3x3 Super Cell.
for i=1:length(xs);
    xdot = (xs(i)<((ncx+ncxdot)/2*A - 0.01)) && (xs(i)>=((ncx-ncxdot)/2*A-0.01));
    ydot = (ys(i)<((ncy+ncydot)/2*A - 0.01)) && (ys(i)>=((ncy-ncydot)/2*A-0.01));
    zdot = (zs(i)<((ncz+nczdot)/2*A - 0.01)) && (zs(i)>=((ncz-nczdot)/2*A-0.01));
    if(xdot && ydot && zdot)
         mI(i) = massGe/massSi; % 1; % 
         atype(i) = 2;
    end
end

% % % % % % % % ***** Si3x3x3 with Ge1x1x1 Super Cell. 
% % % % % for i = (ncx*ncy + ncx + (ncx - ncxdot)/2)*8+1:1:(ncx*ncy + ncx + (ncx - ncxdot)/2 + 1)*8;
% % % % %     mI(i) = massGe/massSi; % 1; % 
% % % % %     atype(i) = 2; 
% % % % % end
% % % % %     mI1 = mI; 
    
    
    
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
% % % % atypepure=zeros(2,na);
% % % % atypepure(1,:)=ones(1,na);
% % % % atypepure(2,:)=2*ones(1,na);  % here "2" means Ge. 

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
Cpsum100 = 0;
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
Nk = [10];    
Nresult = zeros(length(Nk),7);


for ink = 1:length(Nk)
    N = Nk(ink)

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

% % kxx = XX;

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
ksumx_4Gillet = zeros(nbranch,1);
ksum = 0;
% parfor nk = 1:length(kxx);

% fileID = fopen('Fig3_SiGeQD_Gillet_N5_nk4_ksum.csv','w'); % 'ncx','N');

for nk = 1:length(kxx);


%	printf ("Processed %d of '%d'.", nk, N^3);
%	fflush(stdout)
    if(mod(nk,10)==0)
    	disp(nk/length(kxx)) 
    end
% % %         wwmode = zeros(nbranch);
% % %         x = zeros(nbranch);
% % %         y = zeros(nbranch);
% % %         z = zeros(nbranch);


    kx=kxx(nk);
    ky=kyy(nk); % 0; % 
    kz=kzz(nk); % 0; % 
    kall = sqrt(kx^2 + ky^2 + kz^2);
    dk = dkx*(kall/kx);
    
    ninc = ninc + 1;
    
    Dsuper = zeros(naxis*nsupercell, naxis*nsupercell);
Dall = Phiall.*exp(1i*(kx*rxall+ky*ryall+kz*rzall));
for i = 1:na/nsupercell;
Dsuper = Dsuper + Dall(:,(i-1)*naxis*nsupercell+(1:naxis*nsupercell));
end

  [eV,SS]=svd(Dsuper);
  wwI = sqrt(diag(SS));
  ww(nk,:) = wwI;
  
 dkxf = 0.001*1*pi/A1; % 
 % % %### Dispersion Values at a smaller forward. 
    
    Dsuper = zeros(naxis*nsupercell, naxis*nsupercell);
    kxf = (kx+dkxf); % (kx + 0.025*1*pi/(A1)); %(kx+0.01); %  
    kyf = ky; % kxf; 
    Dall = Phiall.*exp(1i*(kxf*rxall+kyf*ryall+kz*rzall)); 
  
    for i = 1:na/nsupercell;
        Dsuper = Dsuper + Dall(:,(i-1)*naxis*nsupercell+(1:naxis*nsupercell));
    end
  [eVf,SS]=svd(Dsuper);
    wwf=(sqrt(diag(SS))); 
    
% % % % % ***Backward Step. **********************    
% %     Dsuper = zeros(naxis*nsupercell, naxis*nsupercell);    
% %     kxb = (kx - dkxf);
% %     kyb = ky; % kxb; 
% %     
% %     Dall = Phiall.*exp(1i*(kxb*rxall+kyb*ryall+kz*rzall));
% %     for i = 1:na/nsupercell;
% %         Dsuper = Dsuper + Dall(:,(i-1)*naxis*nsupercell+(1:naxis*nsupercell));
% %     end 
% %     [UU,SS]=eig(Dsuper);
% %     eVb=UU ; % 
% %     wwb=(sqrt(diag(SS))); 
% % % % % ***Backward Step. END**********************    

    Dall = []; Dsuper = []; UU = []; SS = []; VV = [];

% % % % % % % ### GROUP VELOCITY ***Central***######## 
% %     checkeV = abs((eVb'*eVf));
% %     [M, ip] = max(checkeV');
% %     dkxvg = (kxf - kxb);
% %     vmodec = (wwb - wwf(ip))*sqrt(1/massSi)/(1*dkxvg/sigmaSi);
% %     wwc = (wwb + wwf(ip))/2; % wwI; 
% %     wwI = wwc; 
% %     vmode = vmodec;
% % % % % % % ### GROUP VELOCITY ***Central***######## 

% % % % % % % ### GROUP VELOCITY ***BACKWARD***######## 
% %     checkeV = abs((eVb'*eV));% ./(norm(eVb).*norm(eVf')));
% %     [M, ipb] = max(checkeV');
% %         dkxvg = (kx - kxb);
% %     vmodeb = (wwb - wwI(ipb))*sqrt(1/massSi)/(1*dkxvg/sigmaSi); % (1*0.001*kx/sigmaSi); % (1*0.025*1*pi/A1/sigmaSi); % (1*0.01/sigmaSi); % 
% %     wwI = wwI(ipb);
% %     vmode = vmodeb;
% %     % % % % % ### GROUP VELOCITY ***BACKWARD***######## 


% % % % % ### GROUP VELOCITY ***FORWARD***######## 
    checkeV = abs((eV'*eVf));% ./(norm(eVb).*norm(eVf')));
    [M, ipf] = max(checkeV');
%             clear checkeV M;
%     checkeV = []; M = []; 
    dkxvg = (kxf - kx);
    vmodef = (wwI - wwf(ipf))*sqrt(1/massSi)/(1*dkxvg/sigmaSi); % (1*0.001*kx/sigmaSi); % (1*0.025*1*pi/A1/sigmaSi); % (1*0.01/sigmaSi); % 
% % % % % ### GROUP VELOCITY ***RORWARD***######## 
    vmode = vmodef;


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
    

    ksum = ksum + sum(sum(kmode_holland*(kall/sigmaSi)^2))*8/V;
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
    
    % % % *** Gillet: incoherent scattering and Umklapp scattering. 
   
   %  Uprocess_scat_inv = bu*(wwmode.*wwmode)./(vmode.*vmode);
    gamma = 1.5; theta = 645; % Debye temperature, unit: K
    Uprocess_scat_inv = (h*gamma^2*T/massSi/theta*exp(-theta/3/T))*(wwmode.*wwmode)./(vmode.*vmode);
    dA = 1.5849; dK = 0.24;  % % % data from Gillet-JHT-2009. 
    [inc_scat_inv] = incoherent_scat(ncx,ncxdot,dA,dK,kx,ky,kz,a0,vmode,sigmaSi);
    
    tot_scat_inv_Gillet = Uprocess_scat_inv + inc_scat_inv;
    tot_scat_Gillet = 1./tot_scat_inv_Gillet;
    kmode_holland4_Gillet = vmode.*vmode.*tot_scat_Gillet.*Cpmode; %
    ksumx_4Gillet(nk) = sum(sum(kmode_holland4_Gillet))*8/V;
% % %     ksumx_4Gillet_100(nk) = sum(sum(kmode_holland4_Gillet*(kall/sigmaSi)^2));
% % %     Cpsum100 = Cpsum100 + sum(sum((Cpmode*(kall/sigmaSi)^2)));
    % % % ***Gillet: END
    
    ww(nk,:) = wwmode;
    vv(nk,:) = vmode;
%     ipip(nk,:) = ip;
% % %     clear vmode wwmode imp_scat_inv bound_scat_inv 
%     vmode = []; wwmode = []; imp_scat_inv = []; bound_scat_inv = []; 
%       fprintf(fileID,'%6d %12.8f %12.8f %12.8f %12.8f\n', nk, ksumx2(nk), ksumx_NBS(nk), ksumx_3Part(nk),ksumx_4Gillet(nk));

end
%     fclose(fileID);

% matlabpool('close');
% % % ******* Plot dispersion curve. 
% % % % figure
% % % % for i=[nbranch:-1:nbranch-50,nbranch-51:-30:1]; % 39,38:-1:1]; % 1:8:nbranch; % nbranch
% % % %     plot(XX*ncx*A/pi,ww(:,i)*sqrt(1/massSi)/(2*pi)/1e12,'b-')
% % % %     hold on;
% % % % end
% % % % ylabel('\omega/(2\pi) (THz)')
% % % % xlabel('kd/\pi')
% % % ******* Plot dispersion curve END. 

ksum = sum(sum(ksumx2));
% % % Subroutine: Calculate thermal conductivity.
% Method #1! Be careful about the Volume for Primitive and Super Cell! 
% % % Alos be careful are 8 time difference for Primitive andConventional
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
k3 = ksum; % (8*ksum)/V 
Cp3 = Cpsum*8/V/2329; % 

% %% 
k_NBS = sum(sum(ksumx_NBS)); % sum(sum(ksumx_NBS))*8/V
k_Particle = sum(sum(ksumx_3Part)); % sum(sum(ksumx_3Part))*8/V

k_Gillet = sum(sum(ksumx_4Gillet)); % sum(sum(ksumx_4Gillet))*8/V
% % % Subroutine: End
% % % save('k_N8.mat','k','k2','k3','Cp');
timevalue = toc

Nresult(ink,:) = [Nk(ink),Cp3,k_NBS,k3,k_Particle,k_Gillet, timevalue];
Cp3 = 0; Cpsum = 0; tic

% % % % % % % % ---- Direction: 100 Results!----
% % % % % k = ksum*(dk)/(6*pi^2) % 8/(N*N*N*(A1*sigmaSi)^3); 
% % % % % k_100 = ksum*dk/(3*pi^2)
% % % % % Cp_100 = Cpsum100*dk/(3*pi^2)/2329
% % % % % 
% % % % % k1_100 = ksum*dk*(4*pi/3)/(pi/sqrt(2))^3/8  %  %  % sqrt(2)/(3*pi^2)
% % % % % 
% % % % % kGillet_100 = sum(sum(ksumx_4Gillet_100))*dk/(3*pi^2)
% % % % % 
% % % % % timevalue = toc
% % % % % Nresult(ink,:) = [Nk(ink),Cp_100,k,k_100,k1_100,kGillet_100,timevalue];
% % % % % % % % ----

end
% xlswrite('Fig3_SiGeQD_Gillet_nk8_FullyBZ.xlsx',Nresult)

% % % xlswrite('kxx-ncx6-ncxvoid2-Methok3d1-N5.xlsx',kxx)
% % % xlswrite('ksumx-ncx6-ncxvoid2-Method1-N5.xlsx',ksumx)

% % % % csvwrite('ksumx2-Si6x6x6wGe1x1x1-8-Method1-N15.xlsx',ksumx2)
% % % % csvwrite('ww-Si6x6x6wGe1x1x1-8-Method1-N15.xlsx',ww)
% % % % csvwrite('vv-Si6x6x6wGe1x1x1-8-Method1-N15.xlsx',vv)

% figure
% plot(kxx,wwxx,'k-');
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

