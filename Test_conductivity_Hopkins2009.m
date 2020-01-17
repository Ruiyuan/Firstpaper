% % % # # # DISPERSION CURVE CODE CAN SWITCH BETWEEN PRIMITIVE CELL AND

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
ncx =  2;
ncy = ncx; %5;
ncz = ncy; % 5;

ncxsuper = ncx; % ncx; % ncx; %  ncx; % ncxsuper must >= 2 for Primitive Cell! 
ncxvoid = 0; % 1;

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

[xs,ys,zs,nsupercell]=creatSuperCellwVoid(A,ncx,ncy,ncz,ncell,nbasis,ncxsuper,ncxvoid,flagsupercell);   
na=length(xs); % ncx*ncy*ncz*ncell*nbasis;
% % % figure;
% % % scatter3(xs,ys,zs);


atypepure=zeros(2,na);
atypepure(1,:)=ones(1,na);
atypepure(2,:)=2*ones(1,na);  % here "2" means Ge. 

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

 wwindex=0;
kxx = XX;
nxmax=length(kxx);
% % % nymax=length(kyy);
% % % nzmax=length(kzz);

kindex = nxmax;
nbranch = 1*3; % length(xs)*3; % ncell*nbasis*(ncxsuper^3*3 - ncxvoid);
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
% % % Davis
% % B = 2.1e-19; 
% % C = 180;
% % D = 1.32e-45;
% % E = 1e20; 
% % n=0;
% % % Hopkins 2010
B = 3.73e-19; % unit: sK^-1;
C = 157.3; % unit: K.
D = 9.32e-45; % unit: s^3;
E = 2.3e-3; % unit: m. 

%  matlabmail('ruiyuanma@gmail.com','Alarm','M62 Code Start','townem62@gmail.com', '19580101')

    T = 300

if(ncxvoid == 0)
            L = 1e20;
        else
            L = (ncx/ncxvoid - ncxvoid)*a0;
end
TT = 40:10:400
for ii = 1:length(TT); %   5:5:400 % :4150;  
    T = TT(ii);
    ninc = 0
Cpsum = 0;
ksum = 0;

 ksumx_NBS = zeros(length(XX),1);
ksumx_3Part = zeros(length(XX),1);
ww = zeros(length(XX),nbranch);    
A4 = 1.37E-27; % ?m4?s?1; % 
A3 = -3.53E-17; % ?m3?s?1
A2 = 2.94E-8; % ?m2?s?1; % 
A0 = 8350; % ?m?s?1; % 

% % % B4 = 1.94E-27; % ?m4?s?1
% % % B3 = -3.36E-17;% m3?s?1
% % % B2 = -1.86E-7; % ?m2?s?1
% % % B1 = 6090; % ?m?s?1
B4 = 1.57684126475862e-27; % 1.94E-27; % ?m4?s?1
B3 = -2.21549587094015e-17; % -3.36E-17;% m3?s?1
B2 = -2.96043290448924e-07; % 1.86E-7; % ?m2?s?1
B1 = 6408.76093302284; % ?m?s?1
B0 = -43604080535.7293; % 
TT = 40:10:400; %

for nk = 1:length(XX);
        data = zeros(nbranch);
        vmodex = zeros(nbranch);
        vmodey = zeros(nbranch);
        vmodez = zeros(nbranch);
    if(mod(nk,10) == 0)
        disp(nk/length(XX));
    end
    kx=kxx(nk)/sigmaSi;
    ky= 0; % kx; %  kyy(nk);
    kz= 0; % kzz(nk);
    kall = sqrt(kx^2 + ky^2 + kz^2);
    dk = dkx*(kall/kx);
    ninc = ninc + 1;
   

% % %     Dsuper = zeros(naxis*nsupercell, naxis*nsupercell);
% % % 
% % % Dall = Phiall.*exp(1i*(kx*rxall+ky*ryall+kz*rzall));
% % % 
% % % % % % Dsuper = Dall;
% % % for i = 1:na/nsupercell;
% % % Dsuper = Dsuper + Dall(:,(i-1)*naxis*nsupercell+(1:naxis*nsupercell));
% % % end

%   [eV,SS]=svd(Dsuper);
%   wwI = sqrt(diag(SS));
    kx = kx; % 
    wwI =[A4*kx^4+A3*kx^3+A2*kx^2+A0*kx; B4*kx^4+B3*kx^3+B2*kx^2+B1*kx+B0; B4*kx^4+B3*kx^3+B2*kx^2+B1*kx+B0]; % 

   dkxf = 0.001*kx; % 1*pi/A1/sigmaSi; % 
   kxf = kx + dkxf; 
 % % %### Dispersion Values at a smaller forward. 
    
    wwf= [A4*kxf^4+A3*kxf^3+A2*kxf^2+A0*kxf; B4*kxf^4+B3*kxf^3+B2*kxf^2+B1*kxf; B4*kxf^4+B3*kxf^3+B2*kxf^2+B1*kxf]; %
    kxb = (kx - dkxf);
    wwb= [A4*kxb^4+A3*kxb^3+A2*kxb^2+A0*kxb; B4*kxb^4+B3*kxb^3+B2*kxb^2+B1*kxb; B4*kxb^4+B3*kxb^3+B2*kxb^2+B1*kxb]; %
    vmode = (wwf - wwb)/dkxf/2; % *sqrt(1/massSi)/(1*dkxvg/sigmaSi); %

      wwmode = wwI; % ./(2*pi);; %  % *sqrt(1/massSi); 
%         end
        Cpmode = kB*(h*wwmode/kB/T).^2.*(exp(h*wwmode/kB/T)./((exp(h*wwmode/kB/T)-1).^2));
        Cpsum = Cpsum + sum(sum((Cpmode*(kall)^2))); % kB*(h*data/kB/T).^2.*(exp(h*data/kB/T)./((exp(h*data/kB/T)-1).^2))));
        
            umk_scat_inv = B*(wwmode.*wwmode)*T*exp(-C/T);
    % impurity scattering. 
    imp_scat_inv = D*wwmode.^4;
    % boundary scattering.
    bound_scat_inv = (vmode)/E; %L; %  E; %  (vmode)/E; %0; %  (1-p)/(1+p)*(1/dx + 1/dy)*abs(vmode)/2;
 % partical scattering. from Majumdar-1993
 
    
    tot_scat_inv = umk_scat_inv + imp_scat_inv  + bound_scat_inv; % ;%
    tot_scat = 1./tot_scat_inv;
    kmode_holland = vmode.*vmode.*tot_scat.*Cpmode; %  + wy(nk)*vmodey.*vmodey.*tot_scat.*Cpmode + wz(nk)*vmodez.*vmodez.*tot_scat.*Cpmode; %
    ksum = ksum + sum(sum(kmode_holland.*kall^2));
    ww(nk,:)=wwI;
    
end

k = ksum*(dk)/(6*pi) 
Tresult(ii,:) = [T,k];
end
