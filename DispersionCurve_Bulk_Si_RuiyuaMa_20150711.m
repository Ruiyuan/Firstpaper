clear;
clc;


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
PhiSi=zeros(naxis,naxis,ncell,na);
PhiGe=zeros(naxis,naxis,ncell,na);
atypepure=zeros(2,na);
atypepure(1,:)=ones(1,na);
atypepure(2,:)=2*ones(1,na);  % here "2" means Ge. 

% -----------READ MA DATA-------
load('fmatSi.mat');
load('fmatGe.mat');
% % load('C:\matlabcodes\Dispersion Curve\fmatSiSiGe.mat');
% % load('C:\matlabcodes\Dispersion Curve\fmatSiGeGe.mat');
fmatSi = -fmatSi;
fmatGe = -fmatGe;
% % fmatSiSiGe = -fmatSiSiGe;
% % fmatSiGeGe = -fmatSiGeGe;

% %----------READ ZHAO DATA
%  ASi=dlmread('fmat.Si');
% AGe=dlmread('fmat.Ge');
% ASiSiGe=dlmread('fmat.SiSiGe');
% ASiGeGe=dlmread('fmat.SiGeGe');
% AGeGeSi=dlmread('fmat.GeGeSi');
% AGeSiSi=dlmread('fmat.GeSiSi');
% 
% % global fmatSi fmatGe fmatSiSiGe fmatSiGeGe fmatGeGeSi fmatGeSiSi a0 massSi massGe iota;
% 
% n=0;
% for i=-1:1
% for j=-1:1;
% for k= -1:1;
% if(mod(i+j+k,2)~=0)
% continue;
% end
% n=n+1;
% fmatSi(:,:,i+2,j+2,k+2)=-ASi((n-1)*7+2:n*7,:);
% fmatGe(:,:,i+2,j+2,k+2)=-AGe((n-1)*7+2:n*7,:);
% fmatSiSiGe(:,:,i+2,j+2,k+2)=-ASiSiGe((n-1)*7+2:n*7,:);
% fmatSiGeGe(:,:,i+2,j+2,k+2)=-ASiGeGe((n-1)*7+2:n*7,:);
% fmatGeGeSi(:,:,i+2,j+2,k+2)=-AGeGeSi((n-1)*7+2:n*7,:);
% fmatGeSiSi(:,:,i+2,j+2,k+2)=-AGeSiSi((n-1)*7+2:n*7,:);
% 
% end
% end
% end
% %-------END READ ZHAO DATA--------

% for nmT=1:length(mTmatrix);%-0.3:0.4;
nmT=1; 
    clear transs reff checkk; 
    mT=mTmatrix(nmT); 
    nT=1;
    

%----------positions of bulk material atoms. ---------------%

%--over-----positions of bulk material atoms. -----over----%
%----------positions of interface atoms. ---------------%
   
%----over ---positions of interface atoms. ------over--%

Ksp=1.0;


nmode=naxis*ncell;
wwindex=0;

clear kx;
kxx = (0.001:0.001:1)*2*pi/A; % a0; %  0.01:dk: 2*pi/A; % 
dk = 0.001*2*pi/A/sigmaSi; 
 wwindex=0;
%  matlabpool local 2;
nxmax=length(kxx);
% for kx=ks:dk:kend;
for nx=1:nxmax;
    kx=kxx(nx);
    wwindex=wwindex+1;
    ninc=0;
    kxo(nx,:)=kx;
    ky= 1.0*kx;
    kz= 1.0*kx;  
    clear wwI eVoutI DI kvI; 
    
    kvI=[kx;ky;kz];  % on the [111] plane or direction. 
%     [DLi,D0i,DRi]=DISPERSIONMATRIX(mI,A,kvI,naxis,ncell,nbasis,fmatSi,xs,ys,zs);
    [wwI eVoutII DI]=dispersioncurve3ddiamondspring(kvI,A,ncell,rc,Ksp,naxis,mI,nbasis,forceflag,fmatSi,xs,ys,zs);  
    %     wwI = roundn(wwI, -14);
    wwxx(nx,:)=wwI*sqrt(1/massSi); % /sigmaSi; % anular frequency of Silicon. 
end
kall = kxx/sigmaSi;
wwmode = wwxx; % Angular frequency of Silicon. 
freq = wwxx/(2*pi); % Freqeucny of silicon dispersion curve. 
figure
plot(kall, wwmode);

%##################Group Velocity Calculation. ##################
%### Method 1 ###
kindex = nxmax; % length(kall);
ncell = 2; 
nbranch = 3*ncell;

ww_move_top=zeros(kindex,nbranch);
ww_move_bot=zeros(kindex,nbranch);
dW=zeros(kindex,nbranch);
vmode=zeros(kindex,nbranch);
ww_move_top(1:kindex-1,:)=wwmode(2:kindex,:);
ww_move_bot(2:kindex,:)=wwmode(1:kindex-1,:);
dW(1,:) = ww_move_top(1,:) - wwmode(1,:);
dW(kindex,:) = wwmode(kindex,:) - ww_move_bot(kindex,:); % ww_move_top(kindex,:); %
dW(2:kindex-1,:) = ww_move_top(2:kindex-1,:) - ww_move_bot(2:kindex-1,:); %  - 2*wwmode(2:kindex-1,:);
vmode(1,:)=dW(1,:)/dk;
vmode(2:kindex-1,:)=dW(2:kindex-1,:)/(dk*2);
vmode(kindex,:)=dW(kindex,:)/dk;
figure
plot(kall, vmode)

%##########################################
%## USE DATA from P. E. Hopkin, APL 95, 161902 (2009).
%## not Mingo or Nika data
%##########################################

% dx = 1.08e-09; dy=1.08e-09;
h = 1.054571726e-34;
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
for T= 300% :5:400 % :4150;
    n=n+1;
p=0; % specular parameter of the boundary roughness.
    % umklapp scattering. 
    umk_scat_inv = B*(wwmode.*wwmode)*T*exp(-C/T);
    % impurity scattering. 
    imp_scat_inv = D*wwmode.^4;
    % boundary scattering.
    bound_scat_inv = 0; % vmode/E; % (1-p)/(1+p)*(1/dx + 1/dy)*abs(vmode)/2;
    
tot_scat_inv = umk_scat_inv + imp_scat_inv  + bound_scat_inv; % ;%
tot_scat = 1./tot_scat_inv;
% thermal conductivity equation -- Nika 2012 Paper
kall2=[kall' kall' kall' kall' kall' kall'];
cond_mode =(h*wwmode.*vmode).^2.*tot_scat.*(exp(h*wwmode/kB/T)./((exp(h*wwmode/kB/T)-1).^2)).*kall2.^2*dk/(6*pi*kB*T^2); 
cond_mode_sum = sum(abs(cond_mode(1:kindex,:))); 
cond_tot = sum(cond_mode_sum(1:nbranch));
mode_contri = cond_mode_sum./cond_tot;
% figure
plot(mode_contri)
% End of Nika thermal conductivity.
% thermal conductivity equaiton --- Davis 2013 paper. 
% % % for i = 1: kindex;
% % %     kall_mode(i,:) = kall(i)*ones(1,nbranch);
% % % end
% % % Cph = kB*(h*wwmode/kB/T).^2.*(exp(h*wwmode/kB/T)./((exp(h*wwmode/kB/T)-1).^2));
% % % cond_mode_Davis = Cph.*(vmode.^2).*tot_scat.*kall_mode.^2*4*pi/3/((2*pi)^3)*dk;
% % % cond_tot_Davis = sum(sum(abs(cond_mode_Davis(2:kindex,:))));
%*****************End : Callaway Model for Thermal conductivity**********
%#############Plot thermal conductivity based colorful dispersion curve. 
cond_norm = cond_mode/max(max(cond_mode));
color=colormap(jet);
indexq=floor(cond_norm*63)+1;
% figure
for i=1:kindex;
    
    for j=1:nbranch
        plot(kall(i),freq(i,j),'s','MarkerEdgeColor',color(indexq(i,j),:),'MarkerFaceColor',color(indexq(i,j),:),'Markersize',2);

        hold on; 
    end
end
colorbar
caxis([min(min(cond_mode)) max(max(cond_mode))])
set(gca,'Linewidth',3.0)
ylabel('Phonon Frequency (Hz)')
xlabel('kx')
fprintf('Thermal Conductivity k = %f W/mK \n', sum(sum(cond_mode)))
KvsT(n) = sum(sum(cond_mode));
end
figure;
plot(5:5:400, KvsT)

%%%%%%%% Accumulated frequency dependent thermal conductivity%%%%%%%%%%%%%%%
% Calculate bins of frequency. frqbin
cmknorm = cond_norm;
nbin = 200;
dfreq = (max(max(freq) - 0)/nbin);
Normbin = round(freq./dfreq)+1; 
Kbin = zeros(nbin+1,1);
for i = 1: length(cmknorm(:,1))
    for j = 1: length(cmknorm(1,:))
Kbin(Normbin(i,j)) = Kbin(Normbin(i,j)) + cmknorm(i,j);
    end
end


Freqbin = (0:1:nbin)*dfreq;
figure;
plot(Kbin/sum(sum(cmknorm)),Freqbin);
set(gca,'Linewidth',3.0)
xlabel('Normalized Thermal Conductivity');
ylabel('Frequency')
