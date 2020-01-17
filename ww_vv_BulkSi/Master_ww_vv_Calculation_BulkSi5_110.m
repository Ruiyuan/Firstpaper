% # # # DISPERSION CURVE CODE CAN SWITCH BETWEEN PRIMITIVE CELL AND
% SUPERCELL # # # % % %

clear;
clc;
tic
addpath('S:\S-current-debug-folder\C-H_mode_08192015\CHModel_Primitive');
addpath('S:\S-current-debug-folder\C-H_mode_08192015\CHModel_Primitive\fmatSuperCell');
%----------------------------------------------
    lambdaSi = 21.0; lambdaGe = 31.0; 
    epsSi = 3.473928e-19; epsGe = 3.085e-19; 
    sigmaSi = 2.0951e-10; sigmaGe = 2.0951e-10; 
    massSi = 4.6637e-26; massGe = 1.2057e-25; 
    [eps,isigma,leps,epssig, imass] = SiGein(lambdaSi, lambdaGe, epsSi, epsGe, sigmaSi, sigmaGe, massSi, massGe);
    gamma=1.20d0;
    a0=5.43095e-10;

A = a0/sigmaSi; 
ncell=2;
nbasis=4;
naxis=3;
forceflag=4;  
rc=1.8;

mI=1;
mTmatrix=[massGe/massSi]; 
%---------start--contstruct "force matrix" of pure Si or heavy Si-------
ncx = 2; % 3; %    3; % 5; %  
ncy = ncx; 
ncz = ncx; 

ncxlarge = ncx;  
ncxsuper = ncx;
ncxvoid = 0; 

% % % SWITCH BETWEEN PRIMITIVE CELL AND SUPERCELL!!!
        flagsupercell = 1;

if(ncxvoid>1)
    [xs,ys,zs,nsupercell]=creatSuperCellwMultipleVoid(A,ncx,ncy,ncz,ncell,nbasis,ncxsuper,ncxvoid,flagsupercell);   
else
    [xs,ys,zs,nsupercell]=creatSuperCellwVoid(A,ncx,ncy,ncz,ncell,nbasis,ncxsuper,ncxvoid,flagsupercell);   
end

na=length(xs); 
% % % **** Add Quantum Dots here ****% 
ncxdot = 0; 
ncydot = ncxdot;
nczdot = ncxdot;

mI = ones(1,na); % (massGe/massSi).*ones(1,na); 
atype =ones(1,na); %  2.*ones(1,na); 

% % % ***** Si5x5x5 with embedded Ge3x3x3 Super Cell.
for i=1:length(xs);
    xdot = (xs(i)<((ncx+ncxdot)/2*A - 0.01)) && (xs(i)>=((ncx-ncxdot)/2*A-0.01));
    ydot = (ys(i)<((ncy+ncydot)/2*A - 0.01)) && (ys(i)>=((ncy-ncydot)/2*A-0.01));
    zdot = (zs(i)<((ncz+nczdot)/2*A - 0.01)) && (zs(i)>=((ncz-nczdot)/2*A-0.01));
    if(xdot && ydot && zdot)
         mI(i) = massGe/massSi; % % 1; % 
         atype(i) = 2; % 1; % 
    end
end
   
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
        ncxlarge; 
        ncy = ncxlarge;
        ncz = ncxlarge;
        % % % SWITCH BETWEEN PRIMITIVE CELL AND SUPERCELL!!!
        if(flagsupercell == 1)
            A1 = A*ncxsuper; % *ncx; % ; %
        else 
            A1 = A*1;
        end
        
load('fmatSi.mat');
fmatSi = -fmatSi;

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

h = 1.054571726e-34; % 
kB = 1.3806488e-23;
% % % *** Davis 2011 Data
B = 2.1e-19; 
C = 180;
D = 1.32e-45;
n=0;

Nk = 120; % 30; % 120; %80; % 60; %     80; %20; % 100; % 16; %40; %   [60]; %   
% Nresult = zeros(length(Nk),8);


for ink = 1:length(Nk);
    N = Nk(ink)

if(flagsupercell == 1)
    XX = ((0:1/N:1-1/N)+1/N/2)*1*pi/A1; 
    YY = ((0:1/N:1-1/N)+1/N/2)*1*pi/A1; 
    ZZ = ((0:1/N:1-1/N)+1/N/2)*1*pi/A1; 
nbranch = nsupercell*3; 

else 
    XX = ((0:1/N:1-1/N)+1/N/2)*2*pi/A1; 
    YY = ((0:1/N:1-1/N)+1/N/2)*2*pi/A1; 
    ZZ = ((0:1/N:1-1/N)+1/N/2)*2*pi/A1; 
nbranch = ncell*3; 

end
end

dkx = (XX(2) - XX(1))/sigmaSi;
dky = (YY(2) - YY(1))/sigmaSi;  
dkz = (ZZ(2) - ZZ(1))/sigmaSi;  

Lx = 2*pi/(dkx); Ly = 2*pi/(dky); Lz = 2*pi/(dkz);
V = Lx*Ly*Lz;

 wwindex=0;
nxmax=length(XX);
nymax=length(YY);
nzmax=length(ZZ);

nk = 1;
for nx= 1:nxmax; 
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
nxmax=length(kxx);
nymax=length(kyy);
nzmax=length(kzz);

kindex = nxmax;


    T= 300  
     mypool = parpool(4);

Ptotal = 120; % 450; % 80; % % 
Pendflag = [86:Ptotal]; % TMC62-05: 81:95; % 
%VPC-: 66:80; % TM62-06: 51:65; % VPC-20: 41:50; % VPC-09: 31:40; % VPC-30: 21:30; %  VPC-24: 11:20; % VPC-16: 1:10; %  
for Pendflag1 = 1:length(Pendflag); 
        Pn = Pendflag(Pendflag1);
% % % ksumx2 = zeros(length(kxx),1);
% % % ksumx_NBS = zeros(length(kxx),1);
% % % ksumx_3Part = zeros(length(kxx),1);
% % % ksumx_4Gillet = zeros(length(kxx),1);
ksumx_4Gilletc = zeros(length(kxx),1);
ksumx_tao2 = zeros(length(kxx),1);
ksumx_tao3 = zeros(length(kxx),1);
ksumx_tao4 = zeros(length(kxx),1);

Cpsumnk = zeros(nbranch,1);
Cpsumnkc = zeros(nbranch,1);
timenk = zeros(nbranch,1);
Nresult = zeros(1,11);
ksum = 0;
% [fname1, fname2, kstart, kstop, Pn, kindexendflag] = Towne62paralleloadQDGillet(ncx,Nk,kindex,Pn);
    kindexendflag = 0; 
    bandwidth = kindex/Ptotal;
    kstart = (1)+(Pn(end)-1)*(bandwidth); 
    kstop = (bandwidth)+(Pn(end)-1)*bandwidth;
    if(kstart > kindex)
       fprintf('End of the ncx%d_Nk%d',ncx, kindex); 
       kindexendflag = 1; 
       return;
    else if(kstop > kindex);
            kstop = kindex;
        end
    end
%     fname1 = sprintf('Fig3_QDGillet_N%d_nk_%d_ksum_Vel_cent_diff_P%d.csv', ncx, N, Pn(end));

    fname2 = sprintf('Fig3_QDGillet_N%d_nk_%d_Vel_cent_diff_P%d_dot1_100_DS_Sphere_wotao.csv', ncx, N, Pn(end));


if(kindexendflag == 1);
    continue;
end
disp(Pn)



% ninc = 0; 
wwtotal = zeros(nbranch, length(kxx));
vvtotal = zeros(nbranch, length(kxx));
parfor nk = kstart:kstop; 
% for nk = kstart:kstop; 
%     ww = zeros(nbranch,bandwidth);
%     vv = zeros(nbranch,bandwidth);
    
   if(mod(nk,100)==0)
   	disp(nk/length(kxx)) 
   end

    kx=kxx(nk);
    ky=kyy(nk); 
    kz=kzz(nk); 
    
    % % % Direct Summation in Sphere BZ.
    
%     ninc = ninc + 1;
%    Dsuper = zeros(naxis*nsupercell, naxis*nsupercell);

% Dall = Phiall.*exp(1i*(kx*rxall+ky*ryall+kz*rzall));

% for i = 1:na/nsupercell;
% Dsuper = Dsuper + Dall(:,(i-1)*naxis*nsupercell+(1:naxis*nsupercell));
% end

%   [eV,SS]=svd(Dsuper);
%  wwI = sqrt(diag(SS));
%  ww(nk,:) = wwI;

   dkxf = 0.001*1*pi/A1; % 

    Dsuper = zeros(naxis*nsupercell, naxis*nsupercell);
    kxf = (kx+dkxf); % (kx + 0.025*1*pi/(A1)); %(kx+0.01); %  
    kyf = ky +dkxf; 
    Dall = Phiall.*exp(1i*(kxf*rxall+kyf*ryall+kz*rzall)); 
  
    for i = 1:na/nsupercell;
        Dsuper = Dsuper + Dall(:,(i-1)*naxis*nsupercell+(1:naxis*nsupercell));
    end
   [eVf,SS]=svd(Dsuper);
    wwf=(sqrt(diag(SS))); 

% % % % % ***Backward Step. **********************    
     Dsuper = zeros(naxis*nsupercell, naxis*nsupercell);    
     kxb = (kx - dkxf);
     kyb = ky - dkxf; 
     
     Dall = Phiall.*exp(1i*(kxb*rxall+kyb*ryall+kz*rzall));
     for i = 1:na/nsupercell;
         Dsuper = Dsuper + Dall(:,(i-1)*naxis*nsupercell+(1:naxis*nsupercell));
     end 
     [eVb,SS]=svd(Dsuper);
%     eVb=UU ; % 
     wwb=(sqrt(diag(SS))); 
% % % % % ***Backward Step. END**********************    

    
    Dall = []; Dsuper = []; UU = []; SS = []; VV = [];

% % % % % % % % ### GROUP VELOCITY ***CENTRAL DIFFERENTIAL***######## 
     checkeV = abs((eVb'*eVf));
     [M, ip] = max(checkeV');
		   dkxvg = sqrt((kxf - kxb)^2 + (kyf-kyb)^2);
     vmode = (wwb - wwf(ip))*sqrt(1/massSi)/(1*dkxvg/sigmaSi);
		   wwc = (wwb+wwf(ip))/2;
		   wwmodec = wwc*sqrt(1/massSi);
% % % 		   wwI = wwc; 
     eVf = []; eVb = [];  
% % % % % % % % ### GROUP VELOCITY ***BACKWARD***######## 


% % % % % ### GROUP VELOCITY ***FORWARD***######## 
%    checkeV = abs((eV'*eVf));
%    [M, ipf] = max(checkeV');
%    dkxvg = (kxf - kx);
%    vmodef = (wwI - wwf(ipf))*sqrt(1/massSi)/(1*dkxvg/sigmaSi); % (1*0.001*kx/sigmaSi); % (1*0.025*1*pi/A1/sigmaSi); % (1*0.01/sigmaSi); % 
% % % % % ### GROUP VELOCITY ***RORWARD***######## 
%    vmode = vmodef;
           
% % wwxx(nx,:) = wwI*sqrt(1/massSi);  % wwxx(nz,nx,:) = wwI; 

      wwmode = wwmodec; % wwI*sqrt(1/massSi); 
% % % % % %         Cpmode = kB*(h*wwmode/kB/T).^2.*(exp(h*wwmode/kB/T)./((exp(h*wwmode/kB/T)-1).^2));
% % % % % %         Cpsum = Cpsum + sum(sum((Cpmode))); 
% % % % % %         Cpsumnk(nk) = sum(sum((Cpmode)));
% % % % % % % 		   Cpmodec = kB*(h*wwmodec/kB/T).^2.*(exp(h*wwmodec/kB/T)./((exp(h*wwmodec/kB/T)-1).^2));
% % % % % %         Cpmodec = Cpmode; 
% % % % % % 		   Cpsumnkc(nk) = sum(sum((Cpmodec)));
% % % % % % 
% % % % % %         
% % % % % %             umk_scat_inv = B*(wwmode.*wwmode)*T*exp(-C/T);
% % % % % %     % impurity scattering. 
% % % % % %     imp_scat_inv = D*wwmode.^4;
% % % % % %     % boundary scattering.
% % % % % %     bound_scat_inv = abs(vmode)/L; % E; %  (vmode)/E; %0; %  (1-p)/(1+p)*(1/dx + 1/dy)*abs(vmode)/2;
% % % % % %     % partical scattering. from Majumdar-1993
% % % % % %     rvoid = ncxvoid*a0/2;
% % % % % %     chi = rvoid*sqrt(kx^2+ky^2+kz^2)/sigmaSi;
% % % % % %     sig_void = pi*rvoid^2*(chi^4/(chi^4+1));
% % % % % %     rho_void = (ncxvoid/ncx)^3;
% % % % % %     part_scat_inv = sig_void*rho_void*abs(vmode);
% % % % % %     
% % % % % %         k_reduced = vmode.*vmode.*Cpmodec; %
% % % % % % 
% % % % % % %     tot_scat_inv = umk_scat_inv + imp_scat_inv  + bound_scat_inv; % ;%
% % % % % % %     tot_scat = 1./tot_scat_inv;
% % % % % % %     kmode_holland = k_reduced.*tot_scat; % vmode.*vmode.*tot_scat.*Cpmode; %  + wy(nk)*vmodey.*vmodey.*tot_scat.*Cpmode + wz(nk)*vmodez.*vmodez.*tot_scat.*Cpmode; %
% % % % % % %     ksumx2(nk) = sum(sum(kmode_holland))*8/V;
% % % % % % % % %     % % % *** Don't consider Boundary Scattering. ***
% % % % % % % % %     tot_scat_inv2 = umk_scat_inv + imp_scat_inv;
% % % % % % % % %     tot_scat2 = 1./tot_scat_inv2;
% % % % % % % % %      kmode_holland2_NBS = vmode.*vmode.*tot_scat2.*Cpmode; %
% % % % % % % % %      ksumx_NBS(nk) = sum(sum(kmode_holland2_NBS))*8/V;
% % % % % % % % %     % % % *** END don't consider Boundary scattering. 
% % % % % %     
% % % % % % % % %     % % % *** Consider particle scattering, not Boundary Scattering. 
% % % % % % % % %     tot_scat_inv3 = umk_scat_inv + imp_scat_inv + part_scat_inv;
% % % % % % % % %     tot_scat3 = 1./tot_scat_inv3;
% % % % % % % % %      kmode_holland3_Part = vmode.*vmode.*tot_scat3.*Cpmode; %
% % % % % % % % %      ksumx_3Part(nk) = sum(sum(kmode_holland3_Part))*8/V;
% % % % % % % % %     % % % *** END Particle Scattering, not Boundary scattering. 
% % % % % %     
% % % % % %     % % % *** Gillet: incoherent scattering and Umklapp scattering. 
% % % % % %    
% % % % % %    %  Uprocess_scat_inv = bu*(wwmode.*wwmode)./(vmode.*vmode);
% % % % % %     gamma = 1.5; theta = 645; % Debye temperature, unit: K
% % % % % %     Uprocess_scat_inv = (h*gamma^2*T/massSi/theta*exp(-theta/3/T))*(wwmode.*wwmode)./(vmode.*vmode);
% % % % % %     dA = 1.5849; dK = 0.24;  % % % data from Gillet-JHT-2009. 
% % % % % %     [inc_scat_inv] = incoherent_scat(ncx,ncxdot,dA,dK,kx,ky,kz,a0,vmode,sigmaSi);
% % % % % %     
% % % % % % % % %     tot_scat_inv_Gillet = Uprocess_scat_inv + inc_scat_inv;
% % % % % % % % %     tot_scat_Gillet = 1./tot_scat_inv_Gillet;
% % % % % % % % %     kmode_holland4_Gillet = vmode.*vmode.*tot_scat_Gillet.*Cpmode; %
% % % % % % % % %     ksumx_4Gillet(nk) = sum(sum(kmode_holland4_Gillet))*8/V;
% % % % % %     
% % % % % % 		   Uprocess_scat_invc = (h*gamma^2*T/massSi/theta*exp(-theta/3/T))*(wwmodec.*wwmodec)./(vmode.*vmode);
% % % % % % 		   tot_scat_inv_Gilletc = Uprocess_scat_invc + inc_scat_inv;
% % % % % % 		   tot_scat_Gilletc = 1./tot_scat_inv_Gilletc;
% % % % % % 		   kmode_holland4_Gilletc = k_reduced.*tot_scat_Gilletc; % vmode.*vmode.*tot_scat_Gilletc.*Cpmodec; %
% % % % % % 		   ksumx_4Gilletc(nk) = sum(sum(kmode_holland4_Gilletc))*8/V; % % tao_1; 
% % % % % %            
% % % % % %            tao_3 = 1; 
% % % % % %     if(sqrt(kx^2+ky^2+kz^2)>XX(end))
% % % % % % %         continue;
% % % % % %         tao_2 = 0; 
% % % % % %         tao_4 = 0; 
% % % % % %     else
% % % % % %         tao_2 = tot_scat_Gilletc;
% % % % % %         tao_4 = 1;       
% % % % % %     end
% % % % % %     kmode_tao2 = k_reduced.*tao_2; %  vmode.*vmode.*tao_2.*Cpmodec; %
% % % % % %     ksumx_tao2(nk) = sum(sum(kmode_tao2))*8/V; % 
% % % % % %     kmode_tao3 = k_reduced.*tao_3; %  vmode.*vmode.*tao_3.*Cpmodec; %
% % % % % %     ksumx_tao3(nk) = sum(sum(kmode_tao3))*8/V; % 
% % % % % %     kmode_tao4 = k_reduced.*tao_4; %  vmode.*vmode.*tao_4.*Cpmodec; %
% % % % % %     ksumx_tao4(nk) = sum(sum(kmode_tao4))*8/V; % 

%     ww(:,nk-kstart+1) = wwmode;
%     vv(:,nk-kstart+1) = vmode;
    wwtotal(:,nk) = wwmode;
    vvtotal(:,nk) = vmode;
%       fprintf(fileID,'%6d %12.8f %12.8f %12.8f %12.8f %12.8f\n', nk, Cpsumnk(nk),ksumx2(nk), ksumx_NBS(nk), ksumx_3Part(nk),ksumx_4Gillet(nk));

end
    ww(:,1:bandwidth) = wwtotal(:,kstart:kstop);
    vv(:,1:bandwidth) = vvtotal(:,kstart:kstop);
%%    cd('C:\ww_vv_Calculation_BulkSi');
    fname31 = sprintf('ww_BulkSi%d_Nk_%d_P%d_%d_110.csv', ncx, N, Pn(end),Ptotal);
    save(fname31,'ww','-ascii');
    fname32 = sprintf('vv_BulkSi%d_Nk_%d_P%d_%d_110.csv', ncx, N, Pn(end),Ptotal);
    save(fname32,'vv','-ascii');
%     cd('S:\S-current-debug-folder\C-H_mode_08192015\CHModel_Primitive\test_parallel\ww_vv_Calculation_GeXSi1');
%     fclose(fileID);

end
    delete(gcp('nocreate'))

% ksum = sum(sum(ksumx2));

% k = ksum*1/(N*N*N*(A1*sigmaSi)^3) 
% k2 = 8*ksum*dkx*dky*dkz/(2*pi)^3 
% Cp = Cpsum/(N*N*N*(A1*sigmaSi)^3)/2329 
% 
% Lx = 2*pi/(dkx); Ly = 2*pi/(dky); Lz = 2*pi/(dkz);
% V = Lx*Ly*Lz;
% k3 = ksum; 
% Cp3 = Cpsum*8/V/2329; % 

% %% 
% k_NBS = sum(sum(ksumx_NBS));
% k_Particle = sum(sum(ksumx_3Part)); 

% % % % k_Gillet = sum(sum(ksumx_4Gillet))
% % % k_Gilletc = sum(sum(ksumx_4Gilletc))
% % % k_tao2 = sum(sum(ksumx_tao2))
% % % k_tao3 = sum(sum(ksumx_tao3))
% % % k_tao4 = sum(sum(ksumx_tao4))
% % % 
% % % 
% % % timevalue = toc
% % % 		   Nresult = [Nk(ink),0,0,0,0,0,k_Gilletc,k_tao2,k_tao3,k_tao4, timevalue];
% % % % 		   Nresult = [Nk(ink),Cp3,k_NBS,k3,k_Particle,k_Gillet,k_Gilletc,k_tao2,k_tao3,k_tao4, timevalue];
% % % Cp3 = 0; Cpsum = 0; tic
% % % 
% % % csvwrite(fname2,Nresult);
% % %     end
% % %     
% % %     % % % ****** SUM the Results *****
% % % %  for ink = 1:length(Nk); 
% % % %    N = Nk(ink);
% % % 
% % %     % % % ****** SUM the Results *****
% % %     clear kNP
% % %     % Ptotal = N;
% % %     kNP = zeros(length(Pendflag),11);
% % %     for Pendflag1 = 1:length(Pendflag); % Ptotal; 
% % % 		   Pn = Pendflag(Pendflag1);
% % % 		   fname4 = sprintf('Fig3_QDGillet_N%d_nk_%d_Vel_cent_diff_P%d_dot1_100_DS_Sphere_wotao.csv', ncx, N, Pn(end));
% % %         % Fig3_QDGillet_N%d_nk_%d_Vel_cent_diff_P%d.csv
% % %         f = load(fname4);
% % %         kNP(Pendflag1,:) = f; % zeros(size(f,1),size(f,2));
% % %         Nresult2 = sum(kNP);
% % %     end
% % %     fname4 = sprintf('Fig3_QDGillet_N%d_nk_%d_k_total_cent_dot1_100_DS_Sphere_wotao.csv', ncx, N);
% % %     csvwrite(fname4,Nresult2);
% % % %  end
% % % ****** SUM the Results- The END. 
    

