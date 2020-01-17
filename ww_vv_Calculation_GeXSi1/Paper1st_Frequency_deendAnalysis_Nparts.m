% # # # DISPERSION CURVE CODE CAN SWITCH BETWEEN PRIMITIVE CELL AND
% SUPERCELL # # # % % %
% clear all;
% ww = load('ww.mat');
% vv = load('vv.mat');
% % %% Read all csv files.
ncx = 2;
Nk = 120;
N = Nk;
Ptotal = 20;

Pendflag = 1:Ptotal; %  
nbranch = ncx^3*8*3;
Psize2 = N^3/Ptotal;

% % % ww = zeros(nbranch, N^3);
% % % vv = zeros(nbranch, N^3);
% % % 
% % %     for Pendflag1 = 1:length(Pendflag);
% % %         kstart = 1+(Pendflag1-1)*Psize2;
% % %         kend = Psize2+(Pendflag1-1)*Psize2;
% % %         fname31 = sprintf('ww_Si%dGe1_Nk_%d_P%d_%d.csv', ncx, N, Pendflag1,Ptotal);
% % %         ww(:,kstart:kend) = load(fname31,'-ascii');
% % %         fname32 = sprintf('vv_Si%dGe1_Nk_%d_P%d_%d.csv', ncx, N, Pendflag1,Ptotal);
% % %         vv(:,kstart:kend) = load(fname32,'-ascii');
% % %     end
    
% % %% END of Read all csv files. 
% % % 
% % % 
% % % clear;
% % % clc;
% % % tic


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
% ncx = 5; %  3; %3; % 3; %  3; % 
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
ncxdot = 1; 
ncydot = ncxdot;
nczdot = ncxdot;

mI = ones(1,na); % massGe/massSi*ones(1,na); 
atype = ones(1,na); 

% % % ***** Si5x5x5 with embedded Ge3x3x3 Super Cell.
for i=1:length(xs);
    xdot = (xs(i)<((ncx+ncxdot)/2*A - 0.01)) && (xs(i)>=((ncx-ncxdot)/2*A-0.01));
    ydot = (ys(i)<((ncy+ncydot)/2*A - 0.01)) && (ys(i)>=((ncy-ncydot)/2*A-0.01));
    zdot = (zs(i)<((ncz+nczdot)/2*A - 0.01)) && (zs(i)>=((ncz-nczdot)/2*A-0.01));
    if(xdot && ydot && zdot)
         mI(i) = massGe/massSi; %  1; % 
         atype(i) = 2;
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

% Nk = 30; % 40; %20; %  100; % 60; % 120; %    80; %20; % 100; % 16; %40; %   [60]; %   
% Nresult = zeros(length(Nk),8);

           ncount = 0; 

% for ink = 1:length(Nk);
    N = Nk; % (ink)

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


    T= 300;   

    % % % FOR BIN Calculation; 
nbins = 200;
binEdges = linspace(5.0e10,1.12e+14,nbins+1);
v2tao_1_bin2 = zeros(length(Nk),nbins);
v2_bin2 = zeros(length(Nk),nbins);
tao_1_bin2 = zeros(length(Nk),nbins);
% % % FOR BIN Calculation - - END; 
    
% % % Ptotal = 60; % 1; % 40; % 
% % % Pendflag = 38; % [34:36]; % [39:41,38]; % [45:47]; % [49:50]; % [53:55]; % [57:60]; % [37:42]; % [43:47]; % [48:50]; % [56:60]; %[51:55]; %  [42:55]; % [24:41]; % 1:Ptotal; % [17:20]; % [21:30]; %40:Ptotal,3:30]; %  
% % % 
% % %     for Pendflag1 = 1:length(Pendflag); 
% % %         Pn = Pendflag(Pendflag1);
% % % ksumx2 = zeros(length(kxx),1);
% % % ksumx_NBS = zeros(length(kxx),1);
% % % ksumx_3Part = zeros(length(kxx),1);
% % % ksumx_4Gillet = zeros(length(kxx),1);
ksumx_4Gilletc = zeros(length(kxx),1);
Taosumx_Uprocess = zeros(length(kxx),1);
Taosumx_inc = zeros(length(kxx),1);
Taosumx_UprocINV = zeros(length(kxx),1);
Taosumx_incINV = zeros(length(kxx),1);
Taosumx_Uprocess_Sphere = zeros(length(kxx),1); 
Taosumx_inc_Sphere = zeros(length(kxx),1); 
Taosumx_UprocINV_Sphere = zeros(length(kxx),1); 
Taosumx_incINV_Sphere = zeros(length(kxx),1); 

ksumx_tao2 = zeros(length(kxx),1);
ksumx_tao3 = zeros(length(kxx),1);
ksumx_tao4 = zeros(length(kxx),1);
tao_1_sum = zeros(length(kxx),1);
tao_INV_1_sum = zeros(length(kxx),1);
tao_2_sum = zeros(length(kxx),1);
tao_edge = zeros(length(kxx),1);
v2tao_1_sum = zeros(length(kxx),1);
v2tao_2_sum = zeros(length(kxx),1);
v2tao_edge_sum = zeros(length(kxx),1);
v2_1_sum = zeros(length(kxx),1);
v2_2_sum = zeros(length(kxx),1);
v2_sphere = zeros(length(kxx),1);
v2_edge = zeros(length(kxx),1);

Cpsumnk = zeros(length(kxx),1);
Cpsumnkc = zeros(length(kxx),1);
Cpsumnk_sphere = zeros(length(kxx),1);
Cpsumnk_edge = zeros(length(kxx),1);
timenk = zeros(length(kxx),1);
Nresult = zeros(1,19);
ksum = 0;
%      mypool = parpool(3);
% [fname1, fname2, kstart, kstop, Pn, kindexendflag] = Towne62paralleloadQDGillet(ncx,Nk,kindex,Pn);

%     fname2 = sprintf('Fig3_QDGillet_Ge%dSi1_nk_%d_Vel_cent_diff_P%d_100_DS_Sphere_wotao.csv', ncx, N, Pn(end));



% ww = zeros(nbranch,bandwidth);
% vv = zeros(nbranch,bandwidth);
% nbranch = size(ww,1);
% Psize = size(ww,2);
disp(nbranch)
ninc = 0; 
Pendflag = 1:Ptotal; %
nk_MAX = length(kxx);
  for Pendflag1 = 3; % 16:Ptotal; % 75;
        v2tao_1_bin = zeros(length(Nk),nbranch);
        v2_bin = zeros(length(Nk),nbranch);
        tao_1_bin = zeros(length(Nk),nbranch);
        binIdx_N = zeros(length(Nk),nbranch);

        Pn = Pendflag(Pendflag1);

    kindexendflag = 0; 
    bandwidth = Psize2; % length(kxx)/Ptotal;
    kstart = (1)+(Pn(end)-1)*(bandwidth); 
    kstop = (bandwidth)+(Pn(end)-1)*bandwidth;
    if(kstart > length(kxx))
       fprintf('End of all 3rdIFC pairs'); 
       kindexendflag = 1; 
       return;
    else if(kstop > nk_MAX);
            kstop = nk_MAX;
        end
    end

    fname1 = sprintf('Ave_v2tao_P%d_T%d.csv', Pn(end),Ptotal);
    fname2 = sprintf('Ave_v2_P%d_T%d.csv', Pn(end),Ptotal);
    fname3 = sprintf('Ave_tao_P%d_T%d.csv', Pn(end),Ptotal);
    disp(Pn)
    ww2 = ww.ww; 
    vv2 = vv.vv; 
parfor nk = kstart:kstop; % 1: length(kxx); % 

   if(mod(nk,2000)==0)
   	disp(nk/length(kxx)) 
   end

    kx=kxx(nk);
    ky=kyy(nk); 
    kz=kzz(nk); 
    
    % % % Direct Summation in Sphere BZ.
    
    ninc = ninc + 1;
    wwmodec = ww2(:,nk);% ww.ww(:,nk);
    vmode = vv2(:,nk); % vv.vv(:,nk);

      wwmode = wwmodec; % wwI*sqrt(1/massSi); 
      [Nbin,edges,binIdx] = histcounts(wwmode, [binEdges(1:end-1) Inf]);
        Cpmode = kB*(h*wwmode/kB/T).^2.*(exp(h*wwmode/kB/T)./((exp(h*wwmode/kB/T)-1).^2));
        Cpsum = Cpsum + sum(sum((Cpmode))); 
        Cpsumnk(nk) = sum(sum((Cpmode)));
% 		   Cpmodec = kB*(h*wwmodec/kB/T).^2.*(exp(h*wwmodec/kB/T)./((exp(h*wwmodec/kB/T)-1).^2));
        Cpmodec = Cpmode; 
		   Cpsumnkc(nk) = sum(sum((Cpmodec)));

        
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
    
        k_reduced = vmode.*vmode.*Cpmodec; %


       
   %  Uprocess_scat_inv = bu*(wwmode.*wwmode)./(vmode.*vmode);
    gamma = 1.5; theta = 645; % Debye temperature, unit: K
    Uprocess_scat_inv = (h*gamma^2*T/massSi/theta*exp(-theta/3/T))*(wwmode.*wwmode)./(vmode.*vmode);
    dA = 1.5849; dK = 0.24;  % % % data from Gillet-JHT-2009. 
    [inc_scat_inv] = incoherent_scat(ncx,ncxdot,dA,dK,kx,ky,kz,a0,vmode,sigmaSi);
    
% % %     tot_scat_inv_Gillet = Uprocess_scat_inv + inc_scat_inv;
% % %     tot_scat_Gillet = 1./tot_scat_inv_Gillet;
% % %     kmode_holland4_Gillet = vmode.*vmode.*tot_scat_Gillet.*Cpmode; %
% % %     ksumx_4Gillet(nk) = sum(sum(kmode_holland4_Gillet))*8/V;
    
		   Uprocess_scat_invc = (h*gamma^2*T/massSi/theta*exp(-theta/3/T))*(wwmodec.*wwmodec)./(vmode.*vmode);
		   tot_scat_inv_Gilletc = Uprocess_scat_invc + inc_scat_inv;
		   tot_scat_Gilletc = 1./tot_scat_inv_Gilletc;
		   kmode_holland4_Gilletc = k_reduced.*tot_scat_Gilletc; % vmode.*vmode.*tot_scat_Gilletc.*Cpmodec; %
		   ksumx_4Gilletc(nk) = sum(sum(kmode_holland4_Gilletc))*8/V; % % tao_1; 
           
% % %            Taosumx_Uprocess(nk) = sum(sum((1./Uprocess_scat_invc)))*8/V;
% % %            Taosumx_inc(nk) = sum(sum((1./inc_scat_inv)))*8/V;
% % %            Taosumx_UprocINV(nk) = sum(sum((Uprocess_scat_invc)))*8/V;
% % %            Taosumx_incINV(nk) = sum(sum((inc_scat_inv)))*8/V;
           
           
           tao_1 = tot_scat_Gilletc; 
% % %            tao_INV_1 = tot_scat_inv_Gilletc; 
           tao_3 = 1; 
    if(sqrt(kx^2+ky^2+kz^2)>XX(end))
%         continue;
        tao_2 = 0; 
        tao_2_edge = tot_scat_Gilletc;
        tao_4 = 0; 
% % %         Cpsumnk_sphere(nk) = 0;
% % %         Cpsumnk_edge(nk) = Cpsumnk(nk);
% % %         v2_sphere(nk) = 0; 
% % %         v2_edge(nk) = sum(vmode.*vmode);
% % %         Taosumx_Uprocess_Sphere(nk) = 0; % sum(sum(k_reduced.*(1./Uprocess_scat_invc)))*8/V;
% % %         Taosumx_inc_Sphere(nk) = 0; % sum(sum(k_reduced.*(1./inc_scat_inv)))*8/V;
% % %         Taosumx_UprocINV_Sphere(nk) = 0; % sum(sum(k_reduced.*(1./Uprocess_scat_invc)))*8/V;
% % %         Taosumx_incINV_Sphere(nk) = 0; % sum(sum(k_reduced.*(1./inc_scat_inv)))*8/V;
    else
        ncount = ncount+1;
        tao_2 = tot_scat_Gilletc;
        tao_2_edge = 0; 
        tao_4 = 1; 
% % %         Cpsumnk_sphere(nk) = Cpsumnk(nk);
% % %         Cpsumnk_edge(nk) = 0; 
% % %         v2_sphere(nk) = sum(vmode.*vmode);
% % %         v2_edge(nk) = 0;
% % %         Taosumx_Uprocess_Sphere(nk) = sum(sum((1./Uprocess_scat_invc)))*8/V;
% % %         Taosumx_inc_Sphere(nk) = sum(sum((1./inc_scat_inv)))*8/V;
% % %         Taosumx_UprocINV_Sphere(nk) = sum(sum((1./Uprocess_scat_invc)))*8/V;
% % %         Taosumx_incINV_Sphere(nk) = sum(sum((1./inc_scat_inv)))*8/V;
    end
    tao_1_sum(nk) = sum(tao_1); % sqrt(tao_1'*tao_1);
% % %     tao_INV_1_sum(nk) = sum(tao_INV_1); % sqrt(tao_1'*tao_1); 
% % %     tao_2_sum(nk) = sum(tao_2); % sqrt(tao_2'*tao_2);
% % %     tao_edge(nk) =sum(tao_2_edge);
    v2tao_1_temp = vmode.*vmode.*tao_1;
    v2tao_1_sum(nk) = sum(v2tao_1_temp); % sum(vmode.*vmode.*tao_1);
% % %     v2tao_2_sum(nk) = sum(vmode.*vmode.*tao_2);
% % %     v2tao_edge_sum(nk) = sum(vmode.*vmode.*tao_2_edge);
    v2_temp = vmode.*vmode; 
% %     v2_1_sum(nk) = sum(v2_temp); % sum(vmode.*vmode);
% %     v2_2_sum(nk) = sum(v2_temp); % sum(vmode.*vmode);
    
% % %     kmode_tao2 = k_reduced.*tao_2; %  vmode.*vmode.*tao_2.*Cpmodec; %
% % %     ksumx_tao2(nk) = sum(sum(kmode_tao2))*8/V; % 
% % %     kmode_tao3 = k_reduced.*tao_3; %  vmode.*vmode.*tao_3.*Cpmodec; %
% % %     ksumx_tao3(nk) = sum(sum(kmode_tao3))*8/V; % 
% % %     kmode_tao4 = k_reduced.*tao_4; %  vmode.*vmode.*tao_4.*Cpmodec; %
% % %     ksumx_tao4(nk) = sum(sum(kmode_tao4))*8/V; % 

% %     ww(:,ninc) = wwmode;
% %     vv(:,ninc) = vmode;

%       fprintf(fileID,'%6d %12.8f %12.8f %12.8f %12.8f %12.8f\n', nk, Cpsumnk(nk),ksumx2(nk), ksumx_NBS(nk), ksumx_3Part(nk),ksumx_4Gillet(nk));
    
    tao_1_temp = tao_1; 
   %     for i = 1:nbranch; 
%                      test =    v2tao_1_temp(binIdx(1:nbranch),1);
% %         v2tao_1_bin(nk,binIdx(1:nbranch)) = v2tao_1_temp(1:nbranch),1)';
% %         v2_bin(nk,:) = v2_temp(1:nbranch,1)';
% %         tao_1_bin(nk,:) = tao_1_temp(1:nbranch,1)';
% %         binIdx_N(nk,:) = binIdx(1:nbranch);
        % % -----------------
        v2tao_1_bin(nk,:) = v2tao_1_temp(1:nbranch,1)';
        v2_bin(nk,:) = v2_temp(1:nbranch,1)';
        tao_1_bin(nk,:) = tao_1_temp(1:nbranch,1)';
        binIdx_N(nk,:) = binIdx(1:nbranch);
%     end
%     fname1 = sprintf
end
    for nk = kstart:kstop
        v2tao_1_bin2(nk,binIdx_N(nk,1:nbranch)) = v2tao_1_bin(nk,1:nbranch); % v2tao_1_temp(1:nbranch,1);
        v2_bin2(nk,binIdx_N(nk,1:nbranch)) = v2_bin(nk,1:nbranch); % v2_temp(1:nbranch,1)';
        tao_1_bin2(nk,binIdx_N(nk,1:nbranch)) = tao_1_bin(nk,1:nbranch);
    end
csvwrite(fname1,v2tao_1_bin2);
csvwrite(fname2,v2_bin2);
csvwrite(fname3,tao_1_bin2);
  end
disp('End of these current parts')
% pause; 

% % % % % 
% % % % % Ave_v2_bin = sum(v2_bin,1)*8/V; 
% % % % % Ave_v2tao_1_bin = sum(v2tao_1_bin,1)*8/V;
% % % % % Ave_tao_1_bin = sum(tao_1_bin,1)*8/V;
% % % % % save('Ave_v2_bin.mat','Ave_v2_bin');
% % % % % save('Ave_v2tao_1_bin.mat','Ave_v2tao_1_bin');
% % % % % save('Ave_tao_1_bin.mat','Ave_tao_1_bin');
% % % % % 
% % % % % Ave_Tao_FullBZ = sum(tao_1_sum)*8/V % sum(tao_1_sum)/nbranch/length(kxx) % ./nbranch)/size(tao_1_sum,1)
% % % % % Ave_Tao_INV_FullBZ = sum(tao_INV_1_sum)*8/V % sum(tao_1_sum)/nbranch/length(kxx) % ./nbranch)/size(tao_1_sum,1)
% % % % % 
% % % % % % % Ave_Tao_SphereBZ = sum(tao_2_sum)*8/V % sum(tao_2_sum)/nbranch/ncount
% % % % % % % Ave_Tao_edge = sum(tao_edge)*8/V
% % % % % 
% % % % % Ave_v2Tao_FullBZ = sum(v2tao_1_sum)*8/V % sum(v2tao_1_sum)/nbranch/length(kxx);
% % % % % % % Ave_v2Tao_SphereBZ = sum(v2tao_2_sum)*8/V % sum(v2tao_2_sum)/nbranch/ncount;
% % % % % % % Ave_v2Tao_edge = sum(v2tao_edge_sum)*8/V
% % % % % 
% % % % % Ave_v2_FullBZ = sum(v2_1_sum)*8/V % sum(v2tao_1_sum)/nbranch/length(kxx);
% % % % % Ave_v2_FullBZ2 = sum(v2_2_sum)*8/V % sum(v2tao_2_sum)/nbranch/ncount;
% % % % % % % Ave_v2_SphereBZ = sum(v2_sphere)*8/V
% % % % % % % Ave_v2_edge = sum(v2_edge)*8/V
% % % % % 
% % % % % 
% % % % % 
% % % % % % % % Cpnk = sum(Cpsumnk)*8/V/2329;
% % % % % % % % Cpnkc = sum(Cpsumnkc)*8/V/2329;
% % % % % % % % Cpk_sphereBZ = sum(Cpsumnk_sphere)*8/V/2329;
% % % % % % % % Cpk_edge = sum(Cpsumnk_edge)*8/V/2329;
% % % % % % % % 
% % % % % % % % k_Gilletc = sum(sum(ksumx_4Gilletc))
% % % % % % % % k_tao2 = sum(sum(ksumx_tao2))
% % % % % % % % k_tao3 = sum(sum(ksumx_tao3))
% % % % % % % % k_tao4 = sum(sum(ksumx_tao4))
% % % % % % % % Tao_Uprocess = sum(sum(Taosumx_Uprocess));
% % % % % % % % Tao_inc = sum(sum(Taosumx_inc));
% % % % % % % % Tao_UprocINV = sum(sum(Taosumx_UprocINV));
% % % % % % % % Tao_incINV = sum(sum(Taosumx_incINV));
% % % % % % % % Tao_Uprocess_Sphere = sum(sum(Taosumx_Uprocess_Sphere));
% % % % % % % % Tao_inc_Sphere = sum(sum(Taosumx_inc_Sphere));
% % % % % % % % Tao_UprocINV_Sphere = sum(sum(Taosumx_UprocINV_Sphere));
% % % % % % % % Tao_incINV_Sphere = sum(sum(Taosumx_incINV_Sphere));
% % % % % % % %  k_Gilletc2 = sum(ksumx_4Gilletc/8*V)*1/(N*N*N*(A1*sigmaSi)^3)
% % % % % % % %  k_tao2_2 = sum(ksumx_tao2/8*V)*1/(ncount*(A1*sigmaSi)^3)
% % % % % % % % Nresult = [N,0,0,0,0,0,k_Gilletc,Ave_Tao_FullBZ,Ave_Tao_SphereBZ,Ave_v2Tao_FullBZ, Ave_v2Tao_SphereBZ,Ave_v2_FullBZ,Ave_v2_FullBZ2, Ave_v2_SphereBZ,Cpnk,Cpnkc,Cpk_sphereBZ, Cpk_edge, Ave_v2_edge,Ave_Tao_edge, Ave_v2Tao_edge];
% % % % % % % % [Tao_Uprocess,Tao_inc,Tao_Uprocess_Sphere,Tao_inc_Sphere, Tao_UprocINV,Tao_incINV,Tao_UprocINV_Sphere,Tao_incINV_Sphere]
% % % % % % % % Ave_Tao_INV_FullBZ
% % % % % 
% % % % % figure
% % % % % aj = binEdges(1:end-1);     %# bins lower edge
% % % % % bj = binEdges(2:end);       %# bins upper edge
% % % % % cj = ( aj + bj ) ./ 2;      %# bins center
% % % % % 
% % % % % %# plot histogram
% % % % % bar(cj,Ave_v2tao_1_bin,'hist')
% % % % % % set(gca, 'XTick',binEdges, 'XLim',[binEdges(1):binEdges(end)])
% % % % % xlabel('Bins'), ylabel('Counts'), title('histogram of v2tao')
% % % % % figure(3);
% % % % % hold on;
% % % % % plot(Ave_v2_bin);
% % % % % figure(5)
% % % % % hold on;
% % % % % plot(Ave_tao_1_bin,'k*--');
% % % % % hold on;
