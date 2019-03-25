% # # # DISPERSION CURVE CODE CAN SWITCH BETWEEN PRIMITIVE CELL AND
% SUPERCELL # # # % % %
clear all;
% % %% Read all csv files.
ncx = 4;
Nk = 50;
N = Nk;
Ptotal = 1;

Pendflag = 1:Ptotal; %  
nbranch = ncx^3*8*3;
Psize2 = N^3/Ptotal;

% ww = zeros(nbranch, N^3);
% vv = zeros(nbranch, N^3);
% 
%     for Pendflag1 = 1:length(Pendflag);
%         disp(Pendflag1)
%         kstart = 1+(Pendflag1-1)*Psize2;
%         kend = Psize2+(Pendflag1-1)*Psize2;
%         fname31 = sprintf('ww_Ge%dSi1_Nk_%d_P%d_%d.csv', ncx, N, Pendflag1,Ptotal);
%         ww(:,kstart:kend) = load(fname31,'-ascii');
%         fname32 = sprintf('vv_Ge%dSi1_Nk_%d_P%d_%d.csv', ncx, N, Pendflag1,Ptotal);
%         vv(:,kstart:kend) = load(fname32,'-ascii');
%     end

%     save(ww.mat,'ww','-ascii')
%     save(vv.mat,'vv','-ascii')
% % %% END of Read all csv files. 
% % % 
% % % 
% % % clear;
% % % clc;
% % % tic
% % ww = load('ww.mat','-ASCII');
% % vv = load('vv.mat','-ASCII');

Ptotal = 1;
Pendflag = 1:Ptotal; %  
nbranch = ncx^3*8*3;
Psize2 = N^3/Ptotal;
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

mI = (massGe/massSi).*ones(1,na); 
atype = 2.*ones(1,na); 

% % % ***** Si5x5x5 with embedded Ge3x3x3 Super Cell.
for i=1:length(xs);
    xdot = (xs(i)<((ncx+ncxdot)/2*A - 0.01)) && (xs(i)>=((ncx-ncxdot)/2*A-0.01));
    ydot = (ys(i)<((ncy+ncydot)/2*A - 0.01)) && (ys(i)>=((ncy-ncydot)/2*A-0.01));
    zdot = (zs(i)<((ncz+nczdot)/2*A - 0.01)) && (zs(i)>=((ncz-nczdot)/2*A-0.01));
    if(xdot && ydot && zdot)
         mI(i) = 1; % massGe/massSi; % % 
         atype(i) = 1; % 2;
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

% % % FOR BIN Calculation - - END; 

ksumx_4Gilletc = zeros(length(kxx),1);
ksumx_SBZ = zeros(length(kxx),1);
ksumx_110SBZ = zeros(length(kxx),1);
ksumx_110All = zeros(length(kxx),1);
ksumx_111SBZ = zeros(length(kxx),1);
ksumx_111All = zeros(length(kxx),1);
np_110All = 0; 
np_111All = 0; 
np_110SBZ = 0; 
np_111SBZ = 0; 
np_corner = 0;
tao_1_sum = zeros(length(kxx),1);
tao_INV_1_sum = zeros(length(kxx),1);

v2tao_1_sum = zeros(length(kxx),1);

v2_1_sum = zeros(length(kxx),1);
v2_2_sum = zeros(length(kxx),1);


Cpsumnk = zeros(length(kxx),1);

ksum = 0;

disp(nbranch)
ninc = 0; 
Pendflag = 1:Ptotal; %
nk_MAX =length(kxx);
% ww = ww.ww;
% vv = vv.vv;
  for Pendflag1 = 1:Ptotal; % 1:Ptotal; % 75;
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
    fname4 = sprintf('Ave_k_P%d_T%d.csv', Pn(end),Ptotal);
    
    v2tao_1_bin = zeros(bandwidth,nbins);
v2_bin = zeros(bandwidth,nbins);
tao_1_bin = zeros(bandwidth,nbins);
    k_Gillet_bin = zeros(bandwidth,nbins);
for nk = kstart:kstop; % 1: length(kxx); % 

   if(mod(nk,2000)==0)
   	disp(nk/length(kxx)) 
   end
   
    kx=kxx(nk);
    ky=kyy(nk); 
    kz=kzz(nk); 
        
    % % % Direct Summation in Sphere BZ.
    
    ninc = ninc + 1;
    wwmodec = ww(:,nk);
    vmode = vv(:,nk);

      wwmode = wwmodec; % wwI*sqrt(1/massSi); 
      [Nbin,edges,binIdx] = histcounts(wwmode, [binEdges(1:end-1) Inf]);
        Cpmode = kB*(h*wwmode/kB/T).^2.*(exp(h*wwmode/kB/T)./((exp(h*wwmode/kB/T)-1).^2));
        Cpsum = Cpsum + sum(sum((Cpmode))); 
%         Cpsumnk(nk) = sum(sum((Cpmode)));
% 		   Cpmodec = kB*(h*wwmodec/kB/T).^2.*(exp(h*wwmodec/kB/T)./((exp(h*wwmodec/kB/T)-1).^2));
        Cpmodec = Cpmode; 
% 		   Cpsumnkc(nk) = sum(sum((Cpmodec)));

        
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
           if(kz == ZZ(1) && kx == ky) % kx == ky) % 
               np_110All = np_110All + 1; 
               ksumx_110All(nk) = ksumx_4Gilletc(nk); 
               kpos_110All(np_110All) = nk;
           end
           if(kz == kx && kx == ky)  
                ksumx_111All(nk) = ksumx_4Gilletc(nk); 
                np_111All = np_111All + 1;                 
                kpos_111All(np_111All) = nk;
           end
           if(sqrt(kx^2+ky^2+kz^2)>(XX(end)))
%                continue;
               np_corner = np_corner + 1; 
               ksumx_110SBZ(nk) = 0; % ksumx_4Gilletc(nk);
               ksumx_111SBZ(nk) = 0; % ksumx_4Gilletc(nk); 
               
           else
            ksumx_SBZ(nk) = ksumx_4Gilletc(nk);
                if(kz == ZZ(1) && kx == ky) % kx == ky) % 
                    ksumx_110SBZ(nk) = ksumx_4Gilletc(nk); 
                    np_110SBZ = np_110SBZ + 1;
                    kpos_110SBZ(np_110SBZ) = nk;
%                     ksumx_110All(nk) = ksumx_4Gilletc(nk);                    
                end
                if(kz == kx && kx == ky)
                    ksumx_111SBZ(nk) = ksumx_4Gilletc(nk); 
                    np_111SBZ = np_111SBZ + 1;
                    kpos_111SBZ(np_111SBZ) = nk;                    
%                     ksumx_111All(nk) = ksumx_4Gilletc(nk); 
                end
           end
% % %            Taosumx_Uprocess(nk) = sum(sum((1./Uprocess_scat_invc)))*8/V;
% % %            Taosumx_inc(nk) = sum(sum((1./inc_scat_inv)))*8/V;
% % %            Taosumx_UprocINV(nk) = sum(sum((Uprocess_scat_invc)))*8/V;
% % %            Taosumx_incINV(nk) = sum(sum((inc_scat_inv)))*8/V;
           
           
           tao_1 = tot_scat_Gilletc; 
% % %            tao_INV_1 = tot_scat_inv_Gilletc; 
           tao_3 = 1; 
           
% % %     if(sqrt(kx^2+ky^2+kz^2)>XX(end))
% % %         tao_2 = 0; 
% % %         tao_2_edge = tot_scat_Gilletc;
% % %         tao_4 = 0; 
% % %     else
% % %         ncount = ncount+1;
% % %         tao_2 = tot_scat_Gilletc;
% % %         tao_2_edge = 0; 
% % %         tao_4 = 1; 
% % %     end
    
    tao_1_sum(nk) = sum(tao_1); % sqrt(tao_1'*tao_1);
    v2tao_1_temp = vmode.*vmode.*tao_1;
    v2tao_1_sum(nk) = sum(v2tao_1_temp); % sum(vmode.*vmode.*tao_1);
    v2_temp = vmode.*vmode; 
    tao_1_temp = tao_1; 
               k_Gillet_temp = kmode_holland4_Gilletc; % *8/V;
%        for i = 1:nbranch; 
%         v2tao_1_bin(nk-(kstart-1),binIdx(1:nbranch)) = v2tao_1_temp(1:nbranch,1);
%         v2_bin(nk-(kstart-1),binIdx(1:nbranch)) = v2_temp(1:nbranch,1)';
%         tao_1_bin(nk-(kstart-1),binIdx(1:nbranch)) = tao_1_temp(1:nbranch,1)';
%         k_Gillet_bin(nk-(kstart-1),binIdx(1:nbranch)) = kmode_holland4_Gilletc(1:nbranch,1)';
for i = 1:nbranch; 
        v2tao_1_bin(nk-(kstart-1),binIdx(i)) = v2tao_1_bin(nk-(kstart-1),binIdx(i)) + v2tao_1_temp(i,1);
        v2_bin(nk-(kstart-1),binIdx(i)) = v2_bin(nk-(kstart-1),binIdx(i)) + v2_temp(i,1);
        tao_1_bin(nk-(kstart-1),binIdx(i)) = tao_1_bin(nk-(kstart-1),binIdx(i)) + tao_1_temp(i,1);
        k_Gillet_bin(nk-(kstart-1),binIdx(i)) = k_Gillet_bin(nk-(kstart-1),binIdx(i))+ kmode_holland4_Gilletc(i,1);
end
%     fname1 = sprintf
end
% % csvwrite(fname1,v2tao_1_bin);
% % csvwrite(fname2,v2_bin);
% % csvwrite(fname3,tao_1_bin);
% % csvwrite(fname4,k_Gillet_bin);

  end
  k_Gilletc = sum(sum(ksumx_4Gilletc))
k_SBZ = sum(sum(ksumx_SBZ))
k_110SBZ = sum(sum(ksumx_110SBZ))
k_110All = sum(sum(ksumx_110All))
k_111SBZ = sum(sum(ksumx_111SBZ))
k_111All = sum(sum(ksumx_111All))




% % % disp('End of these current parts')
% % % pause; 
% % % Read the data into file; 
v2tao_1_bin2 = v2tao_1_bin; % zeros(length(kxx),nbins);
v2_bin2 = v2_bin; % zeros(length(kxx),nbins);
tao_1_bin2 = tao_1_bin; % zeros(length(kxx),nbins); 
k_Gillet_bin2 = k_Gillet_bin; % zeros(length(kxx),nbins); 

%     v2tao_1_bin2 = zeros(length(kxx),nbins);
% v2_bin2 = zeros(length(kxx),nbins);
% tao_1_bin2 = zeros(length(kxx),nbins); 
% k_Gillet_bin2 = zeros(length(kxx),nbins); 
% 
% for Pendflag1 = 1:Ptotal; % 1:Ptotal; % 75;
%         Pn = Pendflag(Pendflag1);
%     disp(Pn)
%     kindexendflag = 0; 
%     bandwidth = Psize2; % length(kxx)/Ptotal;
%     kstart = (1)+(Pn(end)-1)*(bandwidth); 
%     kstop = (bandwidth)+(Pn(end)-1)*bandwidth;
%     if(kstart > length(kxx))
%        fprintf('End of all 3rdIFC pairs'); 
%        kindexendflag = 1; 
%        return;
%     else if(kstop > nk_MAX);
%             kstop = nk_MAX;
%         end
%     end
% 
%     fname1 = sprintf('Ave_v2tao_P%d_T%d.csv', Pn(end),Ptotal);
%     fname2 = sprintf('Ave_v2_P%d_T%d.csv', Pn(end),Ptotal);
%     fname3 = sprintf('Ave_tao_P%d_T%d.csv', Pn(end),Ptotal);
%     fname4 = sprintf('Ave_k_P%d_T%d.csv', Pn(end),Ptotal);
%     clear v2tao_temp v2_temp tao_temp k_Gillet_temp
% 
%         v2tao_1_bin2(kstart:kstop,:) = load(fname1,'-ascii'); % v2tao_temp(1:bandwidth,:); %  v2tao_1_temp(1:nbranch,1);
%         v2_bin2(kstart:kstop,:) = load(fname2,'-ascii'); % v2_temp(1:bandwidth,:);
%         tao_1_bin2(kstart:kstop,:) = load(fname3,'-ascii'); % tao_1_temp(1:bandwidth,:);
%          k_Gillet_bin2(kstart:kstop,:) = load(fname4,'-ascii'); % tao_1_temp(1:bandwidth,:);
% 
% end


Ave_v2_bin = sum(v2_bin2,1)*8/V; 
Ave_v2tao_1_bin = sum(v2tao_1_bin2,1)*8/V;
Ave_tao_1_bin = sum(tao_1_bin2,1)*8/V;
Ave_k_Gillet_bin = sum(k_Gillet_bin2)*8/V;
save('Ave_v2_bin.mat','Ave_v2_bin');
save('Ave_v2tao_1_bin.mat','Ave_v2tao_1_bin');
save('Ave_tao_1_bin.mat','Ave_tao_1_bin');
save('Ave_k_Gillet_bin.mat','Ave_k_Gillet_bin');

Ave_Tao_FullBZ = sum(tao_1_sum)*8/V % sum(tao_1_sum)/nbranch/length(kxx) % ./nbranch)/size(tao_1_sum,1)
Ave_Tao_INV_FullBZ = sum(tao_INV_1_sum)*8/V % sum(tao_1_sum)/nbranch/length(kxx) % ./nbranch)/size(tao_1_sum,1)

Ave_v2Tao_FullBZ = sum(v2tao_1_sum)*8/V % sum(v2tao_1_sum)/nbranch/length(kxx);

Ave_v2_FullBZ = sum(v2_1_sum)*8/V % sum(v2tao_1_sum)/nbranch/length(kxx);
Ave_v2_FullBZ2 = sum(v2_2_sum)*8/V % sum(v2tao_2_sum)/nbranch/ncount;
k_Gilletc = sum(sum(ksumx_4Gilletc))
k_SBZ = sum(sum(ksum_SBZ))
sum(sum(k_Gillet_bin))*8/V

k_FBZ = sum(sum(Ave_k_Gillet_bin))
v2tao_FBZ = sum(sum(Ave_v2tao_1_bin))
v2_FBZ = sum(sum(Ave_v2_bin))
tao_FBZ = sum(sum(Ave_tao_1_bin))
[k_Gilletc, K_FBZ,v2tao_FBZ,v2_FBZ,tao_FBZ]












% figure
aj = binEdges(1:end-1);     %# bins lower edge
bj = binEdges(2:end);       %# bins upper edge
cj = ( aj + bj ) ./ 2;      %# bins center

%# plot histogram
figure(1)
% bar(cj,Ave_v2tao_1_bin,'hist')
plot(cj,Ave_v2tao_1_bin,'k-')
hold on
% set(gca, 'XTick',binEdges, 'XLim',[binEdges(1):binEdges(end)])
% xlabel('Bins'), ylabel('Counts'), title('histogram of v2tao')
figure(3);
hold on;
plot(cj, Ave_v2_bin, 'k-');
figure(5)
hold on;
plot(cj,Ave_tao_1_bin,'k-');
hold on;


% * *  *  * * * * ----------------------
% Read data from Folder : S:\S-current-debug-folder\C-H_mode_08192015\CHModel_Primitive\1stPaper-k-frequency-dependent-CorrectUse-6pi2-2016-12-21
k_IsoSi2 = load('k_IsoSi2.mat'); 
k_IsoSi3 = load('k_IsoSi3.mat'); 
k_IsoSi4 = load('k_IsoSi4.mat'); 
k_IsoSi5 = load('k_IsoSi5.mat'); 

Tao_IsoSi2 = load('Tao_IsoSi2.mat'); 
Tao_IsoSi3 = load('Tao_IsoSi3.mat'); 
Tao_IsoSi4 = load('Tao_IsoSi4.mat'); 
Tao_IsoSi5 = load('Tao_IsoSi5.mat'); 


v2_IsoSi2 = load('v2_IsoSi2.mat'); 
v2_IsoSi3 = load('v2_IsoSi3.mat'); 
v2_IsoSi4 = load('v2_IsoSi4.mat'); 
v2_IsoSi5 = load('v2_IsoSi5.mat'); 

v2tao_IsoSi2 = load('v2tao_IsoSi2.mat'); 
v2tao_IsoSi3 = load('v2tao_IsoSi3.mat'); 
v2tao_IsoSi4 = load('v2tao_IsoSi4.mat'); 
v2tao_IsoSi5 = load('v2tao_IsoSi5.mat'); 

k_IsoSi2 = k_IsoSi2.Ave_kGillet_bin; 
k_IsoSi3 = k_IsoSi3.Ave_kGillet_bin; 
k_IsoSi4 = k_IsoSi4.Ave_kGillet_bin; 
k_IsoSi5 = k_IsoSi5.Ave_kGillet_bin; 

Tao_IsoSi2 = Tao_IsoSi2.Ave_tao_1_bin;
Tao_IsoSi3 = Tao_IsoSi3.Ave_tao_1_bin;
Tao_IsoSi4 = Tao_IsoSi4.Ave_tao_1_bin;
Tao_IsoSi5 = Tao_IsoSi5.Ave_tao_1_bin;

v2_IsoSi2 = v2_IsoSi2.Ave_v2_bin;
v2_IsoSi3 = v2_IsoSi3.Ave_v2_bin;
v2_IsoSi4 = v2_IsoSi4.Ave_v2_bin;
v2_IsoSi5 = v2_IsoSi5.Ave_v2_bin;

v2tao_IsoSi2 = v2tao_IsoSi2.Ave_v2tao_1_bin;
v2tao_IsoSi3 = v2tao_IsoSi3.Ave_v2tao_1_bin;
v2tao_IsoSi4 = v2tao_IsoSi4.Ave_v2tao_1_bin;
v2tao_IsoSi5 = v2tao_IsoSi5.Ave_v2tao_1_bin;

    % % % FOR BIN Calculation; 
    
nbins = 200;
binEdges = linspace(5.0e10,1.12e+14,nbins+1);

% figure
aj = binEdges(1:end-1);     %# bins lower edge
bj = binEdges(2:end);       %# bins upper edge
cj = ( aj + bj ) ./ 2;      %# bins center

%# plot histogram
figure(1)
% bar(cj,Ave_v2tao_1_bin,'hist')
plot(cj,Ave_v2tao_1_bin,'k-')
hold on
plot(cj,Iso_Si5,'r--')
legend('FullBZ', 'Iso')
title('v2tao-Si5Ge1')

figure(2)
plot(cj,Ave_v2tao_1_bin-Iso_Si5,'b-.')
title('Difference(FBZ-Iso)-v2tao-Si5Ge1')
figure(4)
plot(cj,(Ave_v2tao_1_bin-Iso_Si5)./Ave_v2tao_1_bin,'b-.')
title('Ratio: (FBZ-Iso)/FBZ -v2tao-Si5Ge1')

% % figure(6)
% % plot(cj,(Iso_Si5)./Ave_v2tao_1_bin,'g-.')
% % title('Ratio: (FBZ-Iso)/FBZ -v2tao-Si5Ge1')


figure(2)
plot(cj,(FBZ_Si2-Iso_Si2)./FBZ_Si2,'k-','LineWidth',2)
hold on
plot(cj,(FBZ_Si3-Iso_Si3)./FBZ_Si3,'b--','LineWidth',2)
hold on
plot(cj,(FBZ_Si4-Iso_Si4)./FBZ_Si4,'m-.','LineWidth',2)
hold on
plot(cj,(FBZ_Si5-Iso_Si5)./FBZ_Si5,'g-','LineWidth',2)
title('Difference Ratio(FBZ-Iso)/FBZ-v2tao-All sample')
legend('Si2Ge1', 'Si3Ge1', 'Si4Ge1', 'Si5Ge1')

D2 = (FBZ_Si2-Iso_Si2)./FBZ_Si2; % FBZ_Si2-Iso_Si2;
D3 = (FBZ_Si3-Iso_Si3)./FBZ_Si3; % FBZ_Si3-Iso_Si3;
D4 = (FBZ_Si4-Iso_Si4)./FBZ_Si4; % FBZ_Si4-Iso_Si4;
D5 = (FBZ_Si5-Iso_Si5)./FBZ_Si5; % FBZ_Si5-Iso_Si5;

D2 = (FBZ_Si2-Iso_Si2)./sum(FBZ_Si2); % FBZ_Si2-Iso_Si2;
D3 = (FBZ_Si3-Iso_Si3)./sum(FBZ_Si3); % FBZ_Si3-Iso_Si3;
D4 = (FBZ_Si4-Iso_Si4)./sum(FBZ_Si4); % FBZ_Si4-Iso_Si4;
D5 = (FBZ_Si5-Iso_Si5)./sum(FBZ_Si5); % FBZ_Si5-Iso_Si5;

for i=1:length(D2);
    AD2(i) = sum(abs(D2(1:i)));
    AD3(i) = sum(abs(D3(1:i)));
    AD4(i) = sum(abs(D4(1:i)));
    AD5(i) = sum(abs(D5(1:i)));
end

for i=1:length(D2);
    AD2(i) = sum((D2(1:i)));
    AD3(i) = sum((D3(1:i)));
    AD4(i) = sum((D4(1:i)));
    AD5(i) = sum((D5(1:i)));
end

figure(4)
plot(cj,AD2,'k-','LineWidth',2)
hold on
plot(cj,AD3,'b--','LineWidth',2)
hold on
plot(cj,AD4,'m-.','LineWidth',2)
hold on
plot(cj,AD5,'g-','LineWidth',2)
title('Accumulative Absolute Diff Ratio.-v2tao-All sample')
legend('Si2Ge1', 'Si3Ge1', 'Si4Ge1', 'Si5Ge1')

% % % 5 - - 12/21/2016 Figure Plots:
D2 = (k_FBZSi2-k_IsoSi2)./sum(k_FBZSi2); % FBZ_Si2-Iso_Si2;
D3 = (k_FBZSi3-k_IsoSi3)./sum(k_FBZSi3); % FBZ_Si2-Iso_Si2;
D4 = (k_FBZSi4-k_IsoSi4)./sum(k_FBZSi4); % FBZ_Si2-Iso_Si2;
D5 = (k_FBZSi5-k_IsoSi5)./sum(k_FBZSi5); % FBZ_Si2-Iso_Si2;

for i=1:length(D2);
    AD2(i) = sum((D2(1:i)));
    AD3(i) = sum((D3(1:i)));
    AD4(i) = sum((D4(1:i)));
    AD5(i) = sum((D5(1:i)));
end

for i=1:length(D2);
    BD2(i) = sum(abs(D2(1:i)));
    BD3(i) = sum(abs(D3(1:i)));
    BD4(i) = sum(abs(D4(1:i)));
    BD5(i) = sum(abs(D5(1:i)));
end

figure(24)
plot(cj,D2,'k-','LineWidth',2)
hold on
plot(cj,D3,'b--','LineWidth',2)
hold on
plot(cj,D4,'m-.','LineWidth',2)
hold on
plot(cj,D5,'g-','LineWidth',2)
title('Frequency dependent.-k-All sample')
legend('Si2Ge1', 'Si3Ge1', 'Si4Ge1', 'Si5Ge1')

figure(25)
plot(cj,AD2,'k-','LineWidth',2)
hold on
plot(cj,AD3,'b--','LineWidth',2)
hold on
plot(cj,AD4,'m-.','LineWidth',2)
hold on
plot(cj,AD5,'g-','LineWidth',2)
title('Accumulative  Diff Ratio.-k-All sample')
legend('Si2Ge1', 'Si3Ge1', 'Si4Ge1', 'Si5Ge1')

figure(26)
plot(cj,BD2,'k-','LineWidth',2)
hold on
plot(cj,BD3,'b--','LineWidth',2)
hold on
plot(cj,BD4,'m-.','LineWidth',2)
hold on
plot(cj,BD5,'g-','LineWidth',2)
title('Accumulative Absolute Diff Ratio.-k-All sample')
legend('Si2Ge1', 'Si3Ge1', 'Si4Ge1', 'Si5Ge1')


