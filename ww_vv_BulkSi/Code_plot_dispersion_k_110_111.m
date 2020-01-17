% Data from: C:\ww_vv_Calculation_GeXSi1\ww_vv_BulkSi5_Nk80\BulkSi5_Paper1st_Frequency_deendAnalysis_Nparts_12_08_2016.m

sum(ksumx_111SBZ(kpos_111SBZ(1:np_111SBZ)))
sum(ksumx_111All(kpos_111All(1:np_111All)))
sum(ksumx_110All(kpos_110All(1:np_110All)))
sum(ksumx_110SBZ(kpos_110SBZ(1:np_110SBZ)))

figure(1)
dist_110SBZ = sqrt(kxx(kpos_110SBZ(1:np_110SBZ)).^2 + kyy(kpos_110SBZ(1:np_110SBZ)).^2)
semilogy(dist_110SBZ, ksumx_110SBZ(kpos_110SBZ(1:np_110SBZ)), 'k*')
hold on
dist_110All = sqrt(kxx(kpos_110All(1:np_110All)).^2 + kyy(kpos_110All(1:np_110All)).^2)
semilogy(dist_110All, ksumx_110All(kpos_110All(1:np_110All)), 'ro')
ylabel('Thermal conductivity (W/mK)')
xlabel('k distantce to origin')
title('BulkSi5 - 110SBZ(black star) vs 110All(red circle)')
saveas(gcf,'k-110-SBZvsAll-BulkSi5.bmp')
saveas(gcf,'k-110-SBZvsAll-BulkSi5.fig')

figure(2)
dist_111SBZ = sqrt(kxx(kpos_111SBZ(1:np_111SBZ)).^2 + kyy(kpos_111SBZ(1:np_111SBZ)).^2 + kzz(kpos_111SBZ(1:np_111SBZ)).^2)
semilogy(dist_111SBZ, ksumx_111SBZ(kpos_111SBZ(1:np_111SBZ)), 'b*')
hold on
dist_111All = sqrt(kxx(kpos_111All(1:np_111All)).^2 + kyy(kpos_111All(1:np_111All)).^2 + kzz(kpos_111All(1:np_111All)).^2)
semilogy(dist_111All, ksumx_111All(kpos_111All(1:np_111All)), 'ro')
ylabel('Thermal conductivity (W/mK)')
xlabel('k distantce to origin')
title('BulkSi5 - 111SBZ(blue star) vs 111All(red circle)')
saveas(gcf,'k-111-SBZvsAll-BulkSi5.bmp')
saveas(gcf,'k-111-SBZvsAll-BulkSi5.fig')

figure(3)
plot(dist_110SBZ, ww(:,kpos_110SBZ(1:np_110SBZ)), 'k*')
hold on
plot(dist_110All, ww(:,kpos_110All(1:np_110All)), 'ro')
xlabel('k distantce to origin')
title('BulkSi5 - 110SBZ(black star) vs 110All(red circle)')
saveas(gcf,'Dispersion-110-SBZvsAll-BulkSi5.bmp')
saveas(gcf,'Dispersion-110-SBZvsAll-BulkSi5.fig')

figure(4)
plot(dist_111SBZ, ww(:,kpos_111SBZ(1:np_111SBZ)), 'b*')
hold on
plot(dist_111All, ww(:,kpos_111All(1:np_111All)), 'ro')
xlabel('k distantce to origin')
title('BulkSi5 - 111SBZ(blue star) vs 111All(red circle)')
saveas(gcf,'Dispersion-111-SBZvsAll-BulkSi5.bmp')
saveas(gcf,'Dispersion-111-SBZvsAll-BulkSi5.fig')

figure(5)
semilogy(dist_110SBZ, tao_1_sum(kpos_110SBZ(1:np_110SBZ)), 'k*')
hold on
semilogy(dist_110All, tao_1_sum(kpos_110All(1:np_110All)), 'ro')
xlabel('k distantce to origin')
title('BulkSi5 - 110SBZ(black star) vs 110All(red circle) - Tao')
saveas(gcf,'Tao-110-SBZvsAll-BulkSi5.bmp')
saveas(gcf,'Tao-110-SBZvsAll-BulkSi5.fig')

figure(6)
semilogy(dist_111SBZ, tao_1_sum(kpos_111SBZ(1:np_111SBZ)), 'b*')
hold on
semilogy(dist_111All, tao_1_sum(kpos_111All(1:np_111All)), 'ro')
xlabel('k distantce to origin')
title('BulkSi5 - 111SBZ(blue star) vs 111All(red circle) - Tao')
saveas(gcf,'Tao-111-SBZvsAll-BulkSi5.bmp')
saveas(gcf,'Tao-111-SBZvsAll-BulkSi5.fig')

figure(7)
plot(dist_110SBZ, v2_temp2(kpos_110SBZ(1:np_110SBZ)), 'k*')
hold on
plot(dist_110All, v2_temp2(kpos_110All(1:np_110All)), 'ro')
xlabel('k distantce to origin')
title('BulkSi5 - 110SBZ(black star) vs 110All(red circle) - v2')
saveas(gcf,'v2-110-SBZvsAll-BulkSi5.bmp')
saveas(gcf,'v2-110-SBZvsAll-BulkSi5.fig')

figure(8)
plot(dist_111SBZ, v2_temp2(kpos_111SBZ(1:np_111SBZ)), 'b*')
hold on
plot(dist_111All, v2_temp2(kpos_111All(1:np_111All)), 'ro')
xlabel('k distantce to origin')
title('BulkSi5 - 111SBZ(blue star) vs 111All(red circle) - v2')
saveas(gcf,'v2-111-SBZvsAll-BulkSi5.bmp')
saveas(gcf,'v2-111-SBZvsAll-BulkSi5.fig')

figure(9)
plot(dist_110SBZ, Cpsumnk(kpos_110SBZ(1:np_110SBZ)), 'k*')
hold on
plot(dist_110All, Cpsumnk(kpos_110All(1:np_110All)), 'ro')
xlabel('k distantce to origin')
title('BulkSi5 - 110SBZ(black star) vs 110All(red circle) - Cp')
saveas(gcf,'Cp-110-SBZvsAll-BulkSi5.bmp')
saveas(gcf,'Cp-110-SBZvsAll-BulkSi5.fig')

figure(10)
plot(dist_111SBZ, Cpsumnk(kpos_111SBZ(1:np_111SBZ)), 'b*')
hold on
plot(dist_111All, Cpsumnk(kpos_111All(1:np_111All)), 'ro')
xlabel('k distantce to origin')
title('BulkSi5 - 111SBZ(blue star) vs 111All(red circle) - Cp')
saveas(gcf,'Cp-111-SBZvsAll-BulkSi5.bmp')
saveas(gcf,'Cp-111-SBZvsAll-BulkSi5.fig')

Para_110_SBZ  = [Cpsumnk(kpos_110SBZ(1:np_110SBZ)), v2_temp2(kpos_110SBZ(1:np_110SBZ)), tao_1_sum(kpos_110SBZ(1:np_110SBZ))]
Para_110_All = [Cpsumnk(kpos_110All(1:np_110All)), v2_temp2(kpos_110All(1:np_110All)), tao_1_sum(kpos_110All(1:np_110All))]; 
Para_111_SBZ = [Cpsumnk(kpos_111SBZ(1:np_111SBZ)), v2_temp2(kpos_111SBZ(1:np_111SBZ)), tao_1_sum(kpos_111SBZ(1:np_111SBZ))] ; % Cp, v2, tao
Para_111_All = [Cpsumnk(kpos_111All(1:np_111All)), v2_temp2(kpos_111All(1:np_111All)), tao_1_sum(kpos_111All(1:np_111All))];
if(sum(kpos_110SBZ(1:np_110SBZ) - kpos_110All(1:np_110SBZ)) ~=0)
    disp('Error! Ratio calcuation are wrong')
    pause;
end
if(sum(kpos_111SBZ(1:np_111SBZ) - kpos_111All(1:np_111SBZ)) ~=0)
    disp('Error! Ratio calcuation are wrong')
    pause;
end
Corner_110 = kpos_110All(1+np_110SBZ:np_110All); %^ setdiff(kpos_110SBZ(1:np_110SBZ),kpos_110All(1:np_110All));
Corner_111 = kpos_111All(1+np_111SBZ:np_110All); %^ setdiff(kpos_111SBZ(1:np_111SBZ),kpos_111All(1:np_111All));
Para_110_Corner = [Cpsumnk(Corner_110), v2_temp2(Corner_110), tao_1_sum(Corner_110)];
Para_111_Corner = [Cpsumnk(Corner_111), v2_temp2(Corner_111), tao_1_sum(Corner_111)];

Mean_110_All = mean(Para_110_All)
Mean_111_All = mean(Para_111_All)
Mean_110_Corner = mean(Para_110_Corner)
Mean_111_Corner = mean(Para_111_Corner)
[Mean_110_Corner, Mean_110_All, Mean_110_Corner./Mean_110_All, Mean_111_Corner, Mean_111_All, Mean_111_Corner./Mean_111_All]

% -- 2018/01/15 -- Add Frequency on Fig. 10(a)
wwsumnk = sum(ww);
%ww2 = ww.*ww; 
% w4sumnk = sum(ww.*ww.*);
ww_110_Corner  = mean(wwsumnk(Corner_110)); % , v2_temp2(kpos_110SBZ(1:np_110SBZ)), tao_1_sum(kpos_110SBZ(1:np_110SBZ))]
ww_110_All = mean(wwsumnk(kpos_110All(1:np_110All))); % , v2_temp2(kpos_110All(1:np_110All)), tao_1_sum(kpos_110All(1:np_110All))]; 
ww_110_Corner./ww_110_All

clear wwsumnk; 
w2sumnk = sum(ww.*ww);
w2_110_Corner  = mean(w2sumnk(Corner_110)); % 
w2_110_All = mean(w2sumnk(kpos_110All(1:np_110All))); %  
w2_110_Corner./w2_110_All

clear w2sumnk; 
w4sumnk = sum(ww.*ww.*ww.*ww);
w4_110_Corner  = mean(w4sumnk(Corner_110)); % 
w4_110_All = mean(w4sumnk(kpos_110All(1:np_110All))); %  
w4_110_Corner./w4_110_All

clear TauUI_110_Corner TauUI_110_All
TauUI_110_Corner = mean([tao_umk_sum(Corner_110), tao_I_sum(Corner_110), tao_Uprocess_sum(Corner_110), wmode_sum(Corner_110), w2mode_sum(Corner_110), w4mode_sum(Corner_110)]);
TauUI_110_All = mean([tao_umk_sum(kpos_110All(1:np_110All)), tao_I_sum(kpos_110All(1:np_110All)), tao_Uprocess_sum(kpos_110All(1:np_110All)), wmode_sum(kpos_110All(1:np_110All)), w2mode_sum(kpos_110All(1:np_110All)), w4mode_sum(kpos_110All(1:np_110All))]); 
TauUI_110_Corner./TauUI_110_All

% Calculate 1/w, 1/w2, 1/w4 ------------
wwsumnk = sum(1./ww);
ww_110_Corner  = mean(wwsumnk(Corner_110)); % , v2_temp2(kpos_110SBZ(1:np_110SBZ)), tao_1_sum(kpos_110SBZ(1:np_110SBZ))]
ww_110_All = mean(wwsumnk(kpos_110All(1:np_110All))); % , v2_temp2(kpos_110All(1:np_110All)), tao_1_sum(kpos_110All(1:np_110All))]; 
ww_110_Corner./ww_110_All

clear wwsumnk; 
w2sumnk = sum(1./(ww.*ww));
w2_110_Corner  = mean(w2sumnk(Corner_110)); % 
w2_110_All = mean(w2sumnk(kpos_110All(1:np_110All))); %  
w2_110_Corner./w2_110_All

clear w2sumnk; 
w4sumnk = sum(1./(ww.*ww.*ww.*ww));
w4_110_Corner  = mean(w4sumnk(Corner_110)); % 
w4_110_All = mean(w4sumnk(kpos_110All(1:np_110All))); %  
w4_110_Corner./w4_110_All