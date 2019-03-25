% Data from: C:\ww_vv_Calculation_GeXSi1\ww_vv_BulkSi3_Nk80\BulkSi3_Paper1st_Frequency_deendAnalysis_Nparts_12_08_2016.m

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
title('BulkSi3 - 110SBZ(black star) vs 110All(red circle)')
saveas(gcf,'k-110-SBZvsAll-BulkSi3.bmp')
saveas(gcf,'k-110-SBZvsAll-BulkSi3.fig')

figure(2)
dist_111SBZ = sqrt(kxx(kpos_111SBZ(1:np_111SBZ)).^2 + kyy(kpos_111SBZ(1:np_111SBZ)).^2 + kzz(kpos_111SBZ(1:np_111SBZ)).^2)
semilogy(dist_111SBZ, ksumx_111SBZ(kpos_111SBZ(1:np_111SBZ)), 'b*')
hold on
dist_111All = sqrt(kxx(kpos_111All(1:np_111All)).^2 + kyy(kpos_111All(1:np_111All)).^2 + kzz(kpos_111All(1:np_111All)).^2)
semilogy(dist_111All, ksumx_111All(kpos_111All(1:np_111All)), 'ro')
ylabel('Thermal conductivity (W/mK)')
xlabel('k distantce to origin')
title('BulkSi3 - 111SBZ(blue star) vs 111All(red circle)')
saveas(gcf,'k-111-SBZvsAll-BulkSi3.bmp')
saveas(gcf,'k-111-SBZvsAll-BulkSi3.fig')

figure(3)
plot(dist_110SBZ, ww(:,kpos_110SBZ(1:np_110SBZ)), 'k*')
hold on
plot(dist_110All, ww(:,kpos_110All(1:np_110All)), 'ro')
xlabel('k distantce to origin')
title('BulkSi3 - 110SBZ(black star) vs 110All(red circle)')
saveas(gcf,'Dispersion-110-SBZvsAll-BulkSi3.bmp')
saveas(gcf,'Dispersion-110-SBZvsAll-BulkSi3.fig')

figure(4)
plot(dist_111SBZ, ww(:,kpos_111SBZ(1:np_111SBZ)), 'b*')
hold on
plot(dist_111All, ww(:,kpos_111All(1:np_111All)), 'ro')
xlabel('k distantce to origin')
title('BulkSi3 - 111SBZ(blue star) vs 111All(red circle)')
saveas(gcf,'Dispersion-111-SBZvsAll-BulkSi3.bmp')
saveas(gcf,'Dispersion-111-SBZvsAll-BulkSi3.fig')