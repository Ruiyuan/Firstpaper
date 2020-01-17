% function [xs ys zs,nsupercell]=creatSuperCellwVoid(A, ncx, ncy, ncz,ncell,nbasis,ncxsuper,ncxvoid,flagsupercell)
clear xs ys zs; 
n1=0;
pvoid = -1;
if(ncxvoid == 1)
    pvoid = median(0:1:ncxsuper-1); 
    xvoid_lower = (ncx/2 - ncxvoid/2 -0.1)*A;
    xvoid_upper = (ncx/2 + ncxvoid/2 -0.1)*A;
    yvoid_lower = (ncy/2 - ncxvoid/2 -0.1)*A;
    yvoid_upper = (ncy/2 + ncxvoid/2 -0.1)*A;
    zvoid_lower = (ncz/2 - ncxvoid/2 -0.1)*A;
    zvoid_upper = (ncz/2 + ncxvoid/2 -0.1)*A;
end
for i=0:ncxsuper-1;
    for j=0:ncxsuper-1;
        for k=0:ncxsuper-1;
% %             if (i==pvoid & j==pvoid & k==pvoid)
% %                 continue;
% %             end
            for n0=1:8;
                n1=n1+1;
                if(rem(n0,8)==1)
                    xs(n1)=i*A;
                    ys(n1)=j*A;
                    zs(n1)=k*A;
                else
                    if(rem(n0,8)==2)
                    xs(n1)=xs(n1-1)+A/4;
                    ys(n1)=ys(n1-1)+A/4;
                    zs(n1)=zs(n1-1)+A/4;
                    else
                    if(rem(n0,8)==3)
                    xs(n1)=i*A+A/2;
                    ys(n1)=j*A+A/2;
                    zs(n1)=k*A;
                    else
                    if(rem(n0,8)==4)
                    xs(n1)=xs(n1-1)+A/4;
                    ys(n1)=ys(n1-1)+A/4;
                    zs(n1)=zs(n1-1)+A/4;
                    else
                    if(rem(n0,8)==5)
                    xs(n1)=i*A;
                    ys(n1)=j*A+A/2;
                    zs(n1)=k*A+A/2;
                    else
                    if(rem(n0,8)==6)
                    xs(n1)=xs(n1-1)+A/4;
                    ys(n1)=ys(n1-1)+A/4;
                    zs(n1)=zs(n1-1)+A/4;
                    else
                    if(rem(n0,8)==7)
                    xs(n1)=i*A+A/2;
                    ys(n1)=j*A;
                    zs(n1)=k*A+A/2;
                    else
                    if(rem(n0,8)==0)
                    xs(n1)=xs(n1-1)+A/4;
                    ys(n1)=ys(n1-1)+A/4;
                    zs(n1)=zs(n1-1)+A/4;
                    end
                    end
                    end
                    end
        
                 end
                 end
                 end
                end

            end
        end
    end
end

n1 = 0; 
for i=1:length(xs)
                        if((xvoid_lower<xs(i))&&(xvoid_upper>xs(i))&&(yvoid_lower<ys(i))&&(yvoid_upper>ys(i))&&(zvoid_lower<zs(i))&&(zvoid_upper>zs(i)))
                            continue; 
                        else 
                            n1=n1+1;
                            xnew(n1) = xs(i);
                            ynew(n1) = ys(i);
                            znew(n1) = zs(i);
                        end
end

xs1=xs; ys1=ys; zs1=zs;
clear xs ys zs;
xs=xnew; ys=ynew; zs=znew; 

Nx = ncx/ncxsuper;
Ny = ncy/ncxsuper;
Nz = ncz/ncxsuper;
if (flagsupercell == 1)
    nsupercell = (ncxsuper^3-ncxvoid)*8;
else
    nsupercell = ncell;
end
Nsc1 = 0;
for i=0:Nx-1; 
    for j=0:Ny-1;
        for k=0:Nz-1;            
            xs((1:nsupercell)+Nsc1*nsupercell) = xs(1:nsupercell)+i*A*ncxsuper;
            ys((1:nsupercell)+Nsc1*nsupercell) = ys(1:nsupercell)+j*A*ncxsuper;
            zs((1:nsupercell)+Nsc1*nsupercell) = zs(1:nsupercell)+k*A*ncxsuper;
            Nsc1 = Nsc1 + 1;
        end
    end
end


% end
rs = [xs;ys;zs]/A;
rlower =[xvoid_lower;yvoid_lower;zvoid_lower]/A;
rupper =[xvoid_upper;yvoid_upper;zvoid_upper]/A;
figure;
plot(rs(1,:),rs(2,:),'o')
hold on
plot(rlower(1),rlower(2),'r*');
plot(rupper(1),rupper(2),'r*');
            