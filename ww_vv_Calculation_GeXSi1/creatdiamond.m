function [xs ys zs]=creatdiamond(A, ncx, ncy, ncz,ncell,nbasis)
n1=0;

for i=0:ncx-1;
    for j=0:ncy-1;
        for k=0:ncz-1;
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
end
            