function [xs ys zs,nsupercell]=creatSuperCellwMultipleVoid(A, ncx, ncy, ncz,ncell,nbasis,ncxsuper,ncxvoid,flagsupercell)

if(ncxvoid>1)
	ncxvoid2 = ncxvoid/ncxvoid;
	ncxsuper2 = ncxsuper/ncxvoid;
	ncx2 = ncx/ncxvoid;
	ncy2 = ncy/ncxvoid;
	ncz2 = ncz/ncxvoid;
end

[x2, y2, z2, nsupercell2] = creatSuperCellwVoid(A, ncx2, ncy2, ncz2, ncell, nbasis, ncxsuper2, ncxvoid2, flagsupercell);

Nx = ncx/ncx2;
Ny = ncy/ncy2;
Nz = ncz/ncz2;
if (flagsupercell == 1)
    nsupercell = nsupercell2*(Nx*Ny*Nz);
else
    nsupercell = ncell;
end
Nsc1 = 0;
for i=0:Nx-1; 
    for j=0:Ny-1;
        for k=0:Nz-1;            
            xs((1:nsupercell2)+Nsc1*nsupercell2) = x2(1:nsupercell2)+i*A*ncx2;
            ys((1:nsupercell2)+Nsc1*nsupercell2) = y2(1:nsupercell2)+j*A*ncy2;
            zs((1:nsupercell2)+Nsc1*nsupercell2) = z2(1:nsupercell2)+k*A*ncz2;
            Nsc1 = Nsc1 + 1;
        end
    end
end
Ncount = Nsc1 - nsupercell;

end      
