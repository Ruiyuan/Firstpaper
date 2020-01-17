The folder is not super clean. 
1. But I used two methods "FBZ" and "Isotropic" to calculate thermal conductivity "k" for bulk Si and Si/Ge QDSL. 
2. Occasionally, you can also see some files are diserpersion curve calculation, which is the first step in thermal conductivity calculation.
3. files with "ww_vv" are to calculate frequency and velocity from dispersion curves, for thermal conductivity in the future. Because it's time consuming and geometry dependent, sometime I splite the "ww_vv" calculation over several files.

Same as last changes. These additional files help calculate dispersion curve and thermal conductivity.
A lot of redundent files, but look at the following ones:
1. "Test_conductivity_CH_Primitive_110direction.m	"
2. "Master_conductivity_CH_SuperCell_SiGe_Gillet_JHT_100.m	"
3. "fmat" are the force constant matrix data, which was calculated and saved from similar codes in the folder "fmatSuperCell"
4. Folder "ww_vv_BulkSi" are used to calculate frequency and velocity for bulk Si. The code in the “ww_vv” calculation need to read data from force constants data from “fmat” type data or codes
5. Briefly, the idea in the codes design are:
	(1). Create atom supercell geometry
              (2). Calculate dispersion curve (dispersion curve code) from force constants (fmat, force constants code)
              (3) Calculate frequency and velocity from dispersion curves using central difference method
              (4) Put the frequency from dispersion curve and velocity at each q point into the FBZ or Isotropic model using some phone relaxation time parameters from literature papers to calculate thermal conducitivity. 

Some of the MATLAB codes needs to read codes in othre folders with an address. You may need to update the address given that the orginal folder is a mess.
