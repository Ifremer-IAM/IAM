Instructions Compilation pour architecture simple :

cd C:\Users\BRI281\Dropbox\These\IAM_Dvt\

%R_pgm_35%/R CMD check IAMMSpTAC --no-multiarch

%R_pgm_35%/R CMD build IAMMSpTAC --no-multiarch

%R_pgm_35%/R CMD INSTALL --build IAMMSpTAC --no-multiarch


Debuggage avec gdb :

R -d gdb -f c:\test.r

Compilation 64 bits

%R_pgm_35%/R CMD INSTALL --build --compile-both IAMMSpTAC_0.1.tar.gz


install.packages("/home1/datahome/fbriton/AMURE/SESSF/IAMMSpTAC_0.1.tar.gz", repos=NULL, type="source")