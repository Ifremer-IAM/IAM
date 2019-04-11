Instructions Compilation pour architecture simple :

cd C:\Users\BRI281\Dropbox\These\IAM_Dvt\

%R_pgm_35%/R CMD check IAM65 --no-multiarch

%R_pgm_35%/R CMD build IAM65 --no-multiarch

%R_pgm_35%/R CMD INSTALL --build IAM65 --no-multiarch


Debuggage avec gdb :

R -d gdb -f c:\test.r

Compilation 64 bits

%R_pgm_35%/R CMD INSTALL --build --compile-both IAM65_0.1.tar.gz


