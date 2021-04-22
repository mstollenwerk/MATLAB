I did a lot of trial and error and so I don't know if what worked in the end only worked because of sth I did before. Anyways, here are the steps:

1. install msys2 and mingw64 using https://stackoverflow.com/questions/30069830/how-to-install-mingw-w64-and-msys2 . GMP and MPFR are already installed alongside with it.
2. install make command in msys2 terminal. 
3. download flint, arb libraries.
4. use mingw64 terminal in msys folder to do
   4.1. run ./configure in flint and arb libraries with options --with-gmp=/c/msys64/mingw64 --with-mpfr=/c/msys2/mingw64 as well as somthing like --dynamic off --- static on. 
   4.2 make 
   4.3 make check [takes forever]
   4.4 make install
5. wrote kummerUc.c function [a lot of work!]
6. in matlab run setenv('MW_MINGW64_LOC', 'C:\msys64\mingw64') to set compiler to the same one I used to build packages in 4.
7. mex mex -R2018a kummerUc.c -IC:\msys64\mingw64\include -IC:\msys64\usr\local\include -LC:\msys64\mingw64\lib -LC:\msys64\usr\local\lib -larb -lflint -lmpfr -lgmp
