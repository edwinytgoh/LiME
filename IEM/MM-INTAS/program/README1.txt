
How to compile the program?

By executing the following commands in the current directory:
1.   ../util/Make_Filelist
2.   ../util/Make_Makefile
3.   make

Explanation:
1. makes a list of the source code files --> creates file 'allfiles.txt'
2. makes the Makefile --> creates 'Makefile'
   The program 'f90deps' should be available in order to run this
   It can be generated in the directory '../util' by 'Make_f90deps'
3. creates the executable --> creates '../exe/MM-INTAS'
