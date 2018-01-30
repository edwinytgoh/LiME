
A general description of the program is given in 'Readme.txt'.


Description of the directories:

./exe, this is where the executable 'MM-INTAS' is put when compiled

./util, contains some tools in order to generate automatically a Makefile

./program

   ./program/README1.txt
     -->  explains how to compile the program

   ./program/src/
     --> subdirectory that contains the f90 source files

   ./program/lib_modules/
     --> subdirectory that contains the f90 source files of some libraries
         (libraries are directly compiled statically with the executable)

   ./program/obj_mod/
     --> subdirectory where the 'obj' and 'mod' files are put after compilation

./test, contains some test cases: Read file './test/README2.txt'!
