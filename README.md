# Lapack_made_easy
This is a collection of fortran 90 module for simplifying the calling lapack subroutines, especially for the first time users, which is also useful when I am doing performance testing with different linear algebra libraries.


Guidelines:
Read the main.f90 file first for the demo of sloving general real/complex, symmetric real/ hermite complex matrix,
the output of eigen vectors can be optionally specified.

Dependence:

F90 compiler and lapack 

Installation:
  1. Before running the code, make sure you have at least one implemention of lapack library( whether it's intel's mkl or gcc's
openblas, or any other implemention for that matter) properly installed! Details of how to install a lapack implemention can be 
found through google pretty easily, if you have trouble doing that, write an email to me.

  2. Once you have lapack installed, simply run the run.sh if you put the static library in your work directory or make a minor
  modification of the running script accroding to your way of linking lapack, then done!
  
  
 Issues:
   On my laptop, the performace of this code compared with my matlab code is extremly show, actually 1000x slower! Updaye will be made when I
   find out why...
