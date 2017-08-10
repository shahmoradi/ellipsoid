# rm -rf *.txt *.png

main: main.f90
	ifort -debug -llapack -lblas -heap-arrays -g -traceback -check all -fp-stack-check -warn interface main.f90 -o main -lnlopt -lm

# ifort /O2 test_paranest.f90 -o main.exe
# ifort /O3 test_paranest.f90 -o main.exe

#@echo off
# TIMEOUT 1 >nul 2>nul
	ulimit -s 10000 
# TIMEOUT 1 >nul 2>nul
	@echo "done compiling main"
# REM TIMEOUT 1 >nul 2>nul

clean:
	rm -rf main *.mod *.obj *.pdb *.ilk *.smod *.txt

lib:
	rm -rf cobyla2.o
	ifort -c cobyla2.f90 

run:
	./main
