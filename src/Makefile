# rm -rf *.txt *.png

# Modiy NLOPTDIR to your desired path for installation
NLOPTDIR = /Users/travisdriver/Documents/Research/NLOPT
# Modiy NLOPTSRC to the path where NLOPT source files exist
NLOPTSRC = ../lib/nlopt-2.4.2/

main: main.f90
	#ifort -debug -llapack -lblas -heap-arrays -g -traceback -check all -fp-stack-check -warn interface main.f90 -o main -lnlopt -lm
	ifort -O2 main.f90 -o main -lnlopt -lm
	#ifort -O2 -mkl main.f90 -o main -L$(NLOPTDIR)/lib -lnlopt -lm
	@echo "DONE with O2 compilation"

debug:
	ifort -debug -mkl -heap-arrays -g -traceback -check all -fp-stack-check -warn interface main.f90 -o main -L$(NLOPTDIR)/lib -lnlopt -lm
	@echo "DONE with debug compilation"
# REM TIMEOUT 1 >nul 2>nul

clean:
	rm -rf main *.mod *.obj *.pdb *.ilk *.smod *.txt
	cd ..; cd out; rm -rf *.txt *.png

nlopt:
	rm -rf $(NLOPTDIR)/*; cd $(NLOPTSRC); chmod 777 configure; ./configure --prefix=$(NLOPTDIR); make; sudo make install

run:
	#ulimit -s unlimited; 
	./main
