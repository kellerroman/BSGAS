HOMEDIR=../../../../

all: run

run: compile grid.h5
	@$(HOMEDIR)/bin/BSGAS

compile:
	@$(MAKE) -C $(HOMEDIR) --no-print-directory

../gridgen:
	@$(MAKE) -C ../ gridgen

grid.h5: ../gridgen
	@$< 1

clean:
	@rm -vf *.o *.mod *.dat *.plt *.png *.vtk *.txt *.bin *.dx stubex *.csv *.h5 *.xdmf fort.*
	@rm -vf bc.cfg
	@rm -vrf *~

.PHONY: ../gridgen
