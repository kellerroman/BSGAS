HOMEDIR=../../../

include $(HOMEDIR)Makefile.global

TOOLS_OBJ=$(HOMEDIR)$(TOOLS_DIR)/obj/
CODE_OBJ=$(HOMEDIR)obj/
INCLUDE += -I$(TOOLS_OBJ) -I$(CODE_OBJ)

all: run

gridgen: gridgen.F90 $(TOOLS_OBJ)/module_gridgen.o
	$(FC) $(FFLAGS) $^ -o $@ $(LIBS) $(INCLUDE)


run: compile grid.h5
	@$(HOMEDIR)/bin/BSGAS

compile:
	@$(MAKE) -C $(HOMEDIR) --no-print-directory

grid.h5: gridgen
	@./$<

compare: compile grid.h5
	@cp config.cfg config.old
	#@echo "grid_out = grid_out.h5" >> config.cfg
	@sed -i "s/wall_move_rest = 2/wall_move_rest = 1/g" config.cfg
	@$(HOMEDIR)/bin/BSGAS
	@mv grid_out.h5 grid_1.h5
	@$(HOMEDIR)/$(TOOLS_DIR)/bin/create_xdmf grid_1.h5 grid_1.xdmf
	@sed -i "s/wall_move_rest = 1/wall_move_rest = 2/g" config.cfg
	@$(HOMEDIR)/bin/BSGAS
	@mv grid_out.h5 grid_2.h5
	@$(HOMEDIR)/$(TOOLS_DIR)/bin/create_xdmf grid_2.h5 grid_2.xdmf
	@mv config.old config.cfg
	

clean:
	@rm -vf *.o *.mod *.dat *.plt *.png *.vtk *.txt *.bin gridgen *.dx stubex *.csv *.h5 *.xdmf fort.*
	@rm -vf bc.cfg
	@rm -vrf *~

.PHONY: compile 
