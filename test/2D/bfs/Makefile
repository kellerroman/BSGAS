HOMEDIR=../../../

include $(HOMEDIR)Makefile.global

TOOLS_OBJ=$(HOMEDIR)$(TOOLS_DIR)/obj/
CODE_OBJ=$(HOMEDIR)obj/
INCLUDE += -I$(TOOLS_OBJ) -I$(CODE_OBJ)

all: gridgen
	@$(MAKE) -C small
	@$(MAKE) -C large
clean:
	@rm -rf gridgen *.mod
	@$(MAKE) -C small clean
	@$(MAKE) -C large clean

gridgen: gridgen.F90 $(TOOLS_OBJ)/module_gridgen.o
	$(FC) $(CFLAGS) $^ -o $@ $(LIBS) $(INCLUDE)

