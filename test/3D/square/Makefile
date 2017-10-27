HOMEDIR=../../../

include $(HOMEDIR)Makefile.global

TOOLS_OBJ=$(HOMEDIR)$(TOOLS_DIR)/obj/
CODE_OBJ=$(HOMEDIR)obj/
INCLUDE += -I$(TOOLS_OBJ) -I$(CODE_OBJ)

all: gridgen
	@$(MAKE) -C 1block
	@$(MAKE) -C turned
	@$(MAKE) -C 2blocks
	@$(MAKE) -C 4blocks
	@$(MAKE) -C 8blocks
clean:
	@rm -rf gridgen
	@$(MAKE) -C 1block clean
	@$(MAKE) -C turned clean
	@$(MAKE) -C 2blocks clean
	@$(MAKE) -C 4blocks clean
	@$(MAKE) -C 8blocks clean

gridgen: gridgen.F90 $(TOOLS_OBJ)/module_gridgen.o
	$(FC) $(FFLAGS) $^ -o $@ $(LIBS) $(INCLUDE)

