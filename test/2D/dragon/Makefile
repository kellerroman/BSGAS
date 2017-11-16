HOMEDIR=../../../

include $(HOMEDIR)Makefile.global

TOOLS_OBJ=$(HOMEDIR)$(TOOLS_DIR)/obj/
CODE_OBJ=$(HOMEDIR)obj/
INCLUDE += -I$(TOOLS_OBJ) -I$(CODE_OBJ)

all: gridgen
	@$(MAKE) -C std
	@$(MAKE) -C 2rows
clean:
	@rm -rf gridgen
	@$(MAKE) -C std clean
	@$(MAKE) -C 2rows clean

gridgen: gridgen.F90 $(TOOLS_OBJ)/module_gridgen.o
	$(FC) $(FFLAGS) $^ -o $@ $(LIBS) $(INCLUDE)

