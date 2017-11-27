HOMEDIR=../../../

include $(HOMEDIR)Makefile.global

TOOLS_OBJ=$(HOMEDIR)$(TOOLS_DIR)/obj/
CODE_OBJ=$(HOMEDIR)obj/
INCLUDE += -I$(TOOLS_OBJ) -I$(CODE_OBJ)

SUBDIRS := ${shell find ./* -maxdepth 0 -type d -print }
SUBCLEAN = $(addsuffix .clean,$(SUBDIRS))
.PHONY: subdirs $(SUBDIRS) clean $(SUBCLEAN)
subdirs: $(SUBDIRS)

$(SUBDIRS): gridgen
	@$(MAKE) -C $@ #--no-print-directory

$(SUBCLEAN): %.clean:
	@$(MAKE) -C $* clean --no-print-directory

clean: $(SUBCLEAN)
	@rm -rf gridgen *.mod

gridgen: gridgen.F90 $(TOOLS_OBJ)/module_gridgen.o
	$(FC) $(FFLAGS) $^ -o $@ $(LIBS) $(INCLUDE)

