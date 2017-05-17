
include ../Makefile.global

all: tools

compile_source:
	@$(MAKE) -C ../ --no-print-directory BSGAS

tools: FORCE #compile_source
	@(cd $(OBJECT_DIR) && $(MAKE) -f $(REALMAKEFILE) --no-print-directory)

clean: FORCE
	@rm -rf $(OBJECT_DIR) $(EXECUTABLE_DIR) *~
	
FORCE:
	@mkdir -p $(OBJECT_DIR) $(EXECUTABLE_DIR)

.PHONY: FORCE
