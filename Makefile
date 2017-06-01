include Makefile.global

all: BSGAS tools

BSGAS: FORCE
	@(cd $(OBJECT_DIR) && $(MAKE) -f $(REALMAKEFILE) --no-print-directory)

tools: 
	$(MAKE) -C $(TOOLS_DIR)

clean:
	@rm -rf $(OBJECT_DIR) $(EXECUTABLE_DIR) inc_Makefile_* src/*.mod
	$(MAKE) -C $(TOOLS_DIR) clean
	$(MAKE) -C test clean

##### CREATES NECESSARY FOLDERS
FORCE:
	@mkdir -p $(OBJECT_DIR) $(EXECUTABLE_DIR)


dep: FORCE
	@(cd $(OBJECT_DIR) && $(MAKE) -f $(REALMAKEFILE) ../inc_Makefile_src --no-print-directory)

test:
	@$(MAKE) -C test
.PHONY: FORCE tools dep test
