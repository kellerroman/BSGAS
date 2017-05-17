include Makefile.global

all: BSGAS

BSGAS: FORCE
	@(cd $(OBJECT_DIR) && $(MAKE) -f $(REALMAKEFILE) --no-print-directory)

tools: 
	$(MAKE) -C $(TOOLS_DIR)

clean:
	@rm -rf $(OBJECT_DIR) $(EXECUTABLE_DIR) inc_Makefile_* src/*.mod

##### CREATES NECESSARY FOLDERS
FORCE:
	@mkdir -p $(OBJECT_DIR) $(EXECUTABLE_DIR)


dep: FORCE
	@(cd $(OBJECT_DIR) && $(MAKE) -f $(REALMAKEFILE) ../inc_Makefile_src --no-print-directory)
.PHONY: FORCE
