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


install_superlu:
	git clone https://github.com/xiaoyeli/superlu.git
	(cd superlu && mkdir build && cd build && cmake .. -DCMAKE_INSTALL_PREFIX=../build && make && make install)

clean_superlu:
	rm -rf superlu

.PHONY: FORCE tools dep test
