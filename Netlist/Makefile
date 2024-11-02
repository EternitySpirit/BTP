NOCOLOR := \033[0m
LIGHTGREEN := \033[1;32m
LIGHTCYAN := \033[1;36m

SRC := src
SRC_NAME := LP.cpp
SOURCE := $(SRC)/$(SRC_NAME)
BIN_NAME := main
BIN_DIR := out

SYSTEM = unix

.SECONDARY:

OS := $(shell uname -s)
S := /

BUILD_TYPE ?= Release

all: clean build copy_bin

##############
##  SOURCE  ##
##############
.PHONY: build
.PHONY: run

ifeq ($(SOURCE), )

build run:
	$(error no source file provided, the "$@" target must be used like so: \
 make $@ SOURCE=relative/path/to/filename.extension)

else

ifeq ( , $(wildcard $(SOURCE)))
$(error File "$(SOURCE)" does not exist !)
endif
build: $(SOURCE) $(SRC)/CMakeLists.txt clean
	@echo "$(LIGHTCYAN)[Building the executable]$(NOCOLOR)"
	cd $(SRC) && \
 cmake \
 -S. -Bbuild \
 -DCMAKE_BUILD_TYPE=$(BUILD_TYPE) \
 -DCMAKE_PREFIX_PATH=$(CURDIR)
	cd $(SRC) && \
 cmake \
 --build build \
 --config $(BUILD_TYPE) \
 -v
	
run: copy_bin
	@echo "$(LIGHTCYAN)[Running the executable]$(NOCOLOR)"
	@echo "$(LIGHTGREEN)[Output of the executable]$(NOCOLOR)"
	./$(BIN_DIR)/main

endif

#############
##  Tests  ##
#############
.PHONY: test
test:
	$(MAKE) run

.PHONY: clean copy_bin
clean:
	@echo "$(LIGHTGREEN)[Cleaning previous build files]$(NOCOLOR)"
	rm -rf $(SRC)/build
copy_bin: build
	@echo "$(LIGHTGREEN)[Copying the executable]$(NOCOLOR)"
	cp $(SRC)/build/bin/main $(BIN_DIR)/

############
##  MISC  ##
############
# Include user makefile
-include Makefile.user

print-%	: ; @echo $* = $($*)
