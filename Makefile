DIR_BITVECTOR := bitvector
DIR_FILES := files
DIR_SDSL := sdsl
DIR_TESTS := tests

override CXX := g++
override CXXFLAGS := -Wall -static -msse4.2 -std=c++11
INCLUDE_FLAGS := \
	-I$(DIR_SDSL)/include
LDFLAGS := \
	-L$(DIR_SDSL)/lib
LDLIBS := \
	-lsdsl \
	-ldivsufsort \
	-ldivsufsort64

ADJUSTMENTS_SDSL := \
	$(DIR_BITVECTOR)/hyb_vector.hpp
CORPUS_NREP := dna/dna.200MB protein/proteins.200MB nlang/english.200MB code/sources.200MB xml/dblp.xml.200MB
CORPUS_REP := para cere influenza world_leaders kernel

CORPUS_LIST := $(DIR_FILES)/list
CORPUS_REP_LIST := $(DIR_FILES)/list_repcorpus

BIN_NAMES := \
	test_select \
	bench_plcp \
	bench_bwt_select
BIN_RELEASE   := $(foreach item,$(BIN_NAMES),$(DIR_TESTS)/$(item))

# All
.PHONY: all
all: $(BIN_RELEASE)

# Libraries
$(DIR_SDSL): $(ADJUSTMENTS_SDSL)
	if [ ! -d $@ ]; then \
		rm -rf $@ $@_src ; \
		git clone --depth 1 --branch v2.1.1 https://github.com/simongog/sdsl-lite.git $@_src ; \
		cp $^ ./$@_src/include/sdsl ; \
		> ./$@_src/build/clean.sh ; \
		./$@_src/install.sh ./$@ ; \
		rm -rf $@_src ; \
	else \
		cp $^ $@/include/sdsl ; \
	fi
	touch $@

# Files
$(CORPUS_LIST) $(CORPUS_REP_LIST):
	mkdir -p $(DIR_FILES)
	>$@

.PHONY: download
download:
	mkdir -p $(DIR_FILES)
	>$(CORPUS_LIST)
	>$(CORPUS_REP_LIST)

	for t in $(CORPUS_NREP) ; do \
		f=$$(echo $$t | cut -d "/" -f 2); \
		if ! [ -f $(DIR_FILES)/$$f ] ; then \
			wget -P $(DIR_FILES) https://pizzachili.dcc.uchile.cl/texts/$$t.gz ; \
			gunzip $(DIR_FILES)/$$f.gz; \
		fi ; \
		echo $$f >> $(CORPUS_LIST) ; \
	done

	for f in $(CORPUS_REP) ; do \
		if ! [ -f $(DIR_FILES)/$$f ] ; then \
			wget -P $(DIR_FILES) https://pizzachili.dcc.uchile.cl/repcorpus/real/$$f.gz ; \
			gunzip $(DIR_FILES)/$$f.gz; \
		fi ; \
		echo $$f >> $(CORPUS_LIST) ; \
		echo $$f >> $(CORPUS_REP_LIST) ; \
	done

# Binaries
$(BIN_RELEASE): override CXXFLAGS += -funroll-loops -O3
$(BIN_RELEASE): override CPPFLAGS += -DNDEBUG

tests/test_select: tests/test_select.cpp sdsl
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(INCLUDE_FLAGS) $< -o $@ $(LDFLAGS) $(LDLIBS)

tests/bench_plcp: tests/bench_lcp.cpp sdsl $(CORPUS_LIST)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(INCLUDE_FLAGS) -DUSE_PLCP $< -o $@ $(LDFLAGS) $(LDLIBS)

tests/bench_bwt_select: tests/bench_bwt_select.cpp sdsl $(CORPUS_LIST)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(INCLUDE_FLAGS) $< -o $@ $(LDFLAGS) $(LDLIBS)

.PHONY: clean
clean:
	rm -f $(BIN_RELEASE) $(DIR_TESTS)/*.o *.sdsl

.PHONY: clean_files
clean_files:
	rm -rf $(DIR_FILES)

.PHONY: clean_all
clean_all: clean_output clean_files
