backup_dir := backup
build_dir := build
bitvector_dir := bitvector
external_dir := external
rsrc_dir := files
test_dir := tests

CXX := g++
CXXFLAGS := -Wall -static -std=c++20 -funroll-loops -O3 -march=native
CPPFLAGS := -DNDEBUG

# FILES
bitvector_files := $(wildcard $(bitvector_dir)/**)
library_files := $(wildcard $(build_dir)/include/**)

#######
# All #
#######

.PHONY: all
all: \
	$(build_dir)/test_serialize \
	$(build_dir)/test_select \
	$(build_dir)/bench_plcp \
	$(build_dir)/bench_bwt_select

################
# BENCH_HEADER #
################

$(build_dir)/bench_header: override CPPFLAGS += -I$(bitvector_dir)
$(build_dir)/bench_header: $(build_dir)/%: $(test_dir)/%.cpp $(bitvector_files)
	mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $< -o $@

################
# BENCH_FOOTER #
################

$(build_dir)/bench_footer: override CPPFLAGS += -I$(bitvector_dir)
$(build_dir)/bench_footer: $(build_dir)/%: $(test_dir)/%.cpp $(bitvector_files)
	mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $< -o $@

#########
# SETUP #
#########

.PHONY: setup
setup: sdsl zombit_sdsl pasta la_vector download

########
# SDSL #
########

sdsl_dir := $(external_dir)/sdsl-lite

.PHONY: sdsl
sdsl:
	mkdir -p $(build_dir)
	-git submodule update --init $(sdsl_dir)
	./$(sdsl_dir)/install.sh $(build_dir)

######################
# ZOMBIT VECTOR SDSL #
######################

zombit_dir := $(external_dir)/zombit-vector
zombit_sdsl_dir := $(zombit_dir)/external/sdsl-lite

.PHONY: zombit_sdsl
zombit_sdsl:
	mkdir -p $(build_dir)/zombit-vector
	-git submodule update --init --recursive $(zombit_dir)
	# Use the original install script, so cmake is called directly. (fixes error on this device)
	cp $(external_dir)/sdsl-lite/install.sh $(zombit_sdsl_dir)/install.sh
	./$(zombit_sdsl_dir)/install.sh $(build_dir)/zombit-vector

#########
# PASTA #
#########

pasta_dir := $(external_dir)/pasta
pasta_utils_dir := $(external_dir)/pasta_utils
tlx_dir := $(external_dir)/tlx

.PHONY: pasta
pasta:
	-git submodule update --init --recursive $(pasta_dir)
	-git submodule update --init --recursive $(pasta_utils_dir)
	-git submodule update --init --recursive $(tlx_dir)

#############
# LA_VECTOR #
#############

la_vector_dir := $(external_dir)/la_vector

.PHONY: la_vector
la_vector:
	-git submodule update --init --recursive $(la_vector_dir)

############
# DOWNLOAD #
############

corpus_nrep := dna/dna.200MB protein/proteins.200MB nlang/english.200MB code/sources.200MB xml/dblp.xml.200MB
corpus_rep := para cere influenza world_leaders kernel

corpus_list := $(rsrc_dir)/list

.PHONY: download
download:
	mkdir -p $(rsrc_dir)
	>$(corpus_list)

	for t in $(corpus_nrep) ; do \
		f=$$(echo $$t | cut -d "/" -f 2); \
		if ! [ -f $(rsrc_dir)/$$f ] ; then \
			wget -P $(rsrc_dir) https://pizzachili.dcc.uchile.cl/texts/$$t.gz ; \
			gunzip $(rsrc_dir)/$$f.gz; \
		fi ; \
		echo $$f >> $(corpus_list) ; \
	done

	for f in $(corpus_rep) ; do \
		if ! [ -f $(rsrc_dir)/$$f ] ; then \
			wget -P $(rsrc_dir) https://pizzachili.dcc.uchile.cl/repcorpus/real/$$f.gz ; \
			gunzip $(rsrc_dir)/$$f.gz; \
		fi ; \
		echo $$f >> $(corpus_list) ; \
	done

###############
# TEST_SELECT #
###############

$(build_dir)/test_select: override CPPFLAGS += -I$(build_dir)/include
$(build_dir)/test_select: override LDFLAGS += -L$(build_dir)/lib
$(build_dir)/test_select: override LDLIBS += -lsdsl -ldivsufsort -ldivsufsort64
$(build_dir)/test_select: $(build_dir)/%: $(test_dir)/%.cpp $(bitvector_files) $(library_files)
	mkdir -p $(dir $@)
	cp $(bitvector_files) $(build_dir)/include/sdsl/
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $< -o $@ $(LDFLAGS) $(LDLIBS)

##################
# TEST_SERIALIZE #
##################

$(build_dir)/test_serialize: override CPPFLAGS += -I$(build_dir)/include
$(build_dir)/test_serialize: override LDFLAGS += -L$(build_dir)/lib
$(build_dir)/test_serialize: override LDLIBS += -lsdsl -ldivsufsort -ldivsufsort64
$(build_dir)/test_serialize: $(build_dir)/%: $(test_dir)/%.cpp $(bitvector_files) $(library_files)
	mkdir -p $(dir $@)
	cp $(bitvector_files) $(build_dir)/include/sdsl/
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $< -o $@ $(LDFLAGS) $(LDLIBS)

##############
# BENCH_PLCP #
##############

$(build_dir)/bench_plcp: override CPPFLAGS += -I$(build_dir)/include -I$(external_dir) 
# $(build_dir)/bench_plcp: override CPPFLAGS += -I$(build_dir)/zombit-vector/include -I$(external_dir)/zombit-vector/include -I$(build_dir)/include -I$(external_dir) 
# $(build_dir)/bench_plcp: override CPPFLAGS += -I$(external_dir)/pasta/include -I$(external_dir)/pasta_utils/include -I$(external_dir)/tlx -I$(build_dir)/include -I$(external_dir) 
$(build_dir)/bench_plcp: override LDFLAGS += -L$(build_dir)/lib
$(build_dir)/bench_plcp: override LDLIBS += -lsdsl -ldivsufsort -ldivsufsort64
$(build_dir)/bench_plcp: $(build_dir)/%: $(test_dir)/%.cpp $(bitvector_files) $(library_files)
	mkdir -p $(dir $@)
	cp $(bitvector_files) $(build_dir)/include/sdsl/
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $< -o $@ $(LDFLAGS) $(LDLIBS)

####################
# BENCH_BWT_SELECT #
####################

$(build_dir)/bench_bwt_select: override CPPFLAGS += -I$(build_dir)/include
# $(build_dir)/bench_bwt_select: override CPPFLAGS += -I$(build_dir)/zombit-vector/include -I$(external_dir)/zombit-vector/include -I$(build_dir)/include -I$(external_dir) 
# $(build_dir)/bench_bwt_select: override CPPFLAGS += -I$(external_dir)/pasta/include -I$(build_dir)/include -I$(external_dir) 
$(build_dir)/bench_bwt_select: override LDFLAGS += -L$(build_dir)/lib
$(build_dir)/bench_bwt_select: override LDLIBS += -lsdsl -ldivsufsort -ldivsufsort64
$(build_dir)/bench_bwt_select: $(build_dir)/%: $(test_dir)/%.cpp $(bitvector_files) $(library_files)
	mkdir -p $(dir $@)
	cp $(bitvector_files) $(build_dir)/include/sdsl/
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $< -o $@ $(LDFLAGS) $(LDLIBS)

#########
# CLEAN #
#########

.PHONY: clean
clean:
	if [ -d $(build_dir) ]; then find $(build_dir) -maxdepth 1 -type f -delete; fi

###########
# NUCLEAR #
###########

.PHONY: nuclear
nuclear:
	rm -rf $(build_dir)

