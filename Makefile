CHPLFLAGS += --fast
#CHPLFLAGS += --local
CHPLFLAGS += --savec build

%: %.chpl
	chpl -o $@ $< $(CHPLFLAGS)
