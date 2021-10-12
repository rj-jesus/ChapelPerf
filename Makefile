CHPLFLAGS += --fast
#CHPLFLAGS += -L${HOME}/.local/lib

%: %.chpl
	chpl -o $@ $< $(CHPLFLAGS)
