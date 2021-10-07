CHPLFLAGS += --fast

%: %.chpl
	chpl -o $@ $< $(CHPLFLAGS)
