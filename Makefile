CFLAGS += -std=c11

CHPLFLAGS += --fast
#CHPLFLAGS += --local
CHPLFLAGS += --savec build

%: %.chpl
	chpl -o $@ $< $(CHPLFLAGS)

lcals: lcals.chpl DataUtils.chpl KernelBase.chpl RunParams.chpl Utils.chpl utils.h
