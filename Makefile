CFLAGS += -std=c11

CHPLFLAGS += --fast
#CHPLFLAGS += --local
CHPLFLAGS += --savec build
CHPLFLAGS += --print-commands

BINDIR ?= bin

VPATH = $(BINDIR):

%: %.chpl
	chpl -o $(BINDIR)/$@ $< $(CHPLFLAGS)

lcals: lcals.chpl DataUtils.chpl KernelBase.chpl RunParams.chpl Utils.chpl utils.h LongDouble.chpl

zzz: zzz.chpl DataUtils.chpl KernelBase.chpl RunParams.chpl Utils.chpl utils.h
