CFLAGS += -std=c11

CHPLFLAGS += --fast
#CHPLFLAGS += --local
CHPLFLAGS += --savec build
CHPLFLAGS += --print-commands

BINDIR ?= bin

VPATH = $(BINDIR):

%: %.chpl
	chpl -o $(BINDIR)/$@ $< $(CHPLFLAGS)

Executor: Executor.chpl \
	DataUtils.chpl \
	Enums.chpl \
	KernelBase.chpl \
	LongDouble.chpl \
	RunParams.chpl \
	Utils.chpl utils.h \
	algorithm.chpl \
	basic.chpl \
	lcals.chpl \

zzz: zzz.chpl
