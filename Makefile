CFLAGS += -std=c11

CHPLFLAGS += --fast
#CHPLFLAGS += --local
CHPLFLAGS += --savec build
CHPLFLAGS += --print-commands

%: %.chpl
	chpl -o $@ $< $(CHPLFLAGS)

chpl-perf: chpl-perf.chpl \
	DataTypes.chpl \
	DataUtils.chpl \
	Enums.chpl \
	Executor.chpl \
	KernelBase.chpl \
	LongDouble.chpl \
	RunParams.chpl \
	Utils.chpl \
	algorithm.chpl \
	apps.chpl \
	basic.chpl \
	lcals.chpl \
	polybench.chpl \
	stream.chpl \

zzz: zzz.chpl
