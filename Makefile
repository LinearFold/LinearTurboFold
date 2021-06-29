################################
# Makefile
#
# author: Sizhen
# edited by: 10/2020
################################

CC=g++
# DEPS=LinearPartitionDouble.h doublebp.cpp energy_parameter.h feature_weight.h intl11.h intl21.h intl22.h utility_v.h utility.h LinearFoldEval.cpp
# DEPS=newbpLog.cpp LinearPartitionLog.h energy_parameter.h feature_weight.h intl11.h intl21.h intl22.h utility_v.h utility.h LinearFoldEval.cpp
# DEPS=LinearPartitionNewprune-C.h energy_parameter.h feature_weight.h intl11.h intl21.h intl22.h utility_v.h utility.h LinearFoldEval.cpp newbp-C.cpp
# DEPS=LinearPartitionNewprune.h energy_parameter.h feature_weight.h intl11.h intl21.h intl22.h utility_v.h utility.h LinearFoldEval.cpp
# DEPS=newbpLog_pf_type.cpp LinearPartitionLog_pf_type.h energy_parameter.h feature_weight.h intl11.h intl21.h intl22.h utility_v.h utility.h LinearFoldEval.cpp
# DEPS= LinearSampling.h energy_parameter.h feature_weight.h intl11.h intl21.h intl22.h utility_v.h utility.h 
DEPS=src/LinearTurboFold.h src/SeqFold.h src/LinearPartition/src/LinearPartition.h src/LinearAlignment/src/LinearAlign.h src/probknot.h \
src/ConfigParser.h src/utils/common_utils.h src/utils/defines.h src/utils/structure_object.h \
src/utils/ansi_string.h src/utils/utils.h src/utils/TProgressDialog.h src/utils/MultiSequence.h src/utils/SafeVector.h src/utils/Sequence.h \
src/utils/structure.h src/utils/rna_library.h \
src/utils/phmm_aln.h src/utils/phmm.h src/utils/p_alignment.h \
src/utils/math/matrix.h src/utils/random.h \
src/ProbabilisticModel.h src/Alignment.h src/utils/GuideTree.h

CFLAGS=-std=c++11 -O3
.PHONY : clean linearturbofold
objects=bin/linearturbofold

linearturbofold: src/LinearTurboFold.cpp $(DEPS) 
		chmod +x linearturbofold
		mkdir -p bin
		$(CC) src/LinearTurboFold.cpp src/SeqFold.cpp src/LinearPartition/src/LinearPartition.cpp src/LinearAlignment/src/LinearAlign.cpp src/probknot.cpp \
		src/ConfigParser.cpp src/utils/common_utils.cpp src/utils/structure_object.cpp \
		src/utils/ansi_string.cpp src/utils/utils.cpp src/utils/TProgressDialog.cpp src/utils/MultiSequence.cpp src/utils/Sequence.cpp \
		src/utils/structure.cpp src/utils/rna_library.cpp \
		src/utils/phmm_aln.cpp src/utils/phmm.cpp src/utils/p_alignment.cpp \
		src/utils/math/matrix.cpp src/utils/random.cpp \
		src/ProbabilisticModel.cpp src/Alignment.cpp src/utils/GuideTree.cpp \
		$(CFLAGS) -Dlv -o bin/linearturbofold

clean:
	-rm $(objects)
