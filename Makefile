#############################################################
# pLink Makefile
#
# 2014.11.25
#   Fan Shengbo
#
#############################################################

include ./make.include

MAKE:= make

all:
	cd ProteomicsSDK && $(MAKE)
	cd ProteinIndex && $(MAKE)
	cd Mass2PepIndex && $(MAKE)
	cd Indexer && $(MAKE)
	cd SpectraIO2.0 && $(MAKE)
	cd PreProcess && $(MAKE)
	cd PeptideGenerator && $(MAKE)
	cd PepPairContainer && $(MAKE)
	cd XLinkPFDIO && $(MAKE)
	cd XLinkScore && $(MAKE)
	cd XLinkEvaluate && $(MAKE)
	cd XLinkResultReport && $(MAKE)
	cd ProteinInfer && $(MAKE)
	cd Flow && $(MAKE)
	cd SearchEngine && $(MAKE)

ifdef IS_WINDOWS
	cd Searcher && $(MAKE)
endif	
	cd XLinkProteinInfer && $(MAKE)
	cd Inferor && $(MAKE)
	cd XLinkPepResultFilter && $(MAKE)
	cd XLinkResultFilter && $(MAKE)
	cd XLinkResultConverter && $(MAKE)
	cd ResultParse && $(MAKE)
	cd builder && $(MAKE)
	cd PSMAnalysis && $(MAKE)

install:
	cd ProteomicsSDK && $(MAKE) install
	cd ProteinIndex && $(MAKE) install
	cd Mass2PepIndex && $(MAKE) install
	cd Indexer && $(MAKE) install
	cd SpectraIO2.0 && $(MAKE) install
	cd PreProcess && $(MAKE) install
	cd PeptideGenerator && $(MAKE) install
	cd PepPairContainer && $(MAKE) install
	cd XLinkPFDIO && $(MAKE) install
	cd XLinkScore && $(MAKE) install
	cd XLinkEvaluate && $(MAKE) install
	cd XLinkResultReport && $(MAKE) install
	cd ProteinInfer && $(MAKE) install
	cd Flow && $(MAKE) install
	cd SearchEngine && $(MAKE) install
ifdef IS_WINDOWS
	cd Searcher && $(MAKE) install
endif	
	cd XLinkProteinInfer && $(MAKE) install
	cd Inferor && $(MAKE) install
	cd XLinkPepResultFilter && $(MAKE) install
	cd XLinkResultFilter && $(MAKE) install
	cd XLinkResultConverter && $(MAKE) install
	cd ResultParse && $(MAKE) install
	cd builder && $(MAKE) install
	cd PSMAnalysis && $(MAKE) install

clean:
	cd ProteomicsSDK && $(MAKE) clean
	cd ProteinIndex && $(MAKE) clean
	cd Mass2PepIndex && $(MAKE) clean
	cd Indexer && $(MAKE) clean
	cd SpectraIO2.0 && $(MAKE) clean
	cd PreProcess && $(MAKE) clean
	cd PeptideGenerator && $(MAKE) clean
	cd PepPairContainer && $(MAKE) clean
	cd XLinkPFDIO && $(MAKE) clean
	cd XLinkScore && $(MAKE) clean
	cd XLinkEvaluate && $(MAKE) clean
	cd XLinkResultReport && $(MAKE) clean
	cd ProteinInfer && $(MAKE) clean
	cd Flow && $(MAKE) clean
	cd SearchEngine && $(MAKE) clean
ifdef IS_WINDOWS
	cd Searcher && $(MAKE) clean
endif	
	cd XLinkProteinInfer && $(MAKE) clean
	cd Inferor && $(MAKE) clean
	cd XLinkPepResultFilter && $(MAKE) clean
	cd XLinkResultFilter && $(MAKE) clean
	cd XLinkResultConverter && $(MAKE) clean
	cd ResultParse && $(MAKE) clean
	cd builder && $(MAKE) clean
	cd PSMAnalysis && $(MAKE) clean
