#############################################################
# FLOW 
# 2014.11.25
#  Fan Shengbo
# 2010.4.13
#   Wang Leheng, Fan Shengbo
# 2008.9.8
#   Chi Hao, Wang Leheng
#
#############################################################

include ../make.include

ifdef IS_LINUX
	SUFFIX     = .so
	OBJS_FLAGS = -fPIC
	LINK_FLAGS =	
endif

ifdef IS_WINDOWS
	SUFFIX     = .dll
	OBJS_FLAGS =
	LINK_FLAGS = -Wl,--enable-auto-import -Wl,--enable-runtime-pseudo-reloc	
endif

SRC    = XLinkOpenFlow.cpp XLinkPepFlow.cpp PepPairList.cpp SpectraWindow.cpp FlowFactory.cpp SpectraSearcher.cpp common.cpp  
LIBS   = -lProteomicsSDK -lProteinIndex -lMass2PepIndex -lPreProcess -lPeptideGenerator -lPepPairContainer -lXLinkPFDIO -lXLinkPeptideScore -lXLinkPeptideEvaluate -lXLinkResultReport -lProteinInfer
TARGET = libFlow$(SUFFIX)

OBJS   = $(SRC:.cpp=.o)

all:	$(TARGET)

clean:
	rm -f $(OBJS) $(TARGET)

install:
	cp $(TARGET) $(INSTALL_DIR)/$(TARGET)
	
$(OBJS): $(SRC)
	$(CXX) $(CXXFLAGS) $(OBJS_FLAGS) -c $(SRC)
	 	
ifneq ($(STATICDLL), yes)
$(TARGET):	$(OBJS)
	$(CXX) $(CXXFLAGS) $(LINK_FLAGS) -L$(BUILD_DIR) -shared -o $(TARGET) $(OBJS) $(LIBS)
	cp $(TARGET) $(BUILD_DIR)/$(TARGET)
else
$(TARGET):	$(OBJS)
	$(AR) $(ARFLAGS) -o $(TARGET) $(OBJS)
	cp $(TARGET) $(BUILD_DIR)/$(TARGET)
endif
