#############################################################
# SearchEngine
#
# 2010.4.13
#   Wang Leheng
# 2010.3.15
#   Wu Yanjie, Chi Hao, Wang Leheng
#
#############################################################

include ../make.include

ifdef IS_LINUX
	SUFFIX     = .so
	PTHREAD    = -lpthread	
	OBJS_FLAGS = -fPIC
	LINK_FLAGS =
endif

ifdef IS_WINDOWS
	SUFFIX     = .dll
	PTHREAD    = -lpthreadGCE2	
	OBJS_FLAGS =
	LINK_FLAGS = -Wl,--enable-auto-import -Wl,--enable-runtime-pseudo-reloc	
endif

SRC    = SearchEngine.cpp
LIBS   = $(PTHREAD) -lProteomicsSDK -lFlow -lMS2Input -lProteinIndex -lProteinInfer
TARGET = libSearchEngine$(SUFFIX)

OBJS   = $(SRC:.cpp=.o)

all:	$(TARGET)

clean:
	rm -f $(OBJS) $(TARGET)	

install:
	cp $(TARGET) $(INSTALL_DIR)/$(TARGET)
	
$(OBJS): $(SRC)
	$(CXX) $(CXXFLAGS) $(OBJS_FLAGS) -c $(SRC)
	 	
#$(TARGET):	$(OBJS)
#	$(CXX) $(CXXFLAGS) $(LINK_FLAGS) -L$(BUILD_DIR) -shared -o $(TARGET) $(OBJS) $(LIBS)
#	cp $(TARGET) $(BUILD_DIR)/$(TARGET)

	
ifneq ($(STATICDLL), yes)
$(TARGET):	$(OBJS)
	$(CXX) $(CXXFLAGS) $(LINK_FLAGS) -L$(BUILD_DIR) -shared -o $(TARGET) $(OBJS) $(LIBS)
	cp $(TARGET) $(BUILD_DIR)/$(TARGET)
else
$(TARGET):	$(OBJS)
	$(AR) $(ARFLAGS) -o $(TARGET) $(OBJS)
	cp $(TARGET) $(BUILD_DIR)/$(TARGET)
endif
	
