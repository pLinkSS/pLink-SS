#############################################################
# XLinkResultFilter 
#
# 2010.4.13
#   Wang Leheng
# 2008.9.8
#   Chi Hao, Wang Leheng
#
#############################################################

include ../make.include

ifdef IS_LINUX
	SUFFIX     =
	LINK_FLAGS =	
endif

ifdef IS_WINDOWS
	SUFFIX     = .exe
	LINK_FLAGS = -Wl,--enable-auto-import -Wl,--enable-runtime-pseudo-reloc		
endif

SRC    = main.cpp FilterConf.cpp XLinkResultFilter.cpp XLinkResultFDRFilter.cpp XLinkResultFilterFactory.cpp
LIBS   = -lProteomicsSDK -lXLinkPepResultFilter -lXLinkResultReport
TARGET = XLinkResultFilter$(SUFFIX)

OBJS   = $(SRC:.cpp=.o)

all:	$(TARGET)

clean:
	rm -f $(OBJS) $(TARGET)

install:
	cp $(TARGET) $(INSTALL_DIR)/$(TARGET)
	
$(OBJS): $(SRC)
	$(CXX) $(CXXFLAGS) -c $(SRC)

$(TARGET):	$(OBJS)
	$(CXX) $(CXXFLAGS) $(LINK_FLAGS) -L$(BUILD_DIR) -o $(TARGET) $(OBJS) $(LIBS)
	cp $(TARGET) $(BUILD_DIR)/$(TARGET)
	

	
