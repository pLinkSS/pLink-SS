#############################################################
# ProteinIndex 
#
# 2010.4.13
#   Zhou Chen, Wang Leheng
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

SRC    = MetaHandler.cpp NTermOnlineDigestion.cpp CTermOnlineDigestion.cpp OnlineDigestionFactory.cpp MSufSort.cpp PepCalcFunc.cpp ProteinReaderFactory.cpp ProteinHandler.cpp  FastaSQReader.cpp FastaProteinDB.cpp FastaParser.cpp DiskRandomSQReader.cpp DiskRandomACReader.cpp DiskProteinWriter.cpp ProteinWriterFactory.cpp 
LIBS   = -lProteomicsSDK
TARGET = libProteinIndex$(SUFFIX)

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
	
