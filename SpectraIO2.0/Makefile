#############################################################
# SpectraIO
#
# 2013.02.17
#   Jeremy Mung
# 2010.4.13
#   Wang Leheng, Fan Shengbo
# 2009.9.14
#   Fan Shengbo
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

INPUT_DEFINITION  = ./MS2Input/
OUTPUT_DEFINITION = ./MS2Output/

INPUT_SRC     = $(INPUT_DEFINITION)DTAInput.cpp $(INPUT_DEFINITION)DTASInput.cpp $(INPUT_DEFINITION)MGFInput.cpp $(INPUT_DEFINITION)MS2TypeInput.cpp $(INPUT_DEFINITION)PKLInput.cpp $(INPUT_DEFINITION)RAWInput.cpp $(INPUT_DEFINITION)DTASingleInput.cpp $(INPUT_DEFINITION)MS2InputFactory.cpp $(INPUT_DEFINITION)SpecUtility.cpp
OUTPUT_SRC    = $(OUTPUT_DEFINITION)MS2OutputFactory.cpp $(OUTPUT_DEFINITION)MGFOutput.cpp $(OUTPUT_DEFINITION)MS2TypeOutput.cpp $(OUTPUT_DEFINITION)DTAOutput.cpp
LIBS          = -lProteomicsSDK
INPUT_TARGET  = libMS2Input$(SUFFIX)
OUTPUT_TARGET = libMS2Output$(SUFFIX)

INPUT_OBJS    = $(INPUT_SRC:.cpp=.o)
OUTPUT_OBJS   = $(OUTPUT_SRC:.cpp=.o)

all:	$(INPUT_TARGET) $(OUTPUT_TARGET)

clean:
	rm -f $(INPUT_OBJS) $(INPUT_TARGET)
	rm -f $(OUTPUT_OBJS) $(OUTPUT_TARGET)

install:
	cp $(INPUT_TARGET) $(INSTALL_DIR)/$(INPUT_TARGET)
	cp $(OUTPUT_TARGET) $(INSTALL_DIR)/$(OUTPUT_TARGET)
	
$(INPUT_OBJS): $(INPUT_SRC)
	$(CXX) $(CXXFLAGS) $(OBJS_FLAGS) -c $(INPUT_SRC) 
	mv *.o $(INPUT_DEFINITION)
	
$(OUTPUT_OBJS): $(OUTPUT_SRC)
	$(CXX) $(CXXFLAGS) $(OBJS_FLAGS) -c $(OUTPUT_SRC) 
	mv *.o $(OUTPUT_DEFINITION)
	 	
#$(INPUT_TARGET):	$(INPUT_OBJS)
#	$(CXX) $(CXXFLAGS) $(LINK_FLAGS) -L$(BUILD_DIR) -shared -o $(INPUT_TARGET) $(INPUT_OBJS) $(LIBS) 
#	cp $(INPUT_TARGET) $(BUILD_DIR)/$(INPUT_TARGET)
#	
#$(OUTPUT_TARGET):	$(OUTPUT_OBJS)
#	$(CXX) $(CXXFLAGS) $(LINK_FLAGS) -L$(BUILD_DIR) -shared -o $(OUTPUT_TARGET) $(OUTPUT_OBJS) $(LIBS)
#	cp $(OUTPUT_TARGET) $(BUILD_DIR)/$(OUTPUT_TARGET)	 
	
	
ifneq ($(STATICDLL), yes)
$(INPUT_TARGET):	$(INPUT_OBJS)
	$(CXX) $(CXXFLAGS) $(LINK_FLAGS) -L$(BUILD_DIR) -shared -o $(INPUT_TARGET) $(INPUT_OBJS) $(LIBS) 
	cp $(INPUT_TARGET) $(BUILD_DIR)/$(INPUT_TARGET)
	
$(OUTPUT_TARGET):	$(OUTPUT_OBJS)
	$(CXX) $(CXXFLAGS) $(LINK_FLAGS) -L$(BUILD_DIR) -shared -o $(OUTPUT_TARGET) $(OUTPUT_OBJS) $(LIBS)
	cp $(OUTPUT_TARGET) $(BUILD_DIR)/$(OUTPUT_TARGET)	 
else
$(INPUT_TARGET):	$(INPUT_OBJS)
	$(AR) $(ARFLAGS) -o $(INPUT_TARGET) $(INPUT_OBJS) $(LIBS) 
	cp $(INPUT_TARGET) $(BUILD_DIR)/$(INPUT_TARGET)
	
$(OUTPUT_TARGET):	$(OUTPUT_OBJS)
	$(AR) $(ARFLAGS) -o $(OUTPUT_TARGET) $(OUTPUT_OBJS) $(LIBS)
	cp $(OUTPUT_TARGET) $(BUILD_DIR)/$(OUTPUT_TARGET)	 
endif
	
