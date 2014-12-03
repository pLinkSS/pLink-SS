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

PATH1 = ./BasicFunction/
PATH2 = ./EnginesParser/
PATH3 = ./ExportResult/
PATH4 = ./DataProcessing/

NET_SRC = ReadWrite.cpp	CommonProcess.cpp	DataAnalysis.cpp StringUtility.cpp	pFindParser.cpp	MascotParser.cpp	SEQUESTParser.cpp	SQTParser.cpp	ExportHTML.cpp	ExportpBuildPre.cpp	ExportForInterface.cpp	ProteinInference.cpp	ExportpepXML.cpp	MergeTopFind.cpp		ExportToIntermediateFile.cpp	FilterUtility.cpp	DataIntegrate.cpp SampleInfor.cpp

SRC =	BasicFunction/ReadWrite.cpp	BasicFunction/CommonProcess.cpp	BasicFunction/DataAnalysis.cpp BasicFunction/StringUtility.cpp	EnginesParser/pFindParser.cpp	EnginesParser/MascotParser.cpp	EnginesParser/SEQUESTParser.cpp	EnginesParser/SQTParser.cpp	ExportResult/ExportHTML.cpp	ExportResult/ExportToIntermediateFile.cpp	ExportResult/ExportpBuildPre.cpp	ExportResult/ProteinInference.cpp	ExportResult/ExportForInterface.cpp	ExportResult/ExportpepXML.cpp	ExportResult/MergeTopFind.cpp	DataProcessing/FilterUtility.cpp	DataProcessing/DataIntegrate.cpp SampleInfor.cpp

LIBS = 
OBJS =	$(NET_SRC:.cpp=.o)
TARGET = libResultParser$(SUFFIX)

all: $(TARGET)
	
clean:
	rm -f $(OBJS) $(TARGET)

install:
	cp $(TARGET) $(INSTALL_DIR)/$(TARGET)
	
$(OBJS): $(SRC)
	$(CXX) $(CXXFLAGS) $(OBJS_FLAGS) -c $(SRC)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) $(LINK_FLAGS) -L$(BUILD_DIR) -shared -o $(TARGET) $(OBJS) $(LIBS)
	cp $(TARGET) $(BUILD_DIR)/$(TARGET)
