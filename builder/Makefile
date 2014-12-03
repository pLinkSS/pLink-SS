#############################################################
# Builder 
#
#  2011.5.16
#     Wang Leheng
#  2008.9.8
#     Chi Hao, Wang Leheng
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

SRC =	builder.cpp

OBJS = $(SRC:.cpp=.o)

LIBS =	$(PTHREAD) -lResultParser

TARGET = builder$(SUFFIX)

all:	$(TARGET)

clean:
	rm -f 	$(OBJS) $(TARGET)

install:
	cp 	$(TARGET) $(INSTALL_DIR)/$(TARGET)
	
$(OBJS): $(SRC)
	$(CXX) $(CXXFLAGS) $(PROFILE) -c $(SRC)
	 	
$(TARGET):	$(OBJS)
	$(CXX) $(CXXFLAGS) $(LINK_FLAGS) -L$(BUILD_DIR) -o $(TARGET) $(OBJS) $(LIBS)
	cp $(TARGET) $(BUILD_DIR)/$(TARGET)
