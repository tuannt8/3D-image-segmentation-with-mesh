# 
TARGET = DSC

BIN_DIR=build
LIBSRCS = $(wildcard ../src/*.cpp) $(wildcard ../CGLA/*.cpp) $(wildcard ../is_mesh/*.cpp) $(wildcard ../DEMO/*.cpp) $(wildcard ../SCGenerator/*.cpp) $(wildcard ../TetGen/*.cpp)
LIBC =  $(wildcard ../SOIL/*.c)
LIBOBJS = $(LIBSRCS:.cpp=.o) $(LIBC:.c=.o) 
DEPENDS = ${OBJECTS:.o=.d}

WARNING= -w


CXX = g++
CXXFLAGS = -Wall -O3 -std=c++11 -fpermissive $(WARNING) -I ../src -I../CGLA -I../SOIL -I../is_mesh -I../ -I/usr/include -I../CImg -I../SCGenerator -I../TetGen -MMD
CFLAGS= $(CXXFLAGS)
OPENGL= -lm -lGL -L/usr/X11R6/lib -lGLU -lglut -lGLEW -lXi -lXmu
BOOST=  -L/usr/lib -lboost_filesystem -lboost_system
IMAGE= -lpng
LDFLAGS = $(OPENGL) $(BOOST) -lpthread $(IMAGE)

CC = g++

all:$(TARGET)
$(TARGET): $(LIBOBJS)
	$(CXX) $(CXXFLAGS) $(CFLAGS) $(LDFLAGS) -o $@ $(LIBOBJS) 

-include ${DEPENDS}

test:
	echo $(LIBOBJS)

clean_common:
	rm ../DEMO/segment_function.o ../DEMO/user_interface.o ../src/DSC.o

.PHONY: clean
clean: 
	rm $(LIBOBJS)

