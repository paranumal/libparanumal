#####################################################################################
#
#The MIT License (MIT)
#
#Copyright (c) 2017 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus
#
#Permission is hereby granted, free of charge, to any person obtaining a copy
#of this software and associated documentation files (the "Software"), to deal
#in the Software without restriction, including without limitation the rights
#to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#copies of the Software, and to permit persons to whom the Software is
#furnished to do so, subject to the following conditions:
#
#The above copyright notice and this permission notice shall be included in all
#copies or substantial portions of the Software.
#
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#SOFTWARE.
#
#####################################################################################

ifndef LIBP_MAKETOP_LOADED
ifeq (,$(wildcard ../make.top))
$(error cannot locate ${PWD}/../make.top)
else
include ../make.top
endif
endif

#includes
INCLUDES=-I${LIBP_LIBS_DIR}/include \
         ${LIBP_INCLUDES}

#defines
DEFINES =${LIBP_DEFINES} \
         -DLIBP_DIR='"${LIBP_DIR}"' \

#.cpp compilation flags
CORE_CXXFLAGS=${CXXFLAGS} ${DEFINES} ${INCLUDES}

#object dependancies
DEPS=$(wildcard $(LIBP_INCLUDE_DIR)/*.h) \
     $(LIBP_INCLUDE_DIR)/settings.hpp \
     $(LIBP_INCLUDE_DIR)/solver.hpp \
     $(LIBP_INCLUDE_DIR)/timeStepper.hpp

SRC =$(wildcard *.cpp)
OBJS=$(SRC:.cpp=.o)

.PHONY: all lib clean realclean silentUpdate

all: lib

lib: libcore.a silentUpdate

libcore.a: $(OBJS)
ifneq (,${verbose})
	ar -cr libcore.a $(OBJS)
else
	@printf "%b" "$(LIB_COLOR)Building library $(@F)$(NO_COLOR)\n";
	@ar -cr libcore.a $(OBJS)
endif

silentUpdate:
	@true

# rule for .cpp files
%.o: %.cpp $(DEPS)
ifneq (,${verbose})
	$(CXX) -o $*.o -c $*.cpp $(CORE_CXXFLAGS)
else
	@printf "%b" "$(OBJ_COLOR)Compiling $(@F)$(NO_COLOR)\n";
	@$(CXX) -o $*.o -c $*.cpp $(CORE_CXXFLAGS)
endif

#cleanup
clean:
	rm -f *.o libcore.a
