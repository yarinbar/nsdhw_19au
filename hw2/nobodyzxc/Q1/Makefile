TARGET   = Q1_main.exe
SOURCES  = Line Q1_main

O_DIR    = O_$(HOSTTYPE)
LIBS     = -lm

CXX      = g++
CXXFLAGS = -g -Wno-deprecated -Wall #-D__STDC__
TARG     = $(TARGET)

SRCS     = ${SOURCES:%=%.cc}
OBJS     = ${SRCS:%.cc=$(O_DIR)/%.o}

FLAGS    = -O3 -Wall -shared -std=c++11 -fPIC

$(TARG) : $(O_DIR) $(OBJS) $(TESTOBJ)
	$(CXX) -o $@ $(CXXFLAGS) $(OBJS) $(TESTOBJ) $(LIBS)

$(O_DIR) :
	mkdir $@

$(O_DIR)/%.o: %.cc
	$(CXX) $(CXXFLAGS) -c -o $@ $<

clean:
	rm -f $(OBJS) $(TESTOBJ) $(TARG)
	rm -rf $(O_DIR) $(BINDOBJ) $(BINDALIAS)

run: $(TARG)
	@./$(TARGET)
