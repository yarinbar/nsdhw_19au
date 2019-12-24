CXX = g++

BIN_DIR = bin
TARGET = $(BIN_DIR)/q1
SRC_DIR = src
OBJ_DIR = obj

SRC = $(wildcard $(SRC_DIR)/*.cc)
OBJ = $(SRC:$(SRC_DIR)/%.cc=$(OBJ_DIR)/%.o)

CPPFLAGS += -Iinclude
CFLAGS += -O3 -Wall
LDFLAGS += -lm

.PHONY: all clean

all: $(SRC) $(TARGET)

# link each object into target file
$(TARGET): $(OBJ)
	@mkdir -p $(BIN_DIR)
	$(CXX) $(LDFLAGS) $^ $(LDLIBS) -o $@

# compile each file into object file
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cc
	@mkdir -p $(OBJ_DIR)
	$(CXX) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

run: all
	@./$(TARGET)

clean:
	$(RM) -r $(OBJ_DIR)
	$(RM) -r $(BIN_DIR)
