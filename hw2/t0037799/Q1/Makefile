CXX = g++
SRC_DIR = src
BUILD_DIR ?= build
CXXFILES = $(wildcard $(SRC_DIR)/*.cc)
OBJFILES = $(CXXFILES:$(SRC_DIR)/%.cc=$(BUILD_DIR)/%.o)
CXXFLAGS = -O2 -Wall -Wextra -Werror

TARGET = $(BUILD_DIR)/a.out

.PHONY: all
all: $(BUILD_DIR) $(TARGET)

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

$(TARGET): $(OBJFILES)
	$(CXX) -o $@ $(OBJFILES)

$(BUILD_DIR)/%.o : $(SRC_DIR)/%.cc
	$(CXX) -o $@ -Iinclude $(CXXFLAGS) -c $^ 

.PHONY: clean
clean:
	rm -r $(BUILD_DIR)

.PHONY: run
run:
	@ ./$(TARGET)
