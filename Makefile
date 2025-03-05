# Compiler and flags
CXX = g++
CXXFLAGS = -std=c++11 -Wall -O2

# Executable name
TARGET = broyden

# Source and object files
SRC = main.cpp broyden.cpp
OBJ = $(SRC:.cpp=.o)

# Default rule: build the executable
all: $(TARGET)

# Compile object files
%.o: %.cpp broyden.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Link object files into the final executable
$(TARGET): $(OBJ)
	$(CXX) $(CXXFLAGS) $(OBJ) -o $(TARGET)

# Run the program
run: $(TARGET)
	./$(TARGET)

# Clean build files
clean:
	rm -f $(OBJ) $(TARGET)

# Clean and rebuild
rebuild: clean all
