# Compilador
CXX = g++

# Flags de compilação
CXXFLAGS = -Wall -Wextra -std=c++17 -Iinclude

# Diretório de build
BUILD_DIR = build

# Nome do executável
TARGET = $(BUILD_DIR)/my_program

# Lista de arquivos fonte
SRCS = src/broyden.cpp src/ns.cpp src/setup.cpp src/main.cpp

# Lista de arquivos objeto (gerados no diretório build)
OBJS = $(patsubst src/%.cpp, $(BUILD_DIR)/%.o, $(SRCS))

# Regra padrão
all: $(TARGET)

# Regra para compilar o executável
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJS)

# Regra para compilar arquivos objeto
$(BUILD_DIR)/%.o: src/%.cpp
	@mkdir -p $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Regra para limpar os arquivos gerados
clean:
	rm -rf $(BUILD_DIR)
	echo "Build files removed."

# Regra para executar o programa
run: $(TARGET)
	./$(TARGET)

# Garante que as regras all, clean e run não sejam confundidas com arquivos
.PHONY: all clean run