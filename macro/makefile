# ROOT 配置
ROOTCFLAGS := $(shell root-config --cflags)
ROOTLIBS   := $(shell root-config --libs)

# 编译器
CXX        := g++
CXXFLAGS   := $(ROOTCFLAGS) -pthread -std=c++14 -m64 -D_GLIBCXX_USE_CXX11_ABI=0
LDFLAGS    := $(ROOTLIBS) -lTMVA -lRooFit -lRooFitCore -lMinuit -lMinuit2

# 路径
SRC_DIR    := /data03/tianye/pbar/syscode/macro
BIN_DIR    := /data03/tianye/pbar/syscode/bin

# 源文件与目标
SOURCES    := $(wildcard $(SRC_DIR)/*.C)
EXECUTABLES:= $(patsubst $(SRC_DIR)/%.C, $(BIN_DIR)/%, $(SOURCES))

# 默认目标
all: | $(BIN_DIR) $(EXECUTABLES)

# 确保 bin 目录存在
$(BIN_DIR):
	mkdir -p $@

# 编译、链接规则
$(BIN_DIR)/%: $(SRC_DIR)/%.C
	$(CXX) $< $(CXXFLAGS) -o $@ $(LDFLAGS)

# 清理
clean:
	rm -f $(EXECUTABLES)

.PHONY: all clean


# # 定义变量
# CXX := g++  # 必须使用 g++ 而不是 gcc
# ROOTCFLAGS := $(shell root-config --cflags)
# ROOTLIBS := $(shell root-config --libs)

# # 正确的编译标志
# CXXFLAGS := -pthread -std=c++14 -m64 -D_GLIBCXX_USE_CXX11_ABI=0
# # 正确的链接标志
# LDFLAGS := $(ROOTLIBS) -lTMVA -lRooFit -lRooFitCore -lMinuit2

# # 目录路径
# SRC_DIR := /data03/tianye/pbar/syscode/macro
# BIN_DIR := /data03/tianye/pbar/syscode/bin

# # 查找所有 .C 文件并定义目标
# SOURCES := $(wildcard $(SRC_DIR)/*.C)
# EXECUTABLES := $(patsubst $(SRC_DIR)/%.C, $(BIN_DIR)/%, $(SOURCES))

# # 默认目标
# all: $(EXECUTABLES)

# # 修正后的编译规则
# $(BIN_DIR)/%: $(SRC_DIR)/%.C
# 	@mkdir -p $(@D)
# 	$(CXX) $(CXXFLAGS) $(ROOTCFLAGS) $< $(LDFLAGS) -o $@

# # 清理
# clean:
# 	rm -f $(EXECUTABLES)

# .PHONY: all clean

# # 定义变量
# CXX := g++  # 使用 g++ 替代 gcc
# ROOTCFLAGS := $(shell root-config --cflags)
# ROOTLIBS := $(shell root-config --libs)
# CXXFLAGS := -std=c++14 -pthread -m64 -D_GLIBCXX_USE_CXX11_ABI=0
# LDFLAGS := -lTMVA -lRooFit -lRooFitCore -lMinuit2 -Wl,-rpath,$(shell root-config --libdir)

# # 目录路径
# SRC_DIR := /data03/tianye/pbar/syscode/macro
# BIN_DIR := /data03/tianye/pbar/syscode/bin

# # 查找所有 .C 文件并定义目标
# SOURCES := $(wildcard $(SRC_DIR)/*.C)
# EXECUTABLES := $(patsubst $(SRC_DIR)/%.C, $(BIN_DIR)/%, $(SOURCES))

# # 默认目标
# all: $(EXECUTABLES)

# # 编译规则
# $(BIN_DIR)/%: $(SRC_DIR)/%.C
# 	@mkdir -p $(@D)  # 确保目标目录存在
# 	$(CXX) $(CXXFLAGS) $(ROOTCFLAGS) $< $(LDFLAGS) $(ROOTLIBS) -o $@

# # 清理
# clean:
# 	rm -f $(EXECUTABLES)

# .PHONY: all clean