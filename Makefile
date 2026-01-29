# OS に応じて voro++ の include/lib を切り替え
# macOS: Homebrew (ARM /opt/homebrew, Intel /usr/local)
# Linux (Ubuntu 等): apt で入れた libvoro++-dev の標準パス
UNAME := $(shell uname -s)
ifeq ($(UNAME),Darwin)
  VORO_ROOT := $(shell ( [ -d /opt/homebrew/opt/voro++ ] && echo /opt/homebrew/opt/voro++ ) || echo /usr/local/opt/voro++)
  VORO_INC := $(VORO_ROOT)/include/voro++
  VORO_LIB := $(VORO_ROOT)/lib
else
  # Ubuntu/Debian: libvoro++1 + libvoro++-dev で /usr/include/voro++, ライブラリはシステム検索
  VORO_INC := /usr/include/voro++
  VORO_LIB :=
endif

easy_voro/voro_volumes: src/voro_volumes.cpp
	g++ -O2 -std=c++17 src/voro_volumes.cpp \
		-I$(VORO_INC) \
		$(if $(VORO_LIB),-L$(VORO_LIB)) \
		-lvoro++ \
		-o easy_voro/voro_volumes

# クラウドボリューム（Dropbox 等）上では .venv をプロジェクト内に作らない
install:
	POETRY_VIRTUALENVS_IN_PROJECT=false poetry install

clean:
	- find . -name Icon\* -delete