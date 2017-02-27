CC = g++
BUILD_PATH = ./build
CFLAGS = -O3 -Wall --std=c++14 -I.

tests : $(BUILD_PATH)/test
	$(BUILD_PATH)/test

$(BUILD_PATH)/test : build
	$(CC) -o $@ $(CFLAGS) ./tests/test.cpp

build :
	mkdir $(BUILD_PATH);

.PHONY: tests build
