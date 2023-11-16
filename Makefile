CC = g++
CFLAGS = -Wall -Werror -Wextra -std=c++17
TFLAGS = -lgtest -lstdc++
OBJECTS = s21_matrix_oop.o
SOURCE_DEC = s21_matrix_oop.cpp

# Проверка на mac m1
ifeq ($(shell uname -p), i386)
	CFLAGS = -arch arm64 -Wall -Werror -Wextra -std=c++17 $(shell pkg-config --cflags gtest)
	TFLAGS = $(shell pkg-config --libs gtest)
endif

LIB = s21_matrix_oop.a

all: s21_matrix_oop.a test

s21_matrix_oop.a:
	$(CC) $(CFLAGS) -c $(SOURCE_DEC)
	ar rcs $(LIB) $(OBJECTS)
	ranlib $(LIB)
	rm *.o

test: s21_matrix_oop.cpp test.cpp s21_matrix_oop.h
	@$(CC) $(CFLAGS) s21_matrix_oop.cpp test.cpp $(TFLAGS) -o test
	@./test

clean:
	rm -rf *.o *.gcno *.a *.gcda test test_output

rebuild:
	$(MAKE) clean
	$(MAKE) all

leaks:
	leaks -atExit -- ./test