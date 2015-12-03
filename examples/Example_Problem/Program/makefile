EXEC = Prob
#
CC = g++
#
CFLAGS = 
#
LIBS = -lm 
#
#
C++FILES = Prob.o update.o

$(EXEC): $(C++FILES)
	@echo "    "
	@echo "Linking..."
	@echo "----------"
	$(CC) -o $(EXEC) $(C++FILES) $(LIBS)

.C.o:
	@echo "   "
	@echo "Compiling Source File --" $<
	@echo "------------------------"
	$(CC) $(C++FLAGS) -c $<

.cpp.o:
	@echo "   "
	@echo "Compiling Numerical Recipe File --" $<
	@echo "----------------------------------"
	$(CC) $(C++FLAGS) -c $<

