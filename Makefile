OPT	= -Wall -O0 -g

ACC_OBJ = acc2tax.o

all:remove_objects $(ACC_OBJ)
	mkdir -p $(BIN); $(CC) $(OPT) -lm -o acc2tax $(ACC_OBJ)

clean:
	rm *.o
	rm -rf acc2tax

remove_objects:
	rm *.o

%.o : %.c
	mkdir -p obj; $(CC) -Iinclude $(OPT) -c $< -o $@

