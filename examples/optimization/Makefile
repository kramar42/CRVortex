CC=gcc
FLAGS=-Ofast -fomit-frame-pointer -mavx
INCLUDE=-I../../include
LINK=-L../../lib -lfft

all: dynamics

sample: sample.c
	$(CC)  $^ -o $@ $(FLAGS) $(INCLUDE) $(LINK)

dynamics: dynamics.c
	$(CC)  $^ -o $@ -Ofast -fomit-frame-pointer -Wall $(INCLUDE) $(LINK)

clean:
	rm -f sample dynamics
