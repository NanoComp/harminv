HARMINV_OBJECTS = harminv-main.o harminv.o

CFLAGS = -g

harminv: $(HARMINV_OBJECTS)
	$(CC) $(CFLAGS) -o harminv $(HARMINV_OBJECTS) -llapack -lf77blas -lcblas -latlas -lg2c -lm

clean:
	rm -rf core $(HARMINV_OBJECTS) harminv
