ARFLAGS = rcs

all: libcubature.a

libcubature.a: hcubature.o
	$(AR) $(ARFLAGS) $@ $<

hcubature.o: hcubature.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	$(RM) -f *.o
	$(RM) -f *.a
