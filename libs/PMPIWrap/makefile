WRAP := ../../3rdParty/wrap/wrap.py
CC   := mpicc

PMPIWrap.a: PMPIWrap.o
	$(AR) -cr libPMPIWrap.a PMPIWrap.o

PMPIWrap.o: PMPIWrap.c
	$(CC) $(CFLAGS) -o PMPIWrap.o -c PMPIWrap.c

clean:
	rm -f PMPIWrap.a PMPIWrap.o
