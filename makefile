l_b = -qopenmp -O3 -xHOST -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_ilp64.a ${MKLROOT}/lib/intel64/libmkl_intel_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -liomp5 -lpthread -lm -ldl

program : main.o MBLAS.o test_output.o
	icpc -o $@ $^ $(l_b)

main.o : main.cpp
	icpc -c $< $(l_b)

MBLAS.o : MBLAS.cpp
	icpc -c $< $(l_b)

test_output.o : test_output.cpp
	icpc -c $< $(l_b)

run : program
	./program

clean:
	rm -f ./program

.PHONY : run clean