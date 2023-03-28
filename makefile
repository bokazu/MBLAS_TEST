l_b =  -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_ilp64.a ${MKLROOT}/lib/intel64/libmkl_intel_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group 
program : main.o MBLAS.o test_output.o
	icpx -DMKL_ILP64  -I"${MKLROOT}/include" -qopenmp -O3 -xHOST -o  $@ $^ $(l_b)

main.o : main.cpp
	icpx -DMKL_ILP64  -I"${MKLROOT}/include" -qopenmp -O3 -xHOST -c  $< $(l_b)

MBLAS.o : MBLAS.cpp
	icpx -DMKL_ILP64  -I"${MKLROOT}/include" -qopenmp -O3 -xHOST  -c $< $(l_b)

test_output.o : test_output.cpp
	icpx -DMKL_ILP64  -I"${MKLROOT}/include" -qopenmp -O3 -xHOST  -c $< $(l_b)

run : program
	./program

clean:
	rm -f ./program

.PHONY : run clean