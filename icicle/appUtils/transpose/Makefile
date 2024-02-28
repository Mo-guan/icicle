test_transpose: test.cu transpose.cu kernels.cu
	nvcc -o test_transpose -I. -I../.. test.cu
	./test_transpose

clean:
	rm -rf test_transpose