gcc allvec_new.c -lm -O3 -o allvec.out
gcc nsvd1* -lm -O3 -o nsvd1.out
gcc nsvd2* -lm -O3 -o nsvd2.out
gcc nsvd3* -lm -O3 -o nsvd3.out
./nsvd1.out tmp/bigfile/main.dat 123 0 1 > tmp1
./nsvd3.out tmp/bigfile/data1.dat tmp/bigfile/data2.dat tmp/bigfile/datav.dat 1 1 1 1 > tmp3
