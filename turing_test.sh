gcc allvec_new.c -lm -O3 -o allvec.out
gcc nsvd1* -lm -O3 -o nsvd1.out
gcc nsvd2* -lm -O3 -o nsvd2.out
gcc nsvd3* -lm -O3 -o nsvd3.out
./nsvd1.out tmp/bigfile/main.dat 123 0 10 > tmp1
./nsvd3.out tmp/bigfile/main.dat 123 0 10 > tmp3
    printf("usage: ./main.out\n path of main.dat\n random seed\n percentage of qualify\n number of training iteration\n");
