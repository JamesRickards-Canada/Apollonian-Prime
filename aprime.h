/*aprime.c methods*/

void thickened(GEN v, GEN B);
GEN thickened_bin(GEN v, unsigned long Bmin, unsigned long binsize, unsigned long nbins);

/*aprime_parallel.c methods*/
GEN parthickenedmult(GEN v, long Bmin, long Bmax, int Nthreads, int load);