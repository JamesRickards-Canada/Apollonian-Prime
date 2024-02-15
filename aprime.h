/*aprime.c methods*/
int sisprime(long p);
GEN primeroots_bin(GEN v, unsigned long Bmin, unsigned long binsize, unsigned long nbins);
void thickened(GEN v, GEN B);
GEN thickened_bin(GEN v, unsigned long Bmin, unsigned long binsize, unsigned long nbins);

/*aprime_parallel.c methods*/
GEN parthickenedmult(GEN vgen, long Bmin, long Bmax, int Nthreads, int load);