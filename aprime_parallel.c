/*Parallel implementation of finding the multiplicities of curvatures in a prime component.*/

/*INCLUSIONS*/
#include <pari.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdatomic.h>

typedef struct _primebins_t {/*For doing depth first search on the prime components.*/
  int myid;
  long Bmin;
  long Bmax;
  atomic_uint *found;/*The counts of the found curvatures, using atomic rather than mutexing since it is done so much.*/
  long Ntasks;/*How many sub-tasks we deal out*/
  long *xinit;/*Initial quadruple x*/
  int swapinit;/*The previous swap done before starting this part of the sequence*/
  int primesinit;/*The initial value of primes (with respect to xinit).*/
  long *todo;/*The smallest task not yet done, NULL if done.*/
} primebins_t;

pthread_mutex_t mutex_primebins;


/*Computes the multiplicites of curvatures in the thickened prime component corresponding to v between Bmin and Bmax (inclusive), using a parallel implementation on Nthreads. The data is saved to a file in ./curvcounts, and if load = 1, we also return it as a Vecsmall (not suggested if Bmax-Bmin >= 10^7).*/
GEN
thickened_mult_parallel(GEN v, long Bmin, long Bmax, int Nthreads, int load)
{
  
  /*COPY PASTED SO FAR*/
  
  if (Nthreads <= 1) pari_err_TYPE("Need to use at least 2 threads", stoi(Nthreads));
  long Bmax = Bmin + (nbins * binsize) - 1;/*How far we have to search.*/
  long *counts = (long *)pari_calloc(nbins * sizeof(long));
  long initialdepth = 8, totstart = 1 << initialdepth, i;
  unsigned long **starts = (unsigned long **)pari_malloc(totstart * sizeof(unsigned long *));/*Number of starting places after the breadth first search.*/
  if (!starts) { printf("Insufficient memory to store the starting sequence.\n"); exit(1); }
  for (i = 0; i < totstart; i++) {
	starts[i] = (unsigned long *)pari_malloc(sizeof(unsigned long) << 1);/*Two places for each.*/
  }
  long stop = 1, j;
  starts[0][0] = 0; starts[0][1] = 1;/*Start with 0, 1, then d_i=d_{i-2}+a_i*d_{i-1}.*/
  for (i = 1; i <= initialdepth; i++) {/*Do the breadth first search to the designated depth.*/
	for (j = 0; j < stop; j++) {
	  starts[stop + j][0] = starts[j][1];
	  starts[stop + j][1] = starts[j][0] + (starts[j][1] << 1);/*Hit it with 2.*/
	  if (starts[stop + j][1] <= Bmax) {
	    long shift = starts[stop + j][1] - Bmin;
	    if (shift >= 0) {/*Update the count*/
	      long bl = shift / binsize;
	      counts[bl]++;
	    }
	  }
	  unsigned long tmp = starts[j][1];
	  starts[j][1] += starts[j][0];
	  starts[j][0] = tmp;/*Hit it with 1.*/
	  if (starts[j][1] <= Bmax) {
	    long shift = starts[j][1] - Bmin;
	    if (shift >= 0) {/*Update the count*/
	      long bl = shift / binsize;
	      counts[bl]++;
	    }
	  }
	}
	stop <<= 1;
  }
  pthread_t thread_id[Nthreads];/*The thread ids*/
  long todo = 0;
  zaremba12_t data[Nthreads];
  for (i = 0; i < Nthreads; i++) {/*Initialize the structures holding our data.*/
    data[i].myid = i;
    data[i].Bmin = Bmin;
    data[i].binsize = binsize;
    data[i].nbins = nbins;
    data[i].Bmax = Bmax;
    data[i].counts = counts;
    data[i].Ntasks = stop;
    data[i].starts = starts;
    data[i].todo = &todo;
  }
  pthread_mutex_init(&mutex_zar12, NULL);/*Make the mutex*/
  for (i = 0; i < Nthreads; i++) pthread_create(&thread_id[i], NULL, zaremba_12_par, (void *)&data[i]);/*Make the threads.*/
  for (i = 0; i < Nthreads; i++) pthread_join(thread_id[i], NULL);/*Wait for them to all finish.*/
  pthread_mutex_destroy(&mutex_zar12);/*Eliminate the mutex*/
  for (i = 0; i < totstart; i++) pari_free(starts[i]);
  pari_free(starts);/*Free the memory.*/
  GEN ct = cgetg(nbins + 1, t_VECSMALL);
  for (i = 0; i < nbins; i++) ct[i + 1] = counts[i];
  pari_free(counts);
  return ct;
}