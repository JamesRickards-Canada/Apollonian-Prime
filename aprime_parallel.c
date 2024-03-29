/*Parallel implementation of finding the multiplicities of curvatures in a prime component.*/

/*INCLUSIONS*/
#include <pari.h>
#include "aprime.h"
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdatomic.h>
#include <time.h>

typedef struct _thickmults_t {/*For doing depth first search on the prime components.*/
  int myid;
  long Bmin;
  long Bmax;
  atomic_uint *found;/*The counts of the found curvatures, using atomic rather than mutexing since it is done so much.*/
  long Ntasks;/*How many sub-tasks we deal out*/
  long **starts;;/*Stores the initial quadruples*/
  int *primesinit;/*Stores the initial value of primes*/
  int *swapsinit;/*Stores the initial value of swaps*/
  long *todo;/*The smallest task not yet done, -1 if done.*/
} thickmults_t;

pthread_mutex_t mutex_thickmults;

static void *thickened_mult_par(void *args);

/*Computes the multiplicites of curvatures in the thickened prime component corresponding to v between Bmin and Bmax (inclusive), using a parallel implementation on Nthreads. The data is saved to a file in ./curvcounts, and if load = 1, we also return it as a Vecsmall (not suggested if Bmax-Bmin >= 10^7).*/
GEN
parthickenedmult(GEN vgen, long Bmin, long Bmax, int Nthreads, int load)
{
  pari_sp av = avma;
  struct timespec start, finish;/*For measuring the real time.*/
  clock_gettime(CLOCK_MONOTONIC, &start);
  if (Nthreads <= 1) pari_err_TYPE("Need to use at least 2 threads", stoi(Nthreads));
  long x[4];
  long iodd = 0, ieven = 2, i;
  for (i = 1; i <= 4; i++) {/*Make the first two odd.*/
    if (Mod2(gel(vgen, i))) x[iodd++] = itos(gel(vgen, i));
    else x[ieven++] = itos(gel(vgen, i));
  }
  if (!sisprime(x[0]) && !sisprime(x[1])) pari_err_TYPE("At least one of the odd numbers must be prime", vgen);
  /*Now we have checked the initial input, time to initialize the multiplicity storing variable.*/
  long nB = Bmax - Bmin + 1;/*Number of multiplicities we are keeping track of*/
  atomic_uint *found = (atomic_uint *)pari_calloc(nB * sizeof(atomic_uint));/*Stores the found curvatures, we won't be going past a uint size.*/
  /*We start with a search to depth 9, making at least 512 different threads to search through. At each stage, there are 2 or 3 numbers we could possibly flip, and we go to depth 9.*/
  long maxdepthp1 = 10, maxtasks = 20000;/*maxtasks > 3^9, and we expect to end up somewhere in the middle. We store maxdepth+1*/
  long **starts = (long **)pari_malloc(maxtasks * sizeof(long *));/*Stores the initial quadruples*/
  int *primesinit = (int *)pari_malloc(maxtasks * sizeof(int));/*Stores the initial value of primes*/
  int *swapsinit = (int *)pari_malloc(maxtasks * sizeof(int));/*Stores the initial value of swaps*/
  /*Now to initialize the initial searching*/
  long *depthseq = (long *)pari_malloc(maxdepthp1 * sizeof(long));/*depthseq[i] tracks the value we swapped away from in the ith iteration.*/
  int *swaps = (int *)pari_malloc(maxdepthp1 * sizeof(int));/*Tracks the sequence of swaps, from index 1 to 4.*/
  int *primes = (int *)pari_calloc(maxdepthp1 * sizeof(int));/*1 if only x[0] is prime, 2 if only x[1] is prime, 3 if both are prime.*/
  for (i = 0; i < maxdepthp1; i++) swaps[i] = -1;/*Initialize to all -1's*/
  if (sisprime(x[0])) primes[0]++;
  if (sisprime(x[1])) primes[0] += 2;
  for (i = 0; i < 4; i++) {/*Do the first curvatures.*/
    if (x[i] < Bmin || x[i] > Bmax) continue;
    long shifted = x[i] - Bmin;
    found[shifted]++;
  }
  /*Let's do the search!*/
  long ind = 1, sind = 0;/*Which depth we are working at, and which index in starts, primesinit, and swapsinit we are at.*/
  long v[4] = {x[0], x[1], x[2], x[3]};/*Initial quadruple.*/
  for (;;) {/*We are coming in trying to swap a circle out. Guaranteed that ind > 0 entering the loop.*/
    long lastind = ind - 1;
    if (ind == maxdepthp1) {/*We reached the end! Add this to our list.*/
      starts[sind] = (long *)pari_malloc(sizeof(long) << 2);/*Store the quadruple.*/
      for (i = 0; i <= 3; i++) starts[sind][i] = v[i];
      primesinit[sind] = primes[lastind];
      swapsinit[sind] = swaps[lastind];
      sind++;/*No need to check if this exceeds maxtasks, we chose maxtasks large enough to be OK.*/
      ind--;/*We assume maxdepth > 0, so we don't enter the for loop with ind=0. It would be silly to have maxdepth=0.*/
      v[swaps[ind]] = depthseq[ind];/*Update our v backwards to the correct thing.*/
      continue;
    }
    /*Continue with the normal depth first search.*/
    int cind = ++swaps[ind];/*Increment the swapping index.*/
    if (cind == 4) {/*Overflowed, go back.*/
      swaps[ind] = -1;
      ind--;
      if (!ind) break; /*All done!!*/
      v[swaps[ind]] = depthseq[ind];/*Update our v backwards to the correct thing.*/
      continue;
    }
    if (cind == swaps[lastind]) continue; /*Same thing twice, so skip it.*/
    if (!cind && primes[lastind] == 1) continue;/*Swapping out the only prime, not allowed!*/
    if (cind == 1 && primes[lastind] == 2) continue;/*Swapping out the only prime, not allowed!*/
    long apbpc = 0;/*Now we can reasonably try a swap.*/
    for (i = 0; i < cind; i++) apbpc += v[i];
    for (i = cind + 1; i < 4; i++) apbpc += v[i];
    long newc = (apbpc << 1) - v[cind];/*2(a+b+c)-d, the new curvature.*/    
    if (newc > Bmax) continue;/*Too big! go back.*/
    long shifted = newc - Bmin;
    if (shifted >= 0) found[shifted]++;/*Update found.*/
    depthseq[ind] = v[cind];/*Store the value we swapped away from.*/
    v[cind] = newc;/*Update v*/
    switch (cind) {/*Update primes*/
      case 0:/*Check if prime*/
        if (sisprime(newc)) primes[ind] = 3;/*1st index must be prime as we swapped the 0th*/
        else primes[ind] = 2;
        break;
      case 1:/*Check if prime*/
        if (sisprime(newc)) primes[ind] = 3;/*0th index must be prime as we swapped the 1st*/
        else primes[ind] = 1;
        break;
      default:
        primes[ind] = primes[lastind];/*Odd numbers did not change.*/
    }
    ind++;
  }
  /*Starting data initialized!! Let's initialize the threads.*/
  pthread_t thread_id[Nthreads];/*The thread ids*/
  long todo = 0;
  thickmults_t data[Nthreads];
  for (i = 0; i < Nthreads; i++) {/*Initialize the structures holding our data.*/
    data[i].myid = i;
    data[i].Bmin = Bmin;
    data[i].Bmax = Bmax;
    data[i].found = found;
    data[i].Ntasks = sind;
    data[i].starts = starts;
    data[i].primesinit = primesinit;
    data[i].swapsinit = swapsinit;
    data[i].todo = &todo;
  }
  pthread_mutex_init(&mutex_thickmults, NULL);/*Make the mutex*/
  for (i = 0; i < Nthreads; i++) pthread_create(&thread_id[i], NULL, thickened_mult_par, (void *)&data[i]);/*Make the threads.*/
  for (i = 0; i < Nthreads; i++) pthread_join(thread_id[i], NULL);/*Wait for them to all finish.*/
  pthread_mutex_destroy(&mutex_thickmults);/*Eliminate the mutex*/
  /*Let's write to the file first.*/
  if (!pari_is_dir("thickcurvcounts")) {
    int s = system("mkdir -p thickcurvcounts");
    if (s == -1) pari_err(e_MISC, "ERROR CREATING DIRECTORY thickcurvcounts");
  }
  char *filename = stack_sprintf("thickcurvcounts/%Pd_%Pd_%Pd_%Pd_%ld-to-%ld.dat", gel(vgen, 1), gel(vgen, 2), gel(vgen, 3), gel(vgen, 4), Bmin, Bmax);
  FILE *F = fopen(filename, "w");/*Created the output file f*/
  for (i = 0; i < nB; i++) pari_fprintf(F, "%d\n", found[i]);
  fclose(F);
  /*Print the time taken.*/
  clock_gettime(CLOCK_MONOTONIC, &finish);
  double elapsed = (finish.tv_sec - start.tv_sec);
  elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
  pari_printf("Total computation time: %fs\n", elapsed);
  /*Free all the memory except for found.*/
  set_avma(av);/*All info we need is not on the stack.*/
  pari_free(primes);
  pari_free(swaps);
  pari_free(depthseq);
  pari_free(swapsinit);
  pari_free(primesinit);
  for (i = 0; i < sind; i++) pari_free(starts[i]);
  pari_free(starts);
  GEN ret = gen_1;
  if (load) {/*Make the return Vecsmall if requested*/
    ret = cgetg(nB + 1, t_VECSMALL);
    for (i = 1; i <= nB; i++) ret[i] = found[i - 1];
  }
  pari_free(found);/*Free found, the last unfreed variable.*/
  return ret;
}

/*Executes the depth first search in parallel.*/
static void *
thickened_mult_par(void *args)
{
  /*Import the given data*/
  thickmults_t *data = (thickmults_t *)args;
  atomic_uint *found = data -> found;
  long Bmin = data -> Bmin, Bmax = data -> Bmax, Ntasks = data -> Ntasks;
  long **starts = data -> starts;
  int *primesinit = data -> primesinit, *swapsinit = data -> swapsinit;
  /*Initialize the required data for the depth first search.*/
  long maxdepth = 200, i;/*Maximal depth, to start.*/
  long *depthseq = (long *)pari_malloc(maxdepth * sizeof(long));/*depthseq[i] tracks the value we swapped away from in the ith iteration.*/
  int *swaps = (int *)pari_malloc(maxdepth * sizeof(int));/*Tracks the sequence of swaps, from index 1 to 4.*/
  int *primes = (int *)pari_malloc(maxdepth * sizeof(int));/*1 if only x[0] is prime, 2 if only x[1] is prime, 3 if both are prime.*/
  /*Time to execute the searches.*/
  for (;;) {/*First, let's figure out which case we are doing.*/
    pthread_mutex_lock(&mutex_thickmults);/*Retrieve the next index*/
    long *next = data -> todo;
    if (*next == -1) { pthread_mutex_unlock(&mutex_thickmults); break; }/*All tasks done or running. Don't forget to unlock!!*/
    long startsind = *next;
    (*next)++;
    if (*next == Ntasks) *next = -1;/*Last one.*/
    pthread_mutex_unlock(&mutex_thickmults);
    /*Initialize the starting data*/
    for (i = 1; i < maxdepth; i++) {
      swaps[i] = -1;/*Initialize to all -1's*/
      primes[i] = 0;/*All 0 to start.*/
    }
    swaps[0] = swapsinit[startsind];/*Only extra one to preset.*/
    primes[0] = primesinit[startsind];/*Same here*/
    long *x = starts[startsind];
    long ind = 1;/*Which depth we are working at.*/
    long v[4] = {x[0], x[1], x[2], x[3]};/*Initial quadruple.*/
    for (;;) {/*We are coming in trying to swap a circle out. Guaranteed that ind > 0 entering the loop.*/
      int cind = ++swaps[ind];/*Increment the swapping index.*/
      if (cind == 4) {/*Overflowed, go back.*/
        swaps[ind] = -1;
        ind--;
        if (!ind) break;/*All done!!*/
        v[swaps[ind]] = depthseq[ind];/*Update our v backwards to the correct thing.*/
        continue;
      }
      long lastind = ind - 1;
      if (cind == swaps[lastind]) continue; /*Same thing twice, so skip it.*/
      if (!cind && primes[lastind] == 1) continue;/*Swapping out the only prime, not allowed!*/
      if (cind == 1 && primes[lastind] == 2) continue;/*Swapping out the only prime, not allowed!*/
      long apbpc = 0;/*Now we can reasonably try a swap.*/
      for (i = 0; i < cind; i++) apbpc += v[i];
      for (i = cind + 1; i < 4; i++) apbpc += v[i];
      long newc = (apbpc << 1) - v[cind];/*2(a+b+c)-d, the new curvature.*/    
      if (newc > Bmax) continue;/*Too big! go back.*/
      long shifted = newc - Bmin;
      if (shifted >= 0) found[shifted]++;/*Update found.*/
      depthseq[ind] = v[cind];/*Store the value we swapped away from.*/
      v[cind] = newc;/*Update v*/
      switch (cind) {/*Update primes*/
        case 0:/*Check if prime*/
          if (sisprime(newc)) primes[ind] = 3;/*1st index must be prime as we swapped the 0th*/
          else primes[ind] = 2;
          break;
        case 1:/*Check if prime*/
          if (sisprime(newc)) primes[ind] = 3;/*0th index must be prime as we swapped the 1st*/
          else primes[ind] = 1;
          break;
        default:
          primes[ind] = primes[lastind];/*Odd numbers did not change.*/
      }
      ind++;
      if (ind == maxdepth) {/*We are going too deep, must pari_reallocate the storage location.*/
        long newdepth = maxdepth << 1;/*Double it.*/
        depthseq = pari_realloc(depthseq, newdepth * sizeof(long));
        swaps = pari_realloc(swaps, newdepth * sizeof(int));
        for (i = maxdepth; i < newdepth; i++) swaps[i] = -1;
        primes = pari_realloc(primes, newdepth * sizeof(int));
        maxdepth = newdepth;
      }
    }
  }
  /*Time to free the allocated memory.*/
  pari_free(primes);
  pari_free(swaps);
  pari_free(depthseq);
  return NULL;
}

