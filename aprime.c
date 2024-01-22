/*Methods dealing with prime components of Apollonian circle packings*/

/*INCLUSIONS*/

#include <pari.h>
#include "apol.h"
#include "aprime.h"

/*STATIC DECLARATIONS*/
static void thickenedcurvatures_execute(long x[], long res[], long nres, long Bmin, long Bmax)
static int sisprime(long p);

/*MAIN BODY*/

/*Finds the multiplicity of all curvatures in the odd prime component corresponding to v, saving them to a file. B can be an integer or a range.*/
void
thickenedcurvatures(GEN v, GEN B)
{
  pari_sp av = avma;
  long Bmin = 0, Bmax = 0, t = typ(B);/*First, sort out Bmin and Bmax*/
  if (t == t_INT) { Bmin = 1; Bmax = itos(B); }
  else if (t == t_VEC || t == t_COL) {
    if (lg(B) < 3) pari_err_TYPE("B must be a positive integer or a range of positive integers", B);
    if (typ(gel(B, 1)) != t_INT) pari_err_TYPE("B must be a positive integer or a range of positive integers", B);
    if (typ(gel(B, 2)) != t_INT) pari_err_TYPE("B must be a positive integer or a range of positive integers", B);
    Bmin = itou(gel(B, 1));
    Bmax = itou(gel(B, 2));
  }
  else if (t == t_VECSMALL) {
    if (lg(B) < 3) pari_err_TYPE("B must be a positive integer or a range of positive integers", B);
    Bmin = B[1];
    Bmax = B[2];
  }
  else pari_err_TYPE("B must be a positive integer or a range of positive integers", B);
  if (Bmin <= 0) Bmin = 1;
  GEN modres = apol_mod24(v);/*Find it modulo 24*/
  long nres = lg(modres) - 1, i;/*6 or 8*/
  long res[nres];
  for (i = 1; i <= nres; i++) res[i - 1] = itos(gel(modres, i));/*Move to C vector.*/
  long newv[4];
  long iodd = 0, ieven = 2, i;
  for (i = 1; i <= 4; i++) {/*Make the first two odd.*/
    if (Mod2(gel(v, i))) { newv[iodd] = itos(gel(v, i)); iodd++; }
    else { newv[ieven] = itos(gel(v, i)); ieven++; }
  }
  if (!sisprime(newv[0]) && !sisprime(newv[1])) pari_err_TYPE("At least one of the odd numbers must be prime", v);
  thickenedcurvatures_execute(newv, res, nres, Bmin, Bmax);
  set_avma(av);
}

/*Executes thickenedcurvatures. Assumes the first two entries of v are odd, with at least one positive prime.*/
static void
thickenedcurvatures_execute(long x[], long res[], long nres, long Bmin, long Bmax)
{
  long Base = Bmin - (Bmin % 24);/*We want to start at a multiple of 24 to not ruin the mod stuff.*/
  long classmax = (Bmax - Base)/ 24 + 1, i;/*Maximal number of curvatures found in each class.*/
  unsigned int **rclass = (unsigned int **)pari_malloc(24 * sizeof(unsigned int *));/*Stores pointers to the individual classes for curvatures.*/
  if (!rclass) {
    printf("Insufficient memory to allocate to store the residue classes.\n");
    exit(1);
  }
  for (i = 0; i < nres; i++) {
    rclass[res[i]] = (unsigned int *)pari_calloc(classmax * sizeof(unsigned int));/*pari_calloc the classes we want, since we want them as 0 to start.*/
    if (!rclass[res[i]]) {
      printf("Insufficient memory to allocate to store the curvatures.\n");
      exit(1);
    }
  }
  long maxdepth = 200;/*Maximal depth, to start.*/
  long *depthseq = (long *)pari_malloc(maxdepth * sizeof(long));/*depthseq[i] tracks the value we swapped away from in the ith iteration.*/
  if (!depthseq) {
    printf("Insufficient memory to allocate to store the depth sequence.\n");
    exit(1);
  }
  int *swaps = (int *)pari_malloc(maxdepth * sizeof(int));/*Tracks the sequence of swaps, from index 1 to 4.*/
  if (!swaps) {
    printf("Insufficient memory to allocate to store the swaps.\n");
    exit(1);
  }
  int *primes = (int **)pari_calloc(maxdepth * sizeof(int));/*1 if only x[0] is prime, 2 if only x[1] is prime, 3 if both are prime.*/
  if (!primes) {
    printf("Insufficient memory to allocate to store which indices are prime.\n");
    exit(1);
  }
  for (i = 0; i < maxdepth; i++) {
    swaps[i] = -1;/*Initialize to all -1's*/
  }
  if (sisprime(x[0])) primes[0]++;
  if (sisprime(x[1])) primes[1] += 2;
  
  
  for (i = swaps[1] + 1; i < 4; i++) {/*Do the first curvatures.*/
    if (x[i] < Bmin || x[i] > Bmax) continue;
    long shifted = x[i] - Base;
    long b = shifted % 24;
    long a = shifted / 24;/*shifted = 24a + b. b gives the residue block, and a gives the place to insert it.*/
    rclass[b][a]++;
  }
  long ind = 1;/*Which depth we are working at.*/
  if (!sym) {/*Symmetries to worry about. More efficient to do this way, since we don't need to check for symmetries ever otherwise.*/
    while (ind > 0) {/*We are coming in trying to swap this circle out.*/
      int cind = ++swaps[ind];/*Increment the swapping index.*/
      if (cind == 4) {/*Overflowed, go back.*/
        swaps[ind] = -1;
        ind--;
        if (ind < sym) sym = 0;/*We moved past the first index swapping 3, so worry about symmetries again.*/
        continue;
      }
      long lastind = ind - 1;
      if (cind == swaps[lastind]) continue; /*Same thing twice, so skip it.*/
      if (!sym) {/*Worry about symmetries.*/
        if (cind == 0) continue;/*We skip the first one.*/
        else if (cind == 2) sym = ind + 1;/*First time we swap out the third one, eliminating symmetries further on in this branch.*/
      }
      long apbpc = 0;/*Now we can reasonably try a swap.*/
      for (i = 0; i < cind; i++) apbpc += depthseq[lastind][i];
      for (i = cind + 1; i < 4; i++) apbpc += depthseq[lastind][i];
      long newc = (apbpc << 1) - depthseq[lastind][cind];/*2(a+b+c)-d, the new curvature.*/
      if (newc > Bmax) {/*Too big! go back.*/
        if (ind < sym) sym = 0;/*Tried flipping out of symmetry here but it's too big.*/
        continue;
      }
      long shifted = newc - Base;
      if (shifted >= 0) {
        long b = shifted % 24;
        long a = shifted / 24;/*shifted = 24a + b. b gives the residue block, and a gives the place to insert it.*/
        rclass[b][a]++;
      }
      for (i = 0; i < cind; i++) depthseq[ind][i] = depthseq[lastind][i];
      depthseq[ind][cind] = newc;
      for (i = cind + 1; i < 4; i++) depthseq[ind][i] = depthseq[lastind][i];/*Add the tuple in.*/
      ind++;
      if (ind == maxdepth) {/*We are going too deep, must pari_reallocate the storage location.*/
        long newdepth = maxdepth << 1;/*Double it.*/
        depthseq = pari_realloc(depthseq, newdepth * sizeof(long *));
        if (!depthseq) {
          printf("Insufficient memory to pari_reallocate the depth sequence.\n");
          exit(1);
        }
        swaps = pari_realloc(swaps, newdepth * sizeof(int));
        if (!swaps) {
          printf("Insufficient memory to pari_reallocate the swaps.\n");
          exit(1);
        }
        for (i = maxdepth; i < newdepth; i++) {
          depthseq[i] = (long *)pari_malloc(sizeof(long) << 2);
          swaps[i] = -1;
        }
        maxdepth = newdepth;
      }
    } 
  }
  else {/*No symmetry to worry about.*/
    while (ind > 0) {/*We are coming in trying to swap this circle out.*/
      int cind = ++swaps[ind];/*Increment the swapping index.*/
      if (cind == 4) {/*Overflowed, go back.*/
        swaps[ind] = -1;
        ind--;
        continue;
      }
      long lastind = ind - 1;
      if (cind == swaps[lastind]) continue; /*Same thing twice, so skip it.*/
      long apbpc = 0;/*Now we can reasonably try a swap.*/
      for (i = 0; i < cind; i++) apbpc += depthseq[lastind][i];
      for (i = cind + 1; i < 4; i++) apbpc += depthseq[lastind][i];
      long newc = (apbpc << 1) - depthseq[lastind][cind];/*2(a+b+c)-d, the new curvature.*/
      if (newc > Bmax) continue;/*Too big! go back.*/
      long shifted = newc - Base;
      if (shifted >= 0) {
        long b = shifted % 24;
        long a = shifted / 24;/*shifted = 24a + b. b gives the residue block, and a gives the place to insert it.*/
        rclass[b][a]++;
      }
      for (i = 0; i < cind; i++) depthseq[ind][i] = depthseq[lastind][i];
      depthseq[ind][cind] = newc;
      for (i = cind + 1; i < 4; i++) depthseq[ind][i] = depthseq[lastind][i];/*Add the tuple in.*/
      ind++;
      if (ind == maxdepth) {/*We are going too deep, must pari_reallocate the storage location.*/
        long newdepth = maxdepth << 1;/*Double it.*/
        depthseq = pari_realloc(depthseq, newdepth * sizeof(long *));
        if (!depthseq) {
          printf("Insufficient memory to pari_reallocate the depth sequence.\n");
          exit(1);
        }
        swaps = pari_realloc(swaps, newdepth * sizeof(int));
        if (!swaps) {
          printf("Insufficient memory to pari_reallocate the swaps.\n");
          exit(1);
        }
        for (i = maxdepth; i < newdepth; i++) {
          depthseq[i] = (long *)pari_malloc(sizeof(long) << 2);
          swaps[i] = -1;
        }
        maxdepth = newdepth;
      }
    }
  }
  /*Time to free some of the allocated memory.*/
  pari_free(swaps);
  for (i = 0; i < maxdepth; i++) pari_free(depthseq[i]);
  pari_free(depthseq);
  if (tofile) curvs_tofile(rclass, Bmin, Bmax, classmax, v, modres);
  if (tofile % 2) {
    for (i = 0; i < lenr; i++) pari_free(rclass[res[i]]);
    pari_free(rclass);/*The last thing to free*/
    return gc_const(av, gen_0);
  }
  set_avma(av);/*Now we make it into a Vecsmall*/
  long maxncur = classmax * lenr + 1, j;
  GEN curvs = vecsmalltrunc_init(maxncur);
  GEN freqs = vecsmalltrunc_init(maxncur);
  long n = Base - 24;
  for (i = 0; i < classmax; i++) {
    n += 24;/*n = Base + 24i*/
    for (j = 0; j < lenr; j++) {
      long n1 = n + res[j];
      if (n1 < Bmin || n1 > Bmax) continue;/*Too big/small*/
      if (!rclass[res[j]][i]) continue;/*Does not occur*/
      vecsmalltrunc_append(curvs, n1);
      vecsmalltrunc_append(freqs, rclass[res[j]][i]);
    }
  }
  for (i = 0; i < lenr; i++) pari_free(rclass[res[i]]);
  pari_free(rclass);/*The last thing to free*/
  return gerepilecopy(av, mkvec2(curvs, freqs));
}

/*Checks if p is prime but allows for long, not just ulong, and returns 0 if p<=0 (even if -p is prime).*/
static int
sisprime(long p)
{
  if (p <= 0) return 0;
  return uisprime(p);
}

