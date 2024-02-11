/*Methods dealing with prime components of Apollonian circle packings*/

/*INCLUSIONS*/

#include <pari.h>
#include "apol.h"
#include "aprime.h"

/*STATIC DECLARATIONS*/
static void thickened_execute(long x[], long res[], long nres, long Bmin, long Bmax, GEN vorig);
static GEN thickened_bin_execute(long x[], unsigned long Bmin, unsigned long binsize, unsigned long nbins, GEN vorig);

/*MAIN BODY*/

/*Finds the multiplicity of all curvatures in the odd prime component corresponding to v, saving them to a file. B can be an integer or a range.*/
void
thickened(GEN v, GEN B)
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
  long iodd = 0, ieven = 2;
  for (i = 1; i <= 4; i++) {/*Make the first two odd.*/
    if (Mod2(gel(v, i))) newv[iodd++] = itos(gel(v, i));
    else newv[ieven++] = itos(gel(v, i));
  }
  if (!sisprime(newv[0]) && !sisprime(newv[1])) pari_err_TYPE("At least one of the odd numbers must be prime", v);
  thickened_execute(newv, res, nres, Bmin, Bmax, v);
  set_avma(av);
}

/*Executes thickened. Assumes the first two entries of v are odd, with at least one positive prime.*/
static void
thickened_execute(long x[], long res[], long nres, long Bmin, long Bmax, GEN vorig)
{
  long Base = Bmin - (Bmin % 24);/*We want to start at a multiple of 24 to not ruin the mod stuff.*/
  long classmax = (Bmax - Base)/ 24 + 1, i;/*Maximal number of curvatures found in each class.*/
  unsigned int **rclass = (unsigned int **)pari_malloc(24 * sizeof(unsigned int *));/*Stores pointers to the individual classes for curvatures.*/
  if (!rclass) {
    pari_printf("Insufficient memory to allocate to store the residue classes.\n");
    exit(1);
  }
  for (i = 0; i < nres; i++) {
    rclass[res[i]] = (unsigned int *)pari_calloc(classmax * sizeof(unsigned int));/*pari_calloc the classes we want, since we want them as 0 to start.*/
    if (!rclass[res[i]]) {
      pari_printf("Insufficient memory to allocate to store the curvatures.\n");
      exit(1);
    }
  }
  long maxdepth = 200;/*Maximal depth, to start.*/
  long *depthseq = (long *)pari_malloc(maxdepth * sizeof(long));/*depthseq[i] tracks the value we swapped away from in the ith iteration.*/
  if (!depthseq) {
    pari_printf("Insufficient memory to allocate to store the depth sequence.\n");
    exit(1);
  }
  int *swaps = (int *)pari_malloc(maxdepth * sizeof(int));/*Tracks the sequence of swaps, from index 1 to 4.*/
  if (!swaps) {
    pari_printf("Insufficient memory to allocate to store the swaps.\n");
    exit(1);
  }
  int *primes = (int *)pari_calloc(maxdepth * sizeof(int));/*1 if only x[0] is prime, 2 if only x[1] is prime, 3 if both are prime.*/
  if (!primes) {
    pari_printf("Insufficient memory to allocate to store which indices are prime.\n");
    exit(1);
  }
  for (i = 0; i < maxdepth; i++) swaps[i] = -1;/*Initialize to all -1's*/
  if (sisprime(x[0])) primes[0]++;
  if (sisprime(x[1])) primes[0] += 2;
  for (i = 0; i < 4; i++) {/*Do the first curvatures.*/
    if (x[i] < Bmin || x[i] > Bmax) continue;
    long shifted = x[i] - Base;
    long b = shifted % 24;
    long a = shifted / 24;/*shifted = 24a + b. b gives the residue block, and a gives the place to insert it.*/
    rclass[b][a]++;
  }
  long ind = 1;/*Which depth we are working at.*/
  long v[4] = {x[0], x[1], x[2], x[3]};/*Initial quadruple.*/
  while (ind > 0) {/*We are coming in trying to swap this circle out.*/
    int cind = ++swaps[ind];/*Increment the swapping index.*/
    if (cind == 4) {/*Overflowed, go back.*/
      swaps[ind] = -1;
      ind--;
      if (!ind) break;
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
    long shifted = newc - Base;
    if (shifted >= 0) {
      long b = shifted % 24;
      long a = shifted / 24;/*shifted = 24a + b. b gives the residue block, and a gives the place to insert it.*/
      rclass[b][a]++;
    }
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
      if (!depthseq) {
        pari_printf("Insufficient memory to reallocate the depth sequence.\n");
        exit(1);
      }
      swaps = pari_realloc(swaps, newdepth * sizeof(int));
      if (!swaps) {
        pari_printf("Insufficient memory to reallocate the swaps.\n");
        exit(1);
      }
      for (i = maxdepth; i < newdepth; i++) swaps[i] = -1;
      primes = pari_realloc(primes, newdepth * sizeof(int));
      if (!primes) {
        pari_printf("Insufficient memory to reallocate the primes.\n");
        exit(1);
      }
      maxdepth = newdepth;
    }
  }
  /*Time to free some of the allocated memory.*/
  pari_free(primes);
  pari_free(swaps);
  pari_free(depthseq);
  /*Save the curvature counts to file.*/
  if (!pari_is_dir("curvcounts")) {
    int s = system("mkdir -p curvcounts");
    if (s == -1) pari_err(e_MISC, "ERROR CREATING DIRECTORY curvcounts");
  }
  char *filestart = stack_sprintf("curvcounts/%Pd_%Pd_%Pd_%Pd_%ld-to-%ld_res-", gel(vorig, 1), gel(vorig, 2), gel(vorig, 3), gel(vorig, 4), Bmin, Bmax);
  for (i = 0; i < nres; i++) {
    long r = res[i];/*The residue.*/
    char *thisfile = stack_sprintf("%s%ld.dat", filestart, r);/*The file name.*/
    FILE *F = fopen(thisfile, "w");/*Created the output file f*/
    long firstind = 0;
    if (Bmin > (Base + r)) firstind = 1;/*Where to start, since we start at Base not Bmin we want to ignore the first one sometimes.*/
    long lastind = classmax - 1, j;
    while ((Base + 24 * lastind + r) > Bmax) lastind--;/*Where to end it.*/
    for (j = firstind; j <= lastind; j++) {
      pari_fprintf(F, "%d\n", rclass[r][j]);
    }
    fclose(F);
  } 
  for (i = 0; i < nres; i++) pari_free(rclass[res[i]]);
  pari_free(rclass);/*The last thing to free*/
}

/*Checks if p is prime but allows for long, not just ulong, and returns 0 if p<=0 (even if -p is prime).*/
int
sisprime(long p)
{
  if (p <= 0) return 0;
  return uisprime(p);
}

/*Finds the number of curvatures between Bmin and Bmin+binsize*nbins-1, saving the counts in blocks of length binsize. Returns [prime counts, thickened counts]. We also save this to a two files.*/
GEN
thickened_bin(GEN v, unsigned long Bmin, unsigned long binsize, unsigned long nbins)
{
  pari_sp av = avma;
  long newv[4];
  long iodd = 0, ieven = 2, i;
  for (i = 1; i <= 4; i++) {/*Make the first two odd.*/
    if (Mod2(gel(v, i))) newv[iodd++] = itos(gel(v, i));
    else newv[ieven++] = itos(gel(v, i));
  }
  if (!sisprime(newv[0]) && !sisprime(newv[1])) pari_err_TYPE("At least one of the odd numbers must be prime", v);
  return gerepilecopy(av, thickened_bin_execute(newv, Bmin, binsize, nbins, v));
}

/*Executes thickened_bin. Assumes the first two entries of v are odd, with at least one positive prime. LEAVES GARBAGE, NOT GEREPILEUPTO SAFE*/
static GEN
thickened_bin_execute(long x[], unsigned long Bmin, unsigned long binsize, unsigned long nbins, GEN vorig)
{
  unsigned long Bmax = Bmin + (nbins * binsize) - 1, i;
  GEN primecounts = const_vecsmall(nbins, 0);
  GEN thickcounts = const_vecsmall(nbins, 0);
  long maxdepth = 200;/*Maximal depth, to start.*/
  long *depthseq = (long *)pari_malloc(maxdepth * sizeof(long));/*depthseq[i] tracks the value we swapped away from in the ith iteration.*/
  if (!depthseq) {
    pari_printf("Insufficient memory to allocate to store the depth sequence.\n");
    exit(1);
  }
  int *swaps = (int *)pari_malloc(maxdepth * sizeof(int));/*Tracks the sequence of swaps, from index 1 to 4.*/
  if (!swaps) {
    pari_printf("Insufficient memory to allocate to store the swaps.\n");
    exit(1);
  }
  int *primes = (int *)pari_calloc(maxdepth * sizeof(int));/*1 if only x[0] is prime, 2 if only x[1] is prime, 3 if both are prime.*/
  if (!primes) {
    pari_printf("Insufficient memory to allocate to store which indices are prime.\n");
    exit(1);
  }
  for (i = 0; i < maxdepth; i++) swaps[i] = -1;/*Initialize to all -1's*/
  for (i = 0; i <= 1; i++) {/*First two, checking primality*/
    if (!sisprime(x[i])) continue;
    primes[0] += (i + 1);
    if (x[i] < Bmin && x[i] > Bmax) continue;/*Ensure we are in the right range.*/
    long binno = 1 + ((x[i] - Bmin) / binsize);
    primecounts[binno]++;
  }
  for (i = 0; i <= 3; i++) {/*Do the first curvatures thickness!.*/
    if (x[i] < Bmin || x[i] > Bmax) continue;
    long binno = 1 + ((x[i] - Bmin) / binsize);
    thickcounts[binno]++;
  }
  long ind = 1;/*Which depth we are working at.*/
  long v[4] = {x[0], x[1], x[2], x[3]};/*Initial quadruple.*/
  while (ind > 0) {/*We are coming in trying to swap this circle out.*/
    int cind = ++swaps[ind];/*Increment the swapping index.*/
    if (cind == 4) {/*Overflowed, go back.*/
      swaps[ind] = -1;
      ind--;
      if (!ind) break;
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
    long binno = newc - Bmin;
    if (binno >= 0) {/*Update this bin in the thickened component.*/
      binno = 1 + (binno / binsize);
      thickcounts[binno]++;
    }
    else binno = 0;
    depthseq[ind] = v[cind];/*Store the value we swapped away from.*/
    v[cind] = newc;/*Update v*/
    switch (cind) {/*Update primes*/
      case 0:/*Check if prime*/
        if (sisprime(newc)) {
          primes[ind] = 3;/*1st index must be prime as we swapped the 0th*/
          if (binno) primecounts[binno]++;
        }
        else primes[ind] = 2;
        break;
      case 1:/*Check if prime*/
        if (sisprime(newc)) {
          primes[ind] = 3;/*0th index must be prime as we swapped the 1st*/
          if (binno) primecounts[binno]++;
        }
        else primes[ind] = 1;
        break;
      default:
        primes[ind] = primes[lastind];/*Odd numbers did not change.*/
    }
    ind++;
    if (ind == maxdepth) {/*We are going too deep, must pari_reallocate the storage location.*/
      long newdepth = maxdepth << 1;/*Double it.*/
      depthseq = pari_realloc(depthseq, newdepth * sizeof(long));
      if (!depthseq) {
        pari_printf("Insufficient memory to reallocate the depth sequence.\n");
        exit(1);
      }
      swaps = pari_realloc(swaps, newdepth * sizeof(int));
      if (!swaps) {
        pari_printf("Insufficient memory to reallocate the swaps.\n");
        exit(1);
      }
      for (i = maxdepth; i < newdepth; i++) swaps[i] = -1;
      primes = pari_realloc(primes, newdepth * sizeof(int));
      if (!primes) {
        pari_printf("Insufficient memory to reallocate the primes.\n");
        exit(1);
      }
      maxdepth = newdepth;
    }
  }
  /*Time to free some of the allocated memory.*/
  pari_free(primes);
  pari_free(swaps);
  pari_free(depthseq);
  /*Save the binned curvature counts to file.*/
  if (!pari_is_dir("curvcounts-binned")) {
    int s = system("mkdir -p curvcounts-binned");
    if (s == -1) pari_err(e_MISC, "ERROR CREATING DIRECTORY curvcounts-binned");
  }
  char *filestart = stack_sprintf("curvcounts-binned/%Pd_%Pd_%Pd_%Pd_from-%lu-size-%lu-nbins-%lu_", gel(vorig, 1), gel(vorig, 2), gel(vorig, 3), gel(vorig, 4), Bmin, binsize, nbins);
  char *filethick = stack_sprintf("%sthick.dat", filestart);
  FILE *Fthick = fopen(filethick, "w");/*Create the output file*/
  for (i = 1; i <= nbins; i++) pari_fprintf(Fthick, "%d\n", thickcounts[i]);
  fclose(Fthick);
  
  char *fileprime = stack_sprintf("%sprime.dat", filestart);
  FILE *Fprime = fopen(fileprime, "w");/*Create the output file*/
  for (i = 1; i <= nbins; i++) pari_fprintf(Fprime, "%d\n", primecounts[i]);
  fclose(Fprime);
  return mkvec2(primecounts, thickcounts);
}


