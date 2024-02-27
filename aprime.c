/*Methods dealing with prime components of Apollonian circle packings*/

/*INCLUSIONS*/

#include <pari.h>
#include "apol.h"
#include "aprime.h"

/*STATIC DECLARATIONS*/
static GEN primeroots_bin_execute(long x[], unsigned long Bmin, unsigned long binsize, unsigned long nbins, GEN vorig);
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
  if (!pari_is_dir("thickcurvcounts")) {
    int s = system("mkdir -p thickcurvcounts");
    if (s == -1) pari_err(e_MISC, "ERROR CREATING DIRECTORY thickcurvcounts");
  }
  char *filestart = stack_sprintf("thickcurvcounts/%Pd_%Pd_%Pd_%Pd_%ld-to-%ld_res-", gel(vorig, 1), gel(vorig, 2), gel(vorig, 3), gel(vorig, 4), Bmin, Bmax);
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

/*In the full Apollonian circle packing, finds the number of curvatures between Bmin and Bmin+binsize*nbins-1, saving the counts in blocks of length binsize. Returns [prime counts, thickened counts]. We also save this to two files.*/
GEN
primeroots_bin(GEN v, unsigned long Bmin, unsigned long binsize, unsigned long nbins)
{
  pari_sp av = avma;
  GEN vred = apol_red(v, 0, 0);/*precision does not matter.*/
  long newv[4];
  long iodd = 0, ieven = 2, i;
  for (i = 1; i <= 4; i++) {/*Make the first two odd.*/
    if (Mod2(gel(vred, i))) newv[iodd++] = itos(gel(vred, i));
    else newv[ieven++] = itos(gel(vred, i));
  }
  return gerepilecopy(av, primeroots_bin_execute(newv, Bmin, binsize, nbins, v));
}

/*Executes primeroots_bin. Assumes the first two entries of v are odd. LEAVES GARBAGE, NOT GEREPILEUPTO SAFE*/
static GEN
primeroots_bin_execute(long x[], unsigned long Bmin, unsigned long binsize, unsigned long nbins, GEN vorig)
{
  unsigned long Bmax = Bmin + (nbins * binsize) - 1, i;
  GEN primerootcounts = const_vecsmall(nbins, 0);/*Counts quadruples with max entry prime and no other primes*/
  GEN primeallcounts = const_vecsmall(nbins, 0);/*Counts quadruples with max entry prime*/
  GEN allcirclecounts = const_vecsmall(nbins, 0);/*Counts total number of curvatures.*/
  long maxdepth = 200;/*Maximal depth, to start.*/
  long *depthseq = (long *)pari_malloc(maxdepth * sizeof(long));/*depthseq[i] tracks the value we swapped away from in the ith iteration.*/
  int *swaps = (int *)pari_malloc(maxdepth * sizeof(int));/*Tracks the sequence of swaps, from index 1 to 4.*/
  int *primes = (int *)pari_calloc(maxdepth * sizeof(int));/*0 if neither x[0] nor x[1] are prime, 1 if only x[0] is prime, 2 if only x[1] is prime, 3 if both are prime.*/
  for (i = 0; i < maxdepth; i++) swaps[i] = -1;/*Initialize to all -1's*/
  for (i = 0; i <= 1; i++) {/*First two, checking primality*/
    if (!sisprime(x[i])) continue;
    primes[0] += (i + 1);/*We aren't counting the reduced quadruple as a prime root, so just initialize primes[0]/*/
  }
  for (i = 0; i <= 3; i++) {/*First four circles*/
    if (x[i] < Bmin || x[i] > Bmax) continue;
    long binno = 1 + ((x[i] - Bmin) / binsize);
    allcirclecounts[binno]++;
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
    long apbpc = 0;/*Now we can reasonably try a swap.*/
    for (i = 0; i < cind; i++) apbpc += v[i];
    for (i = cind + 1; i < 4; i++) apbpc += v[i];
    long newc = (apbpc << 1) - v[cind];/*2(a+b+c)-d, the new curvature.*/    
    if (newc > Bmax) continue;/*Too big! go back.*/
    long binno = newc - Bmin;
    if (binno >= 0) {/*Update this bin in the thickened component.*/
      binno = 1 + (binno / binsize);
      allcirclecounts[binno]++;
    }
    else binno = 0;
    depthseq[ind] = v[cind];/*Store the value we swapped away from.*/
    v[cind] = newc;/*Update v*/
    switch (cind) {/*Update primes*/
      case 0:/*Check if prime*/
        if (sisprime(newc)) {
          if (primes[lastind] >= 2) primes[ind] = 3;/*Both are prime!*/
          else {
            primes[ind] = 1;/*Just this is prime, and we have a prime root.*/
            if (binno) primerootcounts[binno]++;/*Update the prime root as long as it is above Bmin.*/
          }
          if (binno) primeallcounts[binno]++;/*Always update this.*/
        }
        else {
          if (primes[lastind] >= 2) primes[ind] = 2;/*x[1] is prime only*/
          else primes[ind] = 0;/*Neither are prime.*/
        }
        break;
      case 1:/*Check if prime*/
        if (sisprime(newc)) {
          if (primes[lastind] % 2) primes[ind] = 3;/*Both are prime!*/
          else {
            primes[ind] = 2;/*Just x[1] is prime*/
            if (binno) primerootcounts[binno]++;/*Update the prime root as long as it is above Bmin.*/
          }
          if (binno) primeallcounts[binno]++;/*Always update this.*/
        }
        else {
          if (primes[lastind] % 2) primes[ind] = 1;/*x[0] is prime only*/
          else primes[ind] = 0;/*Neither are prime.*/
        }
        break;
      default:
        primes[ind] = primes[lastind];/*Odd numbers did not change, and we did not get a new prime component root..*/
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
  /*Time to free some of the allocated memory.*/
  pari_free(primes);
  pari_free(swaps);
  pari_free(depthseq);
  /*Save the binned curvature counts to files.*/
  if (!pari_is_dir("fullcurvcounts")) {
    int s = system("mkdir -p fullcurvcounts");
    if (s == -1) pari_err(e_MISC, "ERROR CREATING DIRECTORY curvcounts-binned");
  }
  char *filestart = stack_sprintf("fullcurvcounts/%Pd_%Pd_%Pd_%Pd_from-%lu-size-%lu-nbins-%lu-", gel(vorig, 1), gel(vorig, 2), gel(vorig, 3), gel(vorig, 4), Bmin, binsize, nbins);
  char *fileall = stack_sprintf("%sall.dat", filestart);
  FILE *Fall = fopen(fileall, "w");/*Create the output file*/
  for (i = 1; i <= nbins; i++) pari_fprintf(Fall, "%d\n", allcirclecounts[i]);
  fclose(Fall);
  
  char *fileprime = stack_sprintf("%sprimeroot.dat", filestart);
  FILE *Fprime = fopen(fileprime, "w");/*Create the output file*/
  for (i = 1; i <= nbins; i++) pari_fprintf(Fprime, "%d\n", primerootcounts[i]);
  fclose(Fprime);
  
  char *fileprimeall = stack_sprintf("%sprimeall.dat", filestart);
  FILE *Fprimeall = fopen(fileprimeall, "w");/*Create the output file*/
  for (i = 1; i <= nbins; i++) pari_fprintf(Fprimeall, "%d\n", primeallcounts[i]);
  fclose(Fprimeall);
  return mkvec3(primerootcounts, primeallcounts, allcirclecounts);
}

/*Finds the number of curvatures between Bmin and Bmin+binsize*nbins-1 in the corresponding thickened prime component, saving the counts in blocks of length binsize. Returns [prime counts, thickened counts]. We also save this to two files.*/
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
  if (!pari_is_dir("thickcurvcounts-binned")) {
    int s = system("mkdir -p thickcurvcounts-binned");
    if (s == -1) pari_err(e_MISC, "ERROR CREATING DIRECTORY thickcurvcounts-binned");
  }
  char *filestart = stack_sprintf("thickcurvcounts-binned/%Pd_%Pd_%Pd_%Pd_from-%lu-size-%lu-nbins-%lu_", gel(vorig, 1), gel(vorig, 2), gel(vorig, 3), gel(vorig, 4), Bmin, binsize, nbins);
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

/*Reduces v to the prime root it is a part of, or returns 0 if it does not contain a prime.*/
GEN
primerootred(GEN v)
{
  pari_sp av = avma;
  long i;
  GEN w = cgetg(5, t_VECSMALL);
  for (i = 1; i <= 4; i++) w[i] = itos(gel(v, i));
  long oddinds[2], ind = 0;
  for (i = 1; i <= 4; i++) if (w[i] % 2) oddinds[ind++] = i;
  long primes = 0;/*primes=0 if neither odd is prime, 1 if only the first, 2 if only the second, and 3 if both are prime.*/
  if (sisprime(w[oddinds[0]])) primes++;
  if (sisprime(w[oddinds[1]])) primes += 2;
  if (!primes) return gc_const(av, gen_0);
  for (;;) {
    ind = vecsmall_indexmax(w);/*The largest circle, to swap.*/
    long newc = 0;
    for (i = 1; i < ind; i++) newc += w[i];
    for (i = ind + 1; i <= 4; i++) newc += w[i];
    newc <<= 1;
    newc -= w[ind];/*2(a+b+c)-d*/
    if (newc >= w[ind]) break;/*Done, reduced to the root of the whole packing.*/
    if (ind == oddinds[0]) {/*Swapping out the first odd index.*/
      if (primes == 1) break;/*Swapping only prime, not allowed.*/
      if (sisprime(newc)) {/*new curvature prime*/
        if (primes == 2) primes = 3;/*Only case that needs updating.*/
      }
      else {/*new curvature not prime*/
        if (primes == 3) primes = 2;/*Swapping out one of the two primes.*/
      }
    }
    else if (ind == oddinds[1]) {/*Swapping out second odd index.*/
      if (primes == 2) break;/*Swapping only prime, not allowed.*/
      if (sisprime(newc)) {/*new curvature prime*/
        if (primes == 1) primes = 3;/*Only case that needs updating.*/
      }
      else {/*new curvature not prime*/
        if (primes == 3) primes = 1;/*Swapping out one of the two primes.*/
      }
    }
    w[ind] = newc;
  }
  vecsmall_sort(w);
  return gerepileupto(av, gtovec(w));
}



