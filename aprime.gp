\r apol

print("\n\nType '?aprime' for help.\n\n");
addhelp(aprime, "This package is used for computing thickened prime curvature components of Apollonian circle packings. Installed methods:\n\tprimerootred, primeroots_bin, parthickenedmult, thickened, thickened_bin.");
parigp_version = version();
aprime_library = strprintf("./libaprime-%d-%d-%d.so", parigp_version[1], parigp_version[2], parigp_version[3]);

/*aprime.c*/

install(primeroots_bin,"GUUU",,aprime_library);
addhelp(primeroots_bin,"primeroots_bin(v, Bmin, binsize, nbins): finds the number of (prime) curvatures in the ACP corresponding to v between Bmin and Bmin+binsize*nbins-1, divided into nbins bins of length binsize each. We retun the tuple [prime root counts, prime max entry counts, total counts], and also save this to the files ./fullcurvcounts/v[1]_v[2]_v[3]_v[4]_from-$Bmin-size-$binsize-nbins-$nbins-(primeroot/primeall/all).dat");

install(thickened,"vGG");
addhelp(thickened,"thickened(v, B): finds the multiplicities of all curvatures in the thickened odd prime component of the Apollonian circle packing corresponding to the Descartes quadruple v. If v does not contain an prime curvature, raises an error. If B is an integer we search from 1 to B, else B=[Bmin, Bmax] and we search from Bmin to Bmax. Saves the data to files ./thickcurvcounts/v[1]_v[2]_v[3]_v[4]_Bmin-to-Bmax_res-R.dat, one file per each of the 6/8 residue classes R. The first line of the file corresponds to the multiplicity of the smallest integer equivalent to R modulo 24 that is at least Bmin, and each subsequent line corresponds to the next residue in that class. Be aware that if Bmax-Bmin>=5*10^8, then the files start to be larger than 1GB each.");

install(thickened_bin,"GUUU");
addhelp(thickened_bin,"thickened_bin(v, Bmin, binsize, nbins): finds the number of curvatures between Bmin and Bmin+binsize*nbins-1 in the thickened prime component, saving the counts in blocks of length binsize. Returns [prime counts, thickened counts]. We also save this to two files.");

install(primerootred,"G");
addhelp(primerootred,"primerootred(v): returns the prime root quadruple attached to v, or 0 if v contains no odd primes.");

/*aprime_parallel.c*/
install(parthickenedmult,"GLLLD0,L,");
addhelp(parthickenedmult,"parthickenedmult(v, Bmin, Bmax, Nthreads, {load = 0}: finds the multiplicities of all curvatures in the thickened odd prime component of the Apollonian circle packing corresponding to the Descartes quadruple v. If v does not contain an prime curvature, raises an error. We search in the range Bmin to Bmax (inclusive), and do NOT separate out by residue modulo 24 (so 16/18 of the classes modulo 24 will be 0). Saves the data to file ./thickcurvcounts/v[1]_v[2]_v[3]_v[4]_Bmin-to-Bmax.dat, with one multiplicity per line, starting at Bmin. We use a parallel implementation with Nthreads (total number of threads used is one more due to the main thread that runs everything). If load=1, we also return this as a Vecsmall. This is not suggested to be used if Bmax-Bmin>10^7, as the return vector size starts to be prohibitively large.");
