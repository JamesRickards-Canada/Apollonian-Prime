\r apol

print("\n\nType '?aprime' for help.\n\n");
addhelp(aprime, "This package is used for computing thickened prime curvature components of Apollonian circle packings. Installed methods:\n\tthickened, thickened_bin.");
parigp_version = version();
aprime_library = strprintf("./libaprime-%d-%d-%d.so", parigp_version[1], parigp_version[2], parigp_version[3]);

/*aprime.c*/

install(thickened,"vGG",,aprime_library);
addhelp(thickened,"thickened(v, B): finds the multiplicities of all curvatures in the thickened odd prime component of the Apollonian circle packing corresponding to the Descartes quadruple v. If v does not contain an prime curvature, raises an error. If B is an integer we search from 1 to B, else B=[Bmin, Bmax] and we search from Bmin to Bmax. Saves the data to files ./curvcounts/v[1]_v[2]_v[3]_v[4]_Bmin-to-Bmax_res-R.dat, one file per each of the 6/8 residue classes R. The first line of the file corresponds to the multiplicity of the smallest integer equivalent to R modulo 24 that is at least Bmin, and each subsequent line corresponds to the next residue in that class. Be aware that if Bmax-Bmin>=5*10^8, then the files start to be larger than 1GB each.");

install(thickened_bin,"GUUU");
addhelp(thickened_bin,"thickened_bin(v, Bmin, binsize, nbins): Finds the number of curvatures between Bmin and Bmin+binsize*nbins-1 in the thickened prime component, saving the counts in blocks of length binsize. Returns [prime counts, thickened counts]. We also save this to two files.");

