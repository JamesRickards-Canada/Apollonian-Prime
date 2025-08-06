
/*
Tests the hypothesis for the different genera of discriminant -4a^2, where we look at representations by f(x, y)-a.
We test primes equivalent to l modulo m from 1 to nbins*nperbin: binning the results into nbins with nperbin in each bin.
The return is a vector [g1, g2, ..., gk], where there are k genera. Each gk is a vector of length nbins, whose ith entry is the cumulative count up to nperbin*i of primes represented by f(x, y)-a (not counted with multiplicity: 0 or 1 per prime).
*/
testhyp(a, l, m, nbins, nperbin) = {
  my(top, forms, counts, fl, rep);
  top = nbins * nperbin;
  forms = genera(-4 * a^2);
  counts = vector(#forms);
  for (i = 1, #forms,
    counts[i] = vector(nbins);
    fl = forms[i];
    forprime (p = 2, top,
      if (p % m == l,
        rep = 0;
        for (j = 1, #fl,/*See if it is represented.*/
          if (qfbsolve(fl[j], p + a) != 0, rep = 1; break);
        );
        if (rep, counts[i][ceil(p / nperbin)]++);/*Increment*/
      );
    );
    for (j = 2, nbins, counts[i][j] += counts[i][j - 1]);/*Make it cumulative*/
  );
  return(counts);
}

/*
testhyp, but we do forms of discriminant D and test f(x, y)+a and only care about primes, no l mod m condition.
*/
testhyp2(D, a, nbins, nperbin) = {
  my(top, forms, counts, fl, rep);
  top = nbins * nperbin;
  forms = genera(D);
  counts = vector(#forms);
  for (i = 1, #forms,
    counts[i] = vector(nbins);
    fl = forms[i];
    forprime (p = 2, top,
      rep = 0;
      for (j = 1, #fl,/*See if it is represented.*/
        if (qfbsolve(fl[j], p + a) != 0, rep = 1; break);
      );
      if (rep, counts[i][ceil(p / nperbin)]++);/*Increment*/
    );
    for (j = 2, nbins, counts[i][j] += counts[i][j - 1]);/*Make it cumulative*/
  );
  return(counts);
}



/*Lazy way to get the genera. Returns them as a vector, with each component being a vector of 1 genera.*/
genera(D) = {
  my(cg, prin, ords, gens, last, v, ind, q);
  cg = qfbnarrow(D);
  prin = principalgenus(cg);
  ords = cg[2];
  gens = cg[3];
  last = #ords;
  for (i = 1, last,/*Keep only the generators with even order.*/
    if (ords[i] % 2 == 1, last = i - 1; break);
  );
  if (last == 0, return([prin]));
  topow = vector(last, i, [0, 1]);/*The different genera*/
  v = vector(1 << last);/*How many*/
  ind = 1;
  forvec (X = topow,
    v[ind] = qfbpowmulvec(gens, X);
    ind++;
  );/*Now we have the representatives for each genera*/
  for (i = 1, #v,
    q = v[i];
    v[i] = vector(#prin, j, qfbcomp(q, prin[j]));
  );
  return(v);
}

/*Given the output of qfbnarrow(D), this returns the principal genus. It is poorly done / inefficient, but this isn't really an issue.*/
principalgenus(cg) = {
  my(ords, gens, last, topow, L);
  ords = cg[2];
  gens = cg[3];
  last = #ords;
  for (i = 1, last,/*Square the generators if the order is odd.*/
    if (ords[i] == 2, last = i - 1; break);/*Squares to 1, and everything thereafter does as well.*/
    if (ords[i] % 2 == 1, break);/*This and the rest are odd, and we leave them unchanged.*/
    gens[i] = qfbpow(gens[i], 2);/*Square it*/
    ords[i] = ords[i] >> 1;/*Halve it, still bigger than 1.*/
  );
  if (last == 0, return([qfbpow(gens[1], 2)]));/*All order 2, so only element of principal genus is the identity.*/
  topow = vector(last, i, [0, ords[i] - 1]);/*Possible orders*/
  L = List();/*We know the actual length, but I am lazy.*/
  forvec (X = topow,
    listput(L, qfbpowmulvec(gens, X));
  );
  return(Vec(L));
}

/*gens=[q1, ..., qr] and pows=[a1, ..., as], returns q1^a1*...*qr^ar. We must have s<=r (and can have s<r).*/
qfbpowmulvec(gens, pows) = {
  my(q);
  q = qfbpow(gens[1], pows[1]);
  for (i = 2, #pows, q = qfbcomp(q, qfbpow(gens[i], pows[i])));
  return(q);
}
