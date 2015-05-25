@everywhere function loopwinFixsize(nSNP, nReads, rawdata, snpID, bin2geno, nPer,
                                    SNPSize, flagAnnotate, snpWtType, keepIdx,
                                    nPerKeep, kernel, kinship, test,
                                    XtNullBasis, WPreSim, nBlockAscent, nMMmax,
                                    nNullSimPts, tolX, pvalueComputing, device,
                                    nNullSimNewtonIter, rankQPhi, evalPhiAdj,
                                    XPhitNullBasis, yShift, y, KPhiAdj,
                                    weightedW, QPhi, PrePartialSumW,
                                    PreTotalSumW, windowSize, nPreRank, X, ynew)

  QRes = Array(Float64, nPerKeep, rankQPhi);
  nSNPredVec = fill(0, nReads);
  snpSqrtWts = Float64[];
  tmpvec = similar(yShift);
  tmpvecQRes = Array(Float64, rankQPhi);
  yWork = similar(evalPhiAdj);
  partialSumWConst = Array(Float64, nNullSimPts);
  totalSumWConst = Array(Float64, nNullSimPts);
  subXPhitSV = Array(Float64, size(XPhitNullBasis, 2), rankQPhi);
  pSubXPhitSV = pointer(subXPhitSV);
  offset = 0;
  oneConst = ones(nPerKeep);
  geno = Array(Float64, nPer, SNPSize);
  genoOri = Array(Float64, nPer, SNPSize);
  snpIDred = similar(snpID);
  tmpMat = Array(Float64, windowSize, size(XPhitNullBasis, 2));
  VWorkSqrt = Array(Float64, length(evalPhiAdj), windowSize);
  VWorkSqrt2 = Array(Float64, length(evalPhiAdj), windowSize);
  mafOri = Array(Float64, SNPSize);
  maf = Array(Float64, SNPSize);
  sumIsnan = Array(Float64, 1, SNPSize);
  nChrObs = Array(Float64, 1, SNPSize);
  tmpvecKernel = Array(Float64, windowSize);
  if pvalueComputing == "chi2"
    nPtsChi2 = 300;
    partialSumW = Array(Float64, nPtsChi2);
    totalSumW = Array(Float64, nPtsChi2);
    lambda = Array(Float64, 1, nPtsChi2);
    W = Array(Float64, windowSize, nPtsChi2);
    tmpmat0 = Array(Float64, windowSize, nPtsChi2);
    tmpmat1 = Array(Float64, windowSize, nPtsChi2);
    tmpmat2 = Array(Float64, windowSize, nPtsChi2);
    tmpmat3 = Array(Float64, windowSize, nPtsChi2);
    tmpmat4 = Array(Float64, windowSize, nPtsChi2);
    tmpmat5 = Array(Float64, windowSize, nPtsChi2);
    denomvec = Array(Float64, nPtsChi2);
    d1f = Array(Float64, nPtsChi2);
    d2f = Array(Float64, nPtsChi2);
    if test == "eScore"
      simnull = Array(Float64, nPtsChi2);
    else
      simnull = Float64[];
    end
  else
    nPtsChi2 = 300;
    partialSumW = Array(Float64, nNullSimPts);
    totalSumW = Array(Float64, nNullSimPts);
    lambda = Array(Float64, 1, nNullSimPts);
    W = Array(Float64, windowSize, nNullSimPts);
    tmpmat0 = Array(Float64, windowSize, nNullSimPts);
    tmpmat1 = Array(Float64, windowSize, nNullSimPts);
    tmpmat2 = Array(Float64, windowSize, nNullSimPts);
    tmpmat3 = Array(Float64, windowSize, nNullSimPts);
    tmpmat4 = Array(Float64, windowSize, nNullSimPts);
    tmpmat5 = Array(Float64, windowSize, nNullSimPts);
    denomvec = Array(Float64, nNullSimPts);
    d1f = Array(Float64, nNullSimPts);
    d2f = Array(Float64, nNullSimPts);
    if test == "eScore"
      simnull = Array(Float64, nNullSimPts);
    else
      simnull = Float64[];
    end
  end

  results = Array(Any, nReads * iceil(SNPSize / windowSize), 3);
  idx = 0;

  for idxRead = 1 : nReads

    # read the binary data
    if idxRead == nReads
      curSNPSize = nSNP - (idxRead - 1) * SNPSize;
    else
      curSNPSize = SNPSize;
    end

    readgeno!(genoOri, curSNPSize, nPer, SNPSize, idxRead, bin2geno, rawdata,
              offset, flagAnnotate);

    tmpMatIsnan = isnan(genoOri);
    sum!(sumIsnan, tmpMatIsnan);
    for i = 1 : curSNPSize
      nChrObs[1, i] = 2 * (nPer - sumIsnan[1, i]);
    end

    # minor allele frequencies and get rid of mono-allelic SNPs
    counter = 0;
    offset = (idxRead - 1) * SNPSize;
    for i = 1 : curSNPSize
      mafOri[i] = 0.0;
      for j = 1 : nPer
        if !tmpMatIsnan[j, i]
          mafOri[i] += genoOri[j, i];
        end
      end
      mafOri[i] = mafOri[i] / nChrObs[1, i];
      if mafOri[i] != 0
        counter += 1;
        maf[counter] = mafOri[i];
        snpIDred[counter] = snpID[i + offset];
        pGenoOri = pointer(genoOri) + (i - 1) * nPer * sizeof(Float64);
        pGeno = pointer(geno) + (counter - 1) * nPer * sizeof(Float64);
        BLAS.blascopy!(nPer, pGenoOri, 1, pGeno, 1);
      end
    end
    nSNPredVec[idxRead] = curSNPSize - counter;

    if counter == 0
      #println("No valid SNPs at the group which starting from row ",
      #        offset + 1);
      continue;
    end

    # impute missing genotype by expected minor allele counts
    for j = 1 : counter
      for i = 1 : nPer
        if isnan(geno[i, j])
          geno[i, j] = 2 * maf[j];
        end
      end
    end

    # boundary of windows
    if idxRead == nReads
      nSNPUpdate = nSNP - (idxRead - 1) * SNPSize - nSNPredVec[idxRead];
    else
      nSNPUpdate = SNPSize - nSNPredVec[idxRead];
    end
    gLB = [1:windowSize:nSNPUpdate];
    gUB = [gLB[2:end]-1; nSNPUpdate];
    nGrps = length(gLB);
    vc0List = zeros(nGrps);
    vc1List = zeros(nGrps);
    pvalList = zeros(nGrps);

    # looping over windows
    for g = 1 : nGrps

      # get SNP weights
      grpSize = gUB[g] - gLB[g] + 1;
      if isempty(snpWtType)
        # constant weights
        if g > 1 && g < nGrps
          snpSqrtWts = fill!(snpSqrtWts, 1);
        else
          snpSqrtWts = ones(grpSize);
        end
      elseif snpWtType == "beta"
        # weights by beta density evaluated at MAF
        snpSqrtWts = sqrt(pdf(Beta(1, 25), maf[gLB[g] : gUB[g]]));
      elseif snpWtType == "invvar"
        # weights by inverse variance maf*(1-maf)
        snpSqrtWts = 1 ./ sqrt( maf[gLB[g] : gUB[g]] .*
                               (1 - maf[gLB[g] : gUB[g]]) );
      end

      # retrieve group genotypes
      if kernel == "GRM"
        gMAF = maf[gLB[g] : gUB[g]];
        gSNP = geno[keepIdx, gLB[g] : gUB[g]];
        # center and scale genotype matrix by MAF
        BLAS.ger!(-2.0, oneConst, gMAF, gSNP);
        if g == nGrps && grpSize != windowSize
          tmpvecKernel = tmpvecKernel[1 : grpSize];
        end
        for i = 1 : grpSize
          tmpvecKernel[i] = snpSqrtWts[i] / sqrt(2 * gMAF[i] * (1 - gMAF[i]));
        end
        scale!(gSNP, tmpvecKernel);
        #gSNP = gSNP .* (snpSqrtWts ./ sqrt(2 * gMAF .* (1 - gMAF)))';
        if g == nGrps && grpSize != windowSize
          tmpvecKernel = [tmpvecKernel; Array(Float64, windowSize - grpSize)];
        end
      elseif kernel == "IBS1"
        # Table 3 in notes
        gSNP = hcat(geno[keepIdx, gLB[g] : gUB[g]],
                    2 - geno[keepIdx, gLB[g] : gUB[g]]);
        scale!(gSNP, [snpSqrtWts; snpSqrtWts]);
      elseif kernel == "IBS2"
        # Table 2 in notes
        gSNP = hcat(geno[keepIdx, gLB[g] : gUB[g]] .>= 1,
                    geno[keepIdx, gLB[g] : gUB[g]] .<= 1);
        scale!(gSNP, [snpSqrtWts; snpSqrtWts]);
      elseif kernel == "IBS3"
        # Table 1 in notes
        gSNP = hcat(-2 * (geno[keepIdx, gLB[g] : gUB[g]] .<= 1) + 1,
                    2 * (geno[keepIdx, gLB[g] : gUB[g]] .>= 1) - 1,
                    fill(sqrt(2), nPerKeep, gUB[g] - gLB[g] + 1));
        scale!(gSNP, [snpSqrtWts; snpSqrtWts; snpSqrtWts]);
      end

      # testing

      if kinship == "none"

        if test == "eRLRT"
          (b, vc0List[g], vc1List[g], pvalList[g]) =
            vctest(ynew, [], XtNullBasis * gSNP, WPreSim = WPreSim, tests = test,
                   Vform = "half", nBlockAscent = nBlockAscent, nMMmax = nMMmax,
                   nNullSimPts = nNullSimPts, pvalueComputings = pvalueComputing,
                   nNullSimNewtonIter = nNullSimNewtonIter, tolX = tolX,
                   devices = device, PrePartialSumW = PrePartialSumW,
                   PreTotalSumW = PreTotalSumW, partialSumWConst = partialSumWConst,
                   totalSumWConst = totalSumWConst, windowSize = windowSize,
                   partialSumW = partialSumW, totalSumW = totalSumW,
                   lambda = lambda, W = W, nPreRank = nPreRank,
                   tmpmat0 = tmpmat0, tmpmat1 = tmpmat1, tmpmat2 = tmpmat2,
                   tmpmat3 = tmpmat3, tmpmat4 = tmpmat4, tmpmat5 = tmpmat5,
                   denomvec = denomvec, d1f = d1f, d2f = d2f, offset = offset,
                   nPtsChi2 = nPtsChi2, simnull = simnull);
        else
          (b, vc0List[g], vc1List[g], pvalList[g]) =
            vctest(y, X, gSNP, WPreSim = WPreSim, tests = test,
                   Vform = "half", nBlockAscent = nBlockAscent, nMMmax = nMMmax,
                   nNullSimPts = nNullSimPts, pvalueComputings = pvalueComputing,
                   nNullSimNewtonIter = nNullSimNewtonIter, tolX = tolX,
                   devices = device, PrePartialSumW = PrePartialSumW,
                   PreTotalSumW = PreTotalSumW, partialSumWConst = partialSumWConst,
                   totalSumWConst = totalSumWConst, windowSize = windowSize,
                   partialSumW = partialSumW, totalSumW = totalSumW,
                   lambda = lambda, W = W, nPreRank = nPreRank,
                   tmpmat0 = tmpmat0, tmpmat1 = tmpmat1, tmpmat2 = tmpmat2,
                   tmpmat3 = tmpmat3, tmpmat4 = tmpmat4, tmpmat5 = tmpmat5,
                   denomvec = denomvec, d1f = d1f, d2f = d2f, offset = offset,
                   nPtsChi2 = nPtsChi2, simnull = simnull);
        end

      else

        # obtain some orthonormal vectors of space R^n - [QX,QPHI,S]
        # See Golub and Van Loan (1996) Algorithm 12.4.2 on p602
        if g == nGrps && grpSize != windowSize
          tmpMat = tmpMat[1 : grpSize, :];
        end
        BLAS.gemm!('T', 'N', 1.0, gSNP, XPhitNullBasis, 0.0, tmpMat);
        (UXPhitS, svalXPhitS, VXPhitS) = svd(tmpMat, thin = false);
        if g == nGrps && grpSize != windowSize
          tmpMat = [tmpMat; Array(Float64, windowSize - grpSize,
                                  size(XPhitNullBasis, 2))];
        end
        rankXPhitS =
          countnz(svalXPhitS .>
                  max(size(XPhitNullBasis, 2),
                      size(gSNP, 2) * eps(svalXPhitS[1])));
        pXPhitSV = pointer(VXPhitS) +
          rankXPhitS * size(XPhitNullBasis, 2) * sizeof(Float64);
        BLAS.blascopy!(size(XPhitNullBasis, 2) * rankQPhi,
                       pXPhitSV, 1, pSubXPhitSV, 1);
        BLAS.gemm!('N', 'N', 1.0, XPhitNullBasis, subXPhitSV, 0.0, QRes);

        # working Y and working non-trivial variance component
        #yWork = (PhiAdjsvd[:U]' * (yShift + KPhiAdj * (QRes' * y))) ./ sqrt(evalPhiAdj);
        BLAS.blascopy!(length(yShift), yShift, 1, tmpvec, 1);
        BLAS.gemv!('T', 1.0, QRes, y, 0.0, tmpvecQRes);
        BLAS.gemv!('N', 1.0, KPhiAdj, tmpvecQRes, 1.0, tmpvec);
        #yWork = BLAS.gemv('T', 1.0, PhiAdjsvd[:U], tmpvec) ./ sqrt(evalPhiAdj);
        BLAS.gemv!('T', 1.0, weightedW, tmpvec, 0.0, yWork);
        #for i = 1 : length(yWork)
        #  yWork[i] = yWork[i] * InvSqrtevalPhiAdj[i];
        #end
        #VWorkSqrt = PhiAdjsvd[:U]' * (QPhi' * gSNP) ./ sqrt(evalPhiAdj);
        if g == nGrps && grpSize != windowSize
          VWorkSqrt = VWorkSqrt[:, 1 : grpSize];
          VWorkSqrt2 = VWorkSqrt2[:, 1 : grpSize];
        end
        BLAS.gemm!('T', 'N', 1.0, QPhi, gSNP, 0.0, VWorkSqrt2);
        BLAS.gemm!('T', 'N', 1.0, weightedW, VWorkSqrt2, 0.0, VWorkSqrt);
        #scale!(InvSqrtevalPhiAdj, VWorkSqrt);

        # obtain LRT and p-value
        (b, vc0List[g], vc1List[g], pvalList[g]) =
          vctest(yWork, [], VWorkSqrt, WPreSim = WPreSim, tests = test,
                 Vform = "half", nBlockAscent = nBlockAscent, nMMmax = nMMmax,
                 nNullSimPts = nNullSimPts, pvalueComputings = pvalueComputing,
                 nNullSimNewtonIter = nNullSimNewtonIter, tolX = tolX,
                 devices = device, PrePartialSumW = PrePartialSumW,
                 PreTotalSumW = PreTotalSumW, partialSumWConst = partialSumWConst,
                 totalSumWConst = totalSumWConst, windowSize = windowSize,
                 partialSumW = partialSumW, totalSumW = totalSumW,
                 lambda = lambda, W = W, nPreRank = nPreRank,
                 tmpmat0 = tmpmat0, tmpmat1 = tmpmat1, tmpmat2 = tmpmat2,
                 tmpmat3 = tmpmat3, tmpmat4 = tmpmat4, tmpmat5 = tmpmat5,
                 denomvec = denomvec, d1f = d1f, d2f = d2f, offset = offset,
                 nPtsChi2 = nPtsChi2, simnull = simnull);

        # recover VWorkSqrt
        if g == nGrps && grpSize != windowSize
          VWorkSqrt = hcat(VWorkSqrt,
                           Array(Float64, length(evalPhiAdj),
                                 windowSize - grpSize));
          VWorkSqrt2 = hcat(VWorkSqrt2,
                            Array(Float64, length(evalPhiAdj),
                                  windowSize - grpSize));
        end

      end

      idx += 1;
      results[idx, :] = hcat(snpIDred[gLB[g]], snpIDred[gUB[g]], pvalList[g]);

    end

  end

  return results[1:idx, :];

end