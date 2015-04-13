@doc """
gwasvctest(plinkFile::String = "", covFile::String = "", traitFile::String = "",
annotationFile::String = "", outFile::String = "", test::String = "eRLRT",
kinship::String = "GRM", kernel::String = "GRM", pvalueComputing::String = "chi2",
windowSize::Int = 50, infLambda::Float64 = 0.0, MemoryLimit::Int = 200000000,
nNullSimPts::Int = 10000, nMMmax::Int = 0, nBlockAscent::Int = 1000,
nNullSimNewtonIter::Int = 15, tolX::Float64 = 1e-6, snpWtType::String = "")

Test the SNP set effect with presence of additive genetic effect and
enviromental effect.

Optional input name-value pairs:

"plinkFile" - Input Plink files names without extension.

"covFile" - Input file to form the covariate matrix X.

"traitFile" - Input file to form the response vector y.

"annotationFile" - Input file for the annotation of SNPs.

"test"- "eLRT"|"eRLRT"|"eScore"|"none", the requested test, it also
dictates the estimation method, default is "eRLRT".

"kinship" - "GRM"|"MoM"|"theoretical"|"none", indicates method used to
compute the kinship matrix, default is "GRM". You can also use your
own kinship matrix, you can do that by replacing "none" with the file
name in which store the pre-calculated kinship matrix.

"kernel" - "GRM"|"IBS1"|"IBS2"|"IBS3", indicates method used to compute
the kernel matrix S, default is "GRM".

"pvalueComputing" - "MonteCarlo"|"chi2", indicates the method (Monte
Carlo simulaion or Chi square approximation) used to compute the
p-value for eLRT or eRLRT, default is "chi2".

"windowSize" - # SNPs in a SNP set, default is 50.

"infLambda" - Any rational number with float precision. Control the
trade-off for approximating kinship matrix, default is 0.0.

"MemoryLimit" - Max memory used to read and store the raw data, default
is 200000000 (equivalent to 200 megabytes).

"nNullSimPts"- # simulation samples, default is 10000.

"nMMmax" - Max MM iterations, default is 10 (eLRT) or 1000 (eRLRT).

"nBlockAscent" - Max block ascent iterations, default is 1000.

"nNullSimNewtonIter" - Max Newton iteration, default is 15.

"tolX" - Tolerance in change of parameters, default is 1e-6.

"snpWtType" - "beta"|"invvar", add weights to the kernel matrix S,
default is "" (add no weights by default).

Output name-value pair:

"outFile" - Output file which store the P-values for SNP sets effects.
  """ ->
function gwasvctest(args...; covFile::String = "", device::String = "CPU",
                    infLambda::Float64 = 0.0, kernel::String = "GRM",
                    kinship::String = "GRM", nMMmax::Int = 0,
                    nBlockAscent::Int = 1000, nNullSimPts::Int = 10000,
                    nNullSimNewtonIter::Int = 15, outFile::String = "",
                    plinkFile::String = "", snpWtType::String = "",
                    test::String = "eRLRT", tolX::Float64 = 1e-6,
                    traitFile::String = "", MemoryLimit::Int = 200000000,
                    pvalueComputing::String = "chi2", windowSize::Int = 50,
                    annotationFile::String = "")

  #blas_set_num_threads(4);

  ## parse data from files
  plinkBedfile = string(plinkFile, ".bed");
  plinkBimfile = string(plinkFile, ".bim");
  plinkFamfile = string(plinkFile, ".fam");

  # BIM file: chr, rs#, morgan, bp position, allele1, allele2
  bimdata = readdlm(plinkBimfile);
  snpID = bimdata[:, 2];
  snpID = convert(Array{String, 1}, snpID);
  nSNP = length(snpID);

  # FAM file: fam ID, ind ID, father ID, mother ID, sex, phenotype
  famdata = readdlm(plinkFamfile);
  nPer = length(famdata[:, 1]);

  # prepare for reading binary data
  map2geno = [2; NaN; 1; 0];
  bin2geno = Array(Float64, 4, 2 ^ 8);
  countbin = 1;
  for iB4 = 0 : 3
    tmpper4 = map2geno[iB4 + 1];
    for iB3 = 0 : 3
      tmpper3 = map2geno[iB3 + 1];
      for iB2 = 0 : 3
        tmpper2 = map2geno[iB2 + 1];
        for iB1 = 0 : 3
          tmpper1 = map2geno[iB1 + 1];
          bin2geno[1, countbin] = tmpper1;
          bin2geno[2, countbin] = tmpper2;
          bin2geno[3, countbin] = tmpper3;
          bin2geno[4, countbin] = tmpper4;
          countbin += 1;
        end
      end
    end
  end

  # BED file:
  fid = open(plinkBedfile);
  bedheader = mmap_array(Uint8, (3, 1), fid);
  # decide BED file version and snp/individual-major
  # PLINK coding (bits->genotype): 00->1/1, 01->1/2, 10->Missing, 11->2/2
  # We use minor allele counts (0,1,2,nan) coding here
  # geno = 2 => NaN
  # geno = 1 => 1
  # geno = 0 => 2
  # geno = 3 => 0
  if bits(bedheader[1, 1]) == "01101100" && bits(bedheader[2, 1]) == "00011011"
    # v1.0 BED file
    if bits(bedheader[3, 1]) == "00000001"
      # SNP-major
      #rawdata = mmap_array(Uint8, (iceil(nPer / 4), nSNP), fid, 3);
      rawdata = mmap_array(Int8, (iceil(nPer / 4), nSNP), fid, 3);
      SNPSize = ifloor(MemoryLimit / (32 * iceil(nPer / 4)));
    else
      # individual-major
      error("gwasvctest:bedwrongn\n",
            "# Individual-major BED file found!",
            "Please transform to SNP-major BED file using PLINK");
    end
  else
    # v0.99 BED file: individual-major
    error("gwasvctest:bedwrongn\n",
          "# Individual-major BED file found!",
          "Please transform to SNP-major BED file using PLINK");
  end
  close(fid);
  # number of read times
  nReads = iceil(nSNP / SNPSize);

  # check if provide SNP annotation file
  if isempty(annotationFile)
    flagAnnotate = false;
  else
    flagAnnotate = true;
    grpInfo = zeros(Int64, nSNP);
    offsetSize = zeros(Int64, nSNP);
    readAnnotate!(annotationFile, snpID, nSNP, grpInfo, offsetSize);
    grpInfo = grpInfo[grpInfo .> 0];
    nGrp = length(grpInfo);
    offsetSize = offsetSize[1 : nGrp];
  end


  ## determine the additive genetic component

  if kinship == "GRM" || kinship == "MoM"

    genoOri = Array(Float64, nPer, SNPSize);
    geno = Array(Float64, nPer, SNPSize);
    kinMat = zeros(nPer, nPer);
    tmpscalar = 0.0;
    nSNPred = 0;
    offset = 0;
    sumIsnan = Array(Float64, 1, SNPSize);
    nChrObs = Array(Float64, 1, SNPSize);
    mafOri = Array(Float64, SNPSize);
    maf = Array(Float64, SNPSize);
    oneConst = ones(nPer);
    tmpvecKernel = Array(Float64, SNPSize);

    for idxRead = 1 : nReads

      # read the binary data
      if idxRead == nReads
        curSNPSize = nSNP - (idxRead - 1) * SNPSize;
      else
        curSNPSize = SNPSize;
      end

      readgeno!(genoOri, curSNPSize, nPer, SNPSize, idxRead, bin2geno, rawdata,
                offset, false);

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
          pGenoOri = pointer(genoOri) + (i - 1) * nPer * sizeof(Float64);
          pGeno = pointer(geno) + (counter - 1) * nPer * sizeof(Float64);
          BLAS.blascopy!(nPer, pGenoOri, 1, pGeno, 1);
        end
      end
      nSNPred += (curSNPSize - counter);

      # impute missing genotype by expected minor allele counts
      for j = 1 : counter
        for i = 1 : nPer
          if isnan(geno[i, j])
            geno[i, j] = 2 * maf[j];
          end
        end
      end

      if kinship == "GRM"
        # center and scale genotype matrix by minor allele frequencies
        gMAF = maf[1 : counter];
        if counter < SNPSize
          geno = geno[:, 1 : counter];
          tmpvecKernel = tmpvecKernel[1 : counter];
        end
        BLAS.ger!(-2.0, oneConst, gMAF, geno);
        for i = 1 : counter
          tmpvecKernel[i] = 1.0 / sqrt(2 * gMAF[i] * (1 - gMAF[i]));
        end
        scale!(geno, tmpvecKernel);
        # estimate kinship by genetic relation matrix (GRM) from dense markers
        # TODO: Check Yang et al paper the precise formula
        BLAS.gemm!('N', 'T', 1.0, geno, geno, 1.0, kinMat);
        if counter != SNPSize
          geno = [geno Array(Float64, nPer, SNPSize - counter)];
          tmpvecKernel = [tmpvecKernel; Array(Float64, SNPSize - counter)];
        end
      elseif kinship == "MoM"
        if counter < SNPSize
          geno = geno[:, 1 : counter];
        end
        # shift to {-1,0,1} encoding
        geno = geno - 1.0;
        # estimate kinship by method of moment from dense markers
        gMAF = maf[1 : counter];
        tmpscalar += sumabs2(gMAF) + sumabs2(1 - gMAF);
        BLAS.gemm!('N', 'T', 0.5, geno, geno, 1.0, kinMat);
        if counter != SNPSize
          geno = [geno Array(Float64, nPer, SNPSize - counter)];
        end
      end

    end

    if kinship == "GRM"
      scale!(kinMat, 1.0 / (nSNP - nSNPred));
    elseif kinship == "MoM"
      kinMat += (nSNP - nSNPred) / 2.0 - tmpscalar;
      scale!(kinMat, 1.0 / (nSNP - nSNPred - tmpscalar));
    end

  elseif kinship == "theoretical"

    # use theoretical kinship computed from pedigree information
    # use the kinship function provided by Lange
    IndID = famdata[:, 2];
    AugIndID = unique([famdata[:, 2]; famdata[:, 3]; famdata[:, 4]]);
    AugIndID = AugIndID[AugIndID .!= 0];
    AugIndID = sort(AugIndID);
    keepIdx = findin(AugIndID, IndID);
    vals = [1 : length(AugIndID)];
    IndDict = Dict(AugIndID, vals);
    faID = famdata[:, 3];
    moID = famdata[:, 4];
    father = zeros(Int64, length(AugIndID));
    mother = zeros(Int64, length(AugIndID));
    for i = 1 : nPer
      if faID[i] != 0
        father[keepIdx[i]] = IndDict[faID[i]];
        mother[keepIdx[i]] = IndDict[moID[i]];
      end
    end
    kinMat = kinshipcoef(father, mother, 1, length(AugIndID));
    kinMat = kinMat[keepIdx, keepIdx];

  elseif kinship != "none"

    # user provided kinship
    kinMat = readdlm(kinship, ',');
    if size(kinMat, 1) != nPer || size(kinMat, 2) != nPer
      error("gwasvctest:kinshipwrongn\n",
            "# individuals in kinship file does not match plink files");
    end

  end


  ## get responses
  if isempty(traitFile)
    # no trait file provided: retrieve from plink files
    y = famdata[:, 6];
    y = convert(Array{Float64, 1}, y);
  else
    # user supplied trait values (only take the first column) in Plink
    # phenotype file (no headerline) format: FamID IndID Trait1 Trait2 ...
    # Only Trait1 will be read
    y = readdlm(traitFile);
    if size(y, 1) != nPer
      error("gwasvctest:ywrongn\n",
            "# individuals in trait file does not match plink files");
    end
    y = y[:, 3];
    y = convert(Array{Float64, 1}, y);
  end

  ## get covariates
  if isempty(covFile)
    # no covariates provided
    X = ones(nPer, 1);
  else
    # user provided covariates (intercept if exists has to be provided)
    X = readdlm(covFile);
    if size(X, 1) != nPer
      error("gwasvctest:covwrongn\n",
            "# individuals in covariate file does not match plink files");
    end
    X = [ones(nPer, 1) X[:, 3:end]];
    X = convert(Array{Float64, 2}, X);
  end

  ## individuals with misisng phenotype are removed from analysis
  # missing covariate values are imputed by avarages
  # Possible solution for speed up: remove the individuals with missing
  # phenotype in advance

  keepIdx = !isnan(y);
  nPerKeep = countnz(keepIdx);
  y = y[keepIdx];
  X = X[keepIdx, :];
  meanX = zeros(size(X, 2));
  for i = 1:length(meanX)
    meanX[i] = mean(X[:, i][!isnan(X[:, i])]);
  end
  ridx, cidx = ind2sub(size(X), find(isnan(X)));
  X[sub2ind(size(X), find(isnan(X)))] = meanX[cidx];
  kinMat = kinMat[keepIdx, keepIdx];


  ## pre-processing for testing

  # if provide annotation file, avoid large gene set
  if flagAnnotate
    numExceed = countnz(grpInfo .> ifloor(nPerKeep / 4));
    for i = 1 : numExceed
      tmpidx = find(grpInfo .> ifloor(nPerKeep / 4))[1];
      innerItr = iceil(grpInfo[tmpidx] / windowSize);
      tmpgrpInfo = zeros(Int64, innerItr);
      tmpoffsetSize = zeros(Int64, innerItr);
      for j = 1 : innerItr
        if j < innerItr
          tmpgrpInfo[j] = windowSize;
          tmpoffsetSize[j] = offsetSize[tmpidx] + (j - 1) * windowSize;
        else
          tmpgrpInfo[j] = grpInfo[tmpidx] - (j - 1) * windowSize;
          tmpoffsetSize[j] = offsetSize[tmpidx] + (j - 1) * windowSize;
        end
      end
      grpInfo = [grpInfo[1:tmpidx-1]; tmpgrpInfo; grpInfo[tmpidx+1:end]];
      offsetSize = [offsetSize[1:tmpidx-1]; tmpoffsetSize; offsetSize[tmpidx+1:end]];
    end

    windowSize = maximum(grpInfo);
    nReads = length(grpInfo);
  end

  # simulate standard normals and square it
  WPreSim = randn(windowSize, nNullSimPts);
  for idxWc = 1 : nNullSimPts
    for idxWr = 1 : windowSize
      WPreSim[idxWr, idxWc] = WPreSim[idxWr, idxWc] * WPreSim[idxWr, idxWc];
    end
  end

  if kinship == "none"
    # only one non-trivial variance component: SNP set

    # pre-compute basis of N(X') and projected y
    if test == "eRLRT"
      Xsvd = svdfact(X, thin = false);
      rankX = countnz(Xsvd[:S] .> nPerKeep * eps(Xsvd[:S][1]));
      XtNullBasis = Xsvd[:U][:, rankX + 1 : end]';
      ynew = y;
      BLAS.gemv!('N', 1.0, XtNullBasis, y, 0.0, ynew);
    end

  else
    # >=2 non-trivial variance components

    # obtain a basis of X and (possibly) extra variance components
    # Xsvd[:U] -- QX
    Xsvd = svdfact(X, thin = false);
    rankX = countnz(Xsvd[:S] .> nPerKeep * eps(Xsvd[:S][1]));
    XtNullBasis = Xsvd[:U][:, rankX + 1 : end];
    QX = Xsvd[:U][:, 1 : rankX];

    # eigen-decomposition of the additive component
    # Phieig[:vectors] -- evecPhi
    # Phieig[:values] -- evalPhi
    Phieig = eigfact(kinMat);
    rankPhi = countnz(Phieig[:values] .>
                      max(infLambda, nPerKeep * eps(maximum(Phieig[:values]))));
    #evecPhi = Phieig[:vectors];
    #evalPhi = Phieig[:values];
    # enforce first entry of each eigenvector to be >=0
    # for cross-platform reproducibility
    idxEvecPhi = vec(Phieig[:vectors][1, :]) .< 0;
    #Phieig[:vectors][:, idxEvecPhi] = - Phieig[:vectors][:, idxEvecPhi];
    scale!(Phieig[:vectors][:, idxEvecPhi], -1.0);
    #evecPhi[:, idxEvecPhi] = - evecPhi[:, idxEvecPhi];

    # determine a rank for additive component and a rank for snp set S
    if kernel == "GRM"
      rankS = windowSize;
    elseif kernel == "IBS1"
      rankS = windowSize + 1;
    elseif kernel == "IBS2"
      rankS = 2 * windowSize;
    elseif kernel == "IBS3"
      rankS = 2 * windowSize + 1;
    end
    rankPhi = min(rankPhi, ifloor((nPerKeep - rankX - rankS) / 2));

    # low rank approximation
    sortIdx = sortperm(Phieig[:values], rev = true);
    evalPhi = Phieig[:values][sortIdx[1 : rankPhi]];
    evecPhi = Phieig[:vectors][:, sortIdx[1 : rankPhi]];

    # obtain orthonormal basis of (additive component - covariates)
    # See Golub and Van Loan (1996) Algorithm 12.4.3 & Theorem 12.4.2 on p604
    QPhisvd = svdfact(BLAS.gemm('T', 'N', XtNullBasis, evecPhi));
    PhiKeepIdx = (1 - QPhisvd[:S]) .< 1e-6;
    rankQPhi = countnz(PhiKeepIdx);
    #QPhi = evecPhi * QPhisvd[:V][:, 1 : rankQPhi];
    QPhi = Array(Float64, nPerKeep, rankQPhi);
    BLAS.gemm!('N', 'N', 1.0, evecPhi, QPhisvd[:V][:, 1 : rankQPhi], 0.0, QPhi);
    XPhitNullBasis = null([QX evecPhi]');

    # obtain the square root of the rotated Phi
    # PhiAdjsvd[:U] -- W
    # PhiAdjsvd[:S] -- evalPhiAdj
    tmpMat = Array(Float64, rankQPhi, rankPhi);
    #println("rankPhi: ", rankPhi);
    #println("col of evecPhi: ", size(evecPhi, 2));
    BLAS.gemm!('T', 'N', 1.0, QPhi, evecPhi, 0.0, tmpMat);
    scale!(tmpMat, sqrt(evalPhi));
    PhiAdjsvd = svdfact(tmpMat, thin = false);
    #PhiAdjsvd = svdfact(BLAS.gemm('T', 'N', QPhi, evecPhi) .*
    #                    sqrt(evalPhi)', thin = false);
    #W = PhiAdjsvd[:U];
    evalPhiAdj = PhiAdjsvd[:S] .^ 2;
    # enforce first entry of each eigenvector to be >=0
    idxW = vec(PhiAdjsvd[:U][1, :]) .< 0;
    #PhiAdjsvd[:U][:, idxW] = - PhiAdjsvd[:U][:, idxW];
    scale!(PhiAdjsvd[:U][:, idxW], -1.0);
    #KPhiAdj = PhiAdjsvd[:U] .* sqrt(evalPhiAdj' / minimum(evalPhiAdj) - 1);
    KPhiAdj = PhiAdjsvd[:U][:, :];
    scale!(KPhiAdj, sqrt(evalPhiAdj / minimum(evalPhiAdj) - 1));
    InvSqrtevalPhiAdj = similar(evalPhiAdj);
    for i = 1 : length(evalPhiAdj)
      InvSqrtevalPhiAdj[i] = 1 / sqrt(evalPhiAdj[i]);
    end
    weightedW = PhiAdjsvd[:U][:, :];
    scale!(weightedW, InvSqrtevalPhiAdj);

    # precompute shift in Y
    yShift = QPhi' * y;

    # simulate Chi Squares and sums of Chi Squares
    tmpwinSize = min(length(evalPhiAdj), windowSize);
    nPreRank = tmpwinSize;
    PrePartialSumW = Array(Float64, nNullSimPts, nPreRank - 1);
    PreTotalSumW = Array(Float64, nNullSimPts, nPreRank - 1);
    tmpSumVec = Array(Float64, nNullSimPts);
    p1 = pointer(PrePartialSumW);
    BLAS.blascopy!(nNullSimPts,
                   rand(Chisq(length(evalPhiAdj) - tmpwinSize + nPreRank - 1),
                        nNullSimPts),
                   1, p1, 1);
    for i = 1 : nNullSimPts
      pW = pointer(WPreSim) + (i - 1) * windowSize * sizeof(Float64);
      tmpSumVec[i] = BLAS.asum(windowSize - nPreRank + 1, pW, 1);
      PreTotalSumW[i, 1] = tmpSumVec[i] + PrePartialSumW[i, 1];
    end
    for j = 2 : nPreRank - 1
      p1 = pointer(PrePartialSumW) + (j - 1) * nNullSimPts * sizeof(Float64);
      BLAS.blascopy!(nNullSimPts,
                     rand(Chisq(length(evalPhiAdj) - tmpwinSize + nPreRank - j),
                          nNullSimPts),
                     1, p1, 1);
      for i = 1 : nNullSimPts
        tmpSumVec[i] += WPreSim[windowSize - nPreRank + j, i];
        PreTotalSumW[i, j] = tmpSumVec[i] + PrePartialSumW[i, j];
      end
    end

  end

  # prepare output
  if isempty(outFile)
    outFile = string(plinkFile, "-julia.out");
  end
  fid = open(outFile, "w");
  println(fid, "StartSNP,EndSNP,pvalue");
  close(fid);


  ## test group by group

  QRes = Array(Float64, nPerKeep, rankQPhi);
  nSNPredVec = fill(0, nReads);
  nActualSNP = 0;
  snpSqrtWts = Float64[];
  gMAF = Float64[];
  gSNP = Float64[];
  tmpvec = similar(yShift);
  tmpvecQRes = Array(Float64, rankQPhi);
  yWork = similar(evalPhiAdj);
  partialSumWConst = Array(Float64, nNullSimPts);
  totalSumWConst = Array(Float64, nNullSimPts);
  subXPhitSV = Array(Float64, size(XPhitNullBasis, 2), rankQPhi);
  pSubXPhitSV = pointer(subXPhitSV);
  offset = 0;
  oneConst = ones(nPerKeep);
  if !flagAnnotate
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
  else
    geno = Array(Float64, nPer, windowSize);
    genoOri = Array(Float64, nPer, windowSize);
    snpIDred = Array(String, windowSize);
    mafOri = Array(Float64, windowSize);
    maf = Array(Float64, windowSize);
    sumIsnan = Array(Float64, 1, windowSize);
    nChrObs = Array(Float64, 1, windowSize);
  end
  if pvalueComputing == "chi2"
    partialSumW = Array(Float64, 300);
    totalSumW = Array(Float64, 300);
    lambda = Array(Float64, 1, 300);
    W = Array(Float64, windowSize, 300);
    tmpmat0 = Array(Float64, windowSize, 300);
    tmpmat1 = Array(Float64, windowSize, 300);
    tmpmat2 = Array(Float64, windowSize, 300);
    tmpmat3 = Array(Float64, windowSize, 300);
    tmpmat4 = Array(Float64, windowSize, 300);
    tmpmat5 = Array(Float64, windowSize, 300);
    denomvec = Array(Float64, 300);
    d1f = Array(Float64, 300);
    d2f = Array(Float64, 300);
  else
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
  end

  ## no annotation file
  if !flagAnnotate

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
      nActualSNP += counter;

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
      testStats = Stats();

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
            (b, vc0List[g], vc1List[g], testStats) =
              vctest(ynew, [], XtNullBasis * gSNP, WPreSim = WPreSim, tests = test,
                     Vform = "half", nBlockAscent = nBlockAscent, nMMmax = nMMmax,
                     nNullSimPts = nNullSimPts, pvalueComputings = pvalueComputing,
                     nNullSimNewtonIter = nNullSimNewtonIter, tolX = tolX,
                     devices = device);
          else
            (b, vc0List[g], vc1List[g], testStats) =
              vctest(y, X, gSNP, WPreSim = WPreSim, tests = test,
                     Vform = "half", nBlockAscent = nBlockAscent, nMMmax = nMMmax,
                     nNullSimPts = nNullSimPts, pvalueComputings = pvalueComputing,
                     nNullSimNewtonIter = nNullSimNewtonIter, tolX = tolX,
                     devices = device);
          end
          pvalList[g] = testStats.vc1_pvalue;

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
          (b, vc0List[g], vc1List[g], testStats) =
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
                   denomvec = denomvec, d1f = d1f, d2f = d2f);
          pvalList[g] = testStats.vc1_pvalue;

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

      end

      # write output
      fid = open(outFile, "a");
      for g = 1 : nGrps
        println(fid, snpIDred[gLB[g]], ",", snpIDred[gUB[g]], ",", pvalList[g]);
      end
      close(fid);

      # display progress
      if nReads > 1
        if idxRead < nReads
          println(SNPSize * idxRead / nSNP * 100, " %");
        else
          println("100 %");
        end
      end

    end

    # summary
    println("There are ", nSNP, " SNPs in total.");
    println(nActualSNP, " SNPs were analysed.");

    ## provide annotation file
  else

    # define some variables
    vc0List = 0.0;
    vc1List = 0.0;
    pvalList = 0.0;
    testStats = Stats();
    if nReads > 1
      percentVec = zeros(Int64, 10);
      for i = 1 : 10
        percentVec[i] = iceil(nReads * i / 10);
      end
      percentIdx = 1;
    end

    # loop over genes
    for idxRead = 1 : nReads

      # read the binary data
      curSNPSize = grpInfo[idxRead];
      readgeno!(genoOri, curSNPSize, nPer, SNPSize, idxRead, bin2geno, rawdata,
                offset, flagAnnotate, offsetSize = offsetSize);

      tmpMatIsnan = isnan(genoOri);
      sum!(sumIsnan, tmpMatIsnan);
      for i = 1 : curSNPSize
        nChrObs[1, i] = 2 * (nPer - sumIsnan[1, i]);
      end

      # minor allele frequencies and get rid of mono-allelic SNPs
      counter = 0;
      offset = offsetSize[idxRead];
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
      nActualSNP += counter;
      if counter == 0
        #println("No valid SNPs at the group which starting from row ",
        #        offsetSize[idxRead] + 1);
        # display progress
        if nReads > 1
          if idxRead == percentVec[percentIdx]
            println(percentIdx, "0 %");
            percentIdx += 1;
          end
        end
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

      # get SNP weights
      grpSize = counter;
      if isempty(snpWtType)
        # constant weights
        snpSqrtWts = ones(grpSize);
      elseif snpWtType == "beta"
        # weights by beta density evaluated at MAF
        snpSqrtWts = sqrt(pdf(Beta(1, 25), maf[1 : grpSize]));
      elseif snpWtType == "invvar"
        # weights by inverse variance maf*(1-maf)
        snpSqrtWts = 1 ./ sqrt( maf[1 : grpSize] .*
                               (1 - maf[1 : grpSize]) );
      end

      # retrieve group genotypes
      if kernel == "GRM"
        gMAF = maf[1 : grpSize];
        gSNP = geno[keepIdx, 1 : grpSize];
        # center and scale genotype matrix by MAF
        BLAS.ger!(-2.0, oneConst, gMAF, gSNP);
        tmpvecKernel = zeros(grpSize);
        for i = 1 : grpSize
          tmpvecKernel[i] = snpSqrtWts[i] / sqrt(2 * gMAF[i] * (1 - gMAF[i]));
        end
        scale!(gSNP, tmpvecKernel);
        #gSNP = gSNP .* (snpSqrtWts ./ sqrt(2 * gMAF .* (1 - gMAF)))';
      elseif kernel == "IBS1"
        # Table 3 in notes
        gSNP = hcat(geno[keepIdx, 1 : grpSize],
                    2 - geno[keepIdx, 1 : grpSize]);
        scale!(gSNP, [snpSqrtWts; snpSqrtWts]);
      elseif kernel == "IBS2"
        # Table 2 in notes
        gSNP = hcat(geno[keepIdx, 1 : grpSize] .>= 1,
                    geno[keepIdx, 1 : grpSize] .<= 1);
        scale!(gSNP, [snpSqrtWts; snpSqrtWts]);
      elseif kernel == "IBS3"
        # Table 1 in notes
        gSNP = hcat(-2 * (geno[keepIdx, 1 : grpSize] .<= 1) + 1,
                    2 * (geno[keepIdx, 1 : grpSize] .>= 1) - 1,
                    fill(sqrt(2), nPerKeep, grpSize));
        scale!(gSNP, [snpSqrtWts; snpSqrtWts; snpSqrtWts]);
      end

      # testing

      if kinship == "none"

        if test == "eRLRT"
          (b, vc0List, vc1List, testStats) =
            vctest(ynew, [], XtNullBasis * gSNP, WPreSim = WPreSim, tests = test,
                   Vform = "half", nBlockAscent = nBlockAscent, nMMmax = nMMmax,
                   nNullSimPts = nNullSimPts, pvalueComputings = pvalueComputing,
                   nNullSimNewtonIter = nNullSimNewtonIter, tolX = tolX,
                   devices = device);
        else
          (b, vc0List, vc1List, testStats) =
            vctest(y, X, gSNP, WPreSim = WPreSim, tests = test,
                   Vform = "half", nBlockAscent = nBlockAscent, nMMmax = nMMmax,
                   nNullSimPts = nNullSimPts, pvalueComputings = pvalueComputing,
                   nNullSimNewtonIter = nNullSimNewtonIter, tolX = tolX,
                   devices = device);
        end
        pvalList = testStats.vc1_pvalue;

      else

        # obtain some orthonormal vectors of space R^n - [QX,QPHI,S]
        # See Golub and Van Loan (1996) Algorithm 12.4.2 on p602
        tmpMat = Array(Float64, grpSize, size(XPhitNullBasis, 2));
        BLAS.gemm!('T', 'N', 1.0, gSNP, XPhitNullBasis, 0.0, tmpMat);
        (UXPhitS, svalXPhitS, VXPhitS) = svd(tmpMat, thin = false);
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
        VWorkSqrt = Array(Float64, length(evalPhiAdj), size(gSNP, 2));
        VWorkSqrt2 = Array(Float64, length(evalPhiAdj), size(gSNP, 2));
        BLAS.gemm!('T', 'N', 1.0, QPhi, gSNP, 0.0, VWorkSqrt2);
        BLAS.gemm!('T', 'N', 1.0, weightedW, VWorkSqrt2, 0.0, VWorkSqrt);
        #scale!(InvSqrtevalPhiAdj, VWorkSqrt);

        # obtain LRT and p-value
        (b, vc0List, vc1List, testStats) =
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
                 denomvec = denomvec, d1f = d1f, d2f = d2f, offset = offset);
        pvalList = testStats.vc1_pvalue;

      end

      # display progress
      if nReads > 1
        if idxRead == percentVec[percentIdx]
          println(percentIdx, "0 %");
          percentIdx += 1;
        end
      end

      # write output
      fid = open(outFile, "a");
      println(fid, snpIDred[1], ",", snpIDred[counter], ",", pvalList);
      close(fid);

    end

    # summary
    println("There are ", nSNP, " SNPs in total.");
    println(nActualSNP, " SNPs were analysed.");

  end

end
