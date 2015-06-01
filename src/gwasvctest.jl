@doc """
gwasvctest(plinkFile::String = "", covFile::String = "", traitFile::String = "",
annotationFile::String = "", outFile::String = "", test::String = "eRLRT",
kinship::String = "GRM", kernel::String = "GRM", pvalueComputing::String = "chi2",
windowSize::Int = 50, infLambda::Float64 = 0.0, MemoryLimit::Int = 200000000,
nNullSimPts::Int = 10000, nMMmax::Int = 0, nBlockAscent::Int = 1000,
nNullSimNewtonIter::Int = 15, tolX::Float64 = 1e-4, snpWtType::String = "")

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

"tolX" - Tolerance in change of parameters, default is 1e-4.

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
                    test::String = "eRLRT", tolX::Float64 = 1e-4,
                    traitFile::String = "", MemoryLimit::Int = 200000000,
                    pvalueComputing::String = "chi2", windowSize::Int = 50,
                    annotationFile::String = "")

  ## parse data from files
  plinkBedfile = string(plinkFile, ".bed");
  plinkBimfile = string(plinkFile, ".bim");
  plinkFamfile = string(plinkFile, ".fam");

  # BIM file: chr, rs#, morgan, bp position, allele1, allele2
  bimdata = readdlm(plinkBimfile, String);
  chrID = bimdata[:, 1];
  snpID = bimdata[:, 2];
  nSNP = length(snpID);
  snpPos = bimdata[:, 4];

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
    geneName = Array(String, nSNP);
    readAnnotate!(annotationFile, snpID, nSNP, grpInfo, offsetSize, geneName);
    grpInfo = grpInfo[grpInfo .> 0];
    nGrp = length(grpInfo);
    offsetSize = offsetSize[1 : nGrp];
    geneName = geneName[1 : nGrp];
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
  ## missing covariate values are imputed by avarages

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
  if kinship != "none"
    kinMat = kinMat[keepIdx, keepIdx];
  end


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
      geneName = vcat(geneName[1:tmpidx-1], repmat([geneName[tmpidx]], innerItr, 1),
                      geneName[tmpidx+1:end]);
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
    Xsvd = svdfact(X, thin = false);
    rankX = countnz(Xsvd[:S] .> nPerKeep * eps(Xsvd[:S][1]));
    XtNullBasis = Xsvd[:U][:, rankX + 1 : end]';
    ynew = zeros(nPerKeep - rankX);
    BLAS.gemv!('N', 1.0, XtNullBasis, y, 0.0, ynew);

    rankQPhi = 0;
    evalPhiAdj = Float64[];
    XPhitNullBasis = [Float64[] Float64[]];
    yShift = Float64[];
    KPhiAdj = [Float64[] Float64[]];
    weightedW = [Float64[] Float64[]];
    QPhi = [Float64[] Float64[]];
    PrePartialSumW = [Float64[] Float64[]];
    PreTotalSumW = [Float64[] Float64[]];
    nPreRank = min(nPerKeep, windowSize);
    tmpn = nPerKeep - rankX;

  else
    # >=2 non-trivial variance components

    ynew = Float64[];

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

    # prepare for simulating Chi Squares and sums of Chi Squares
    nPreRank = min(length(evalPhiAdj), windowSize);
    tmpn = length(evalPhiAdj);

  end

  # simulate Chi Squares and sums of Chi Squares
  PrePartialSumW = Array(Float64, nNullSimPts, nPreRank+1);
  PreTotalSumW = Array(Float64, nNullSimPts, nPreRank+1);
  tmpSumVec = Array(Float64, nNullSimPts);
  p1 = pointer(PrePartialSumW);
  BLAS.blascopy!(nNullSimPts, rand(Chisq(tmpn), nNullSimPts), 1, p1, 1);
  for i = 1 : nNullSimPts
    #pW = pointer(WPreSim) + (i - 1) * windowSize * sizeof(Float64);
    #tmpSumVec[i] = BLAS.asum(1, pW, 1);
    tmpSumVec[i] = 0.0;
    PreTotalSumW[i, 1] = tmpSumVec[i] + PrePartialSumW[i, 1];
  end
  for j = 1 : nPreRank
    p1 = pointer(PrePartialSumW) + j * nNullSimPts * sizeof(Float64);
    BLAS.blascopy!(nNullSimPts, rand(Chisq(tmpn - j), nNullSimPts), 1, p1, 1);
    for i = 1 : nNullSimPts
      tmpSumVec[i] += WPreSim[j, i];
      PreTotalSumW[i, j+1] = tmpSumVec[i] + PrePartialSumW[i, j+1];
    end
  end

  # prepare output
  if isempty(outFile)
    outFile = string(plinkFile, "-julia.out");
  end
  fid = open(outFile, "w");
  if !flagAnnotate
    println(fid, "StartSNP,EndSNP,pvalue");
  else
    println(fid, "StartSNP,EndSNP,Chr,StartPos,EndPos,GeneName,#SNPs,pvalue");
  end
  close(fid);


  #### test group by group

  ## no annotation file
  if !flagAnnotate

    # preparation for parallel
    if nprocs() == 1
      ncores = 1;
    else
      ncores = nprocs() - 1;
    end
    nSNPParallel = Array(Int, ncores, 1);
    nReadsParallel = Array(Int, ncores, 1);
    rawdataParallel = Array(Matrix{Int8}, ncores, 1);
    snpIDParallel = Array(Vector{String}, ncores, 1);
    fill!(nSNPParallel, int(nSNP / ncores));
    nSNPParallel[ncores] = nSNP - sum(nSNPParallel[1 : ncores-1]);
    totaloffset = 0;
    for i = 1 : ncores
      if i > 1
        totaloffset += nSNPParallel[i-1];
      end
      nReadsParallel[i] = iceil(nSNPParallel[i] / SNPSize);
      rawdataParallel[i] =
        rawdata[:, totaloffset+1 : totaloffset+nSNPParallel[i]];
      snpIDParallel[i] = snpID[totaloffset+1 : totaloffset+nSNPParallel[i]];
    end
    # loop over windows
    outresults = pmap(loopwinFixsize, nSNPParallel, nReadsParallel,
                      rawdataParallel, snpIDParallel, {bin2geno for i = 1:ncores},
                      {nPer for i = 1:ncores}, {SNPSize for i = 1:ncores},
                      {flagAnnotate for i = 1:ncores},
                      {snpWtType for i = 1:ncores}, {keepIdx for i = 1:ncores},
                      {nPerKeep for i = 1:ncores}, {kernel for i = 1:ncores},
                      {kinship for i = 1:ncores}, {test for i = 1:ncores},
                      {XtNullBasis for i = 1:ncores}, {WPreSim for i = 1:ncores},
                      {nBlockAscent for i = 1:ncores}, {nMMmax for i = 1:ncores},
                      {nNullSimPts for i = 1:ncores}, {tolX for i = 1:ncores},
                      {pvalueComputing for i = 1:ncores},
                      {device for i = 1:ncores},
                      {nNullSimNewtonIter for i = 1:ncores},
                      {rankQPhi for i = 1:ncores}, {evalPhiAdj for i = 1:ncores},
                      {XPhitNullBasis for i = 1:ncores},
                      {yShift for i = 1:ncores}, {y for i = 1:ncores},
                      {KPhiAdj for i = 1:ncores}, {weightedW for i = 1:ncores},
                      {QPhi for i = 1:ncores}, {PrePartialSumW for i = 1:ncores},
                      {PreTotalSumW for i = 1:ncores},
                      {windowSize for i = 1:ncores},
                      {nPreRank for i = 1:ncores}, {X for i = 1:ncores},
                      {ynew for i = 1:ncores});

    # write output
    fid = open(outFile, "a");
    for i = 1 : ncores
      for j = 1 : size(outresults[i], 1)
        println(fid, outresults[i][j, 1], ",", outresults[i][j, 2], ",",
                outresults[i][j, 3]);
      end
    end
    close(fid);

  ## provide annotation file
  else

    # preparation for parallel
    if nprocs() == 1
      ncores = 1;
    else
      ncores = nprocs() - 1;
    end
    nReadsParallel = Array(Int, ncores, 1);
    rawdataParallel = Array(Matrix{Int8}, ncores, 1);
    grpInfoParallel = Array(Vector{Int}, ncores, 1);
    offsetSizeParallel = Array(Vector{Int}, ncores, 1);
    snpIDParallel = Array(Vector{String}, ncores, 1);
    chrIDParallel = Array(Vector{String}, ncores, 1);
    snpPosParallel = Array(Vector{String}, ncores, 1);
    geneNameParallel = Array(Vector{String}, ncores, 1);
    fill!(nReadsParallel, int(nReads / ncores));
    nReadsParallel[ncores] = nReads - sum(nReadsParallel[1 : ncores-1]);
    totaloffset = 0;
    for i = 1 : ncores
      if i == 1
        tmpoffset = 0;
      else
        tmpoffset = sum(nReadsParallel[1:i-1]);
        totaloffset += sum(grpInfoParallel[i-1]);
      end
      grpInfoParallel[i] = grpInfo[tmpoffset+1 : tmpoffset+nReadsParallel[i]];
      geneNameParallel[i] = geneName[tmpoffset+1 : tmpoffset+nReadsParallel[i]];
      offsetSizeParallel[i] =
        offsetSize[tmpoffset+1 : tmpoffset+nReadsParallel[i]] -
        offsetSize[tmpoffset+1];
      rawdataParallel[i] =
        rawdata[:, totaloffset+1 : totaloffset+sum(grpInfoParallel[i])];
      snpIDParallel[i] =
        snpID[totaloffset+1 : totaloffset+sum(grpInfoParallel[i])];
      chrIDParallel[i] =
        chrID[totaloffset+1 : totaloffset+sum(grpInfoParallel[i])];
      snpPosParallel[i] =
        snpPos[totaloffset+1 : totaloffset+sum(grpInfoParallel[i])];
    end

    # loop over genes
    outresults = pmap(loopwinAnnot, nReadsParallel, rawdataParallel,
                      grpInfoParallel, offsetSizeParallel, snpIDParallel,
                      {bin2geno for i = 1:ncores}, {nPer for i = 1:ncores},
                      {SNPSize for i = 1:ncores}, {flagAnnotate for i = 1:ncores},
                      {snpWtType for i = 1:ncores}, {keepIdx for i = 1:ncores},
                      {nPerKeep for i = 1:ncores}, {kernel for i = 1:ncores},
                      {kinship for i = 1:ncores}, {test for i = 1:ncores},
                      {XtNullBasis for i = 1:ncores}, {WPreSim for i = 1:ncores},
                      {nBlockAscent for i = 1:ncores}, {nMMmax for i = 1:ncores},
                      {nNullSimPts for i = 1:ncores}, {tolX for i = 1:ncores},
                      {pvalueComputing for i = 1:ncores},
                      {device for i = 1:ncores},
                      {nNullSimNewtonIter for i = 1:ncores},
                      {rankQPhi for i = 1:ncores}, {evalPhiAdj for i = 1:ncores},
                      {XPhitNullBasis for i = 1:ncores},
                      {yShift for i = 1:ncores}, {y for i = 1:ncores},
                      {KPhiAdj for i = 1:ncores}, {weightedW for i = 1:ncores},
                      {QPhi for i = 1:ncores}, {PrePartialSumW for i = 1:ncores},
                      {PreTotalSumW for i = 1:ncores},
                      {windowSize for i = 1:ncores},
                      {nPreRank for i = 1:ncores}, {X for i = 1:ncores},
                      snpPosParallel, chrIDParallel, geneNameParallel,
                      {ynew for i = 1:ncores});

    # write output
    fid = open(outFile, "a");
    for i = 1 : ncores
      for j = 1 : size(outresults[i], 1)
        println(fid, outresults[i][j, 1], ",", outresults[i][j, 2], ",",
                outresults[i][j, 3], ",", outresults[i][j, 4], ",",
                outresults[i][j, 5], ",", outresults[i][j, 6], ",",
                outresults[i][j, 7], ",", outresults[i][j, 8]);
      end
    end
    close(fid);

  end

end
