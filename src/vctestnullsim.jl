function vctestnullsim(teststat, evalV, evalAdjV, n, rankX, WPreSim;
                       device::String = "CPU", nSimPts::Int = 10000,
                       nNewtonIter::Int = 15, test::String = "eLRT",
                       pvalueComputing::String = "chi2",
                       PrePartialSumW::Array{Float64, 2} = [Float64[] Float64[]],
                       PreTotalSumW::Array{Float64, 2} = [Float64[] Float64[]],
                       partialSumWConst::Array{Float64, 1} = Float64[],
                       totalSumWConst::Array{Float64, 1} = Float64[],
                       windowSize::Int = 50,
                       partialSumW::Array{Float64, 1} = Float64[],
                       totalSumW::Array{Float64, 1} = Float64[],
                       lambda::Array{Float64, 2} = [Float64[] Float64[]],
                       W::Array{Float64, 2} = [Float64[] Float64[]],
                       nPreRank::Int = 20,
                       tmpmat0::Array{Float64, 2} = [Float64[] Float64[]],
                       tmpmat1::Array{Float64, 2} = [Float64[] Float64[]],
                       tmpmat2::Array{Float64, 2} = [Float64[] Float64[]],
                       tmpmat3::Array{Float64, 2} = [Float64[] Float64[]],
                       tmpmat4::Array{Float64, 2} = [Float64[] Float64[]],
                       tmpmat5::Array{Float64, 2} = [Float64[] Float64[]],
                       denomvec::Array{Float64, 1} = Float64[],
                       d1f::Array{Float64, 1} = Float64[],
                       d2f::Array{Float64, 1} = Float64[], offset::Int = 0)
  # VCTESTNULLSIM Simulate null distribution for testing zero var. component
  #
  # VCTESTNULLSIM(evalV,evalAdjV,n,rankX,WPreSim) simulate the null distributions
  # of various tests for testing $H_0:simga_1^2=0$ in the variance components
  # model $Y \sim N(X \beta, \sigma_0^2 I + \sigma_1^2 V)$.
  #
  #    OUTPUT:
  #        pvalue: p-value of the given test

  # test statistics = 0 means p-value=1
  if teststat <= 0
    return pvalue = 1;
  end

  # obtain rank information
  rankAdjV = length(evalAdjV);
  rankV = length(evalV);

  # check dimensions
  if n < rankX
    error("vctestnullsim:nLessThanRankX\n", "n should be larger than rankX");
  end
  if rankAdjV > n - rankX
    error("vctestnullsim:wrongRankAdjV\n", "rankAdjv should be <= n-rankX");
  end

  #W = Float64[];

  if device == "CPU"

    # Preparations
    nPtsChi2 = 300;
    if size(WPreSim, 1) < rankAdjV
      newSim = randn(rankAdjV - size(WPreSim, 1), nSimPts);
      if isempty(WPreSim)
        WPreSim = newSim;
      else
        WPreSim = [WPreSim, newSim];
      end
    end

    if isempty(PrePartialSumW) || isempty(PreTotalSumW)
      partialSumWConst = rand(Chisq(n - rankX - rankAdjV), nSimPts);
      totalSumWConst = similar(partialSumWConst);
      for i = 1 : nSimPts
        pW = pointer(WPreSim) + (i - 1) * windowSize * sizeof(Float64);
        totalSumWConst[i] = BLAS.asum(rankAdjV, pW, 1) + partialSumWConst[i];
      end
    else
      if rankAdjV >= nPreRank
        partialSumWConst = rand(Chisq(n - rankX - rankAdjV), nSimPts);
        for i = 1 : nSimPts
          pW = pointer(WPreSim) + (i - 1) * windowSize * sizeof(Float64);
          totalSumWConst[i] = BLAS.asum(rankAdjV, pW, 1) + partialSumWConst[i];
        end
      else
        ppw = pointer(partialSumWConst);
        ptw = pointer(totalSumWConst);
        pppw = pointer(PrePartialSumW) +
          (rankAdjV - 1) * nSimPts * sizeof(Float64);
        pptw = pointer(PreTotalSumW) +
          (rankAdjV - 1) * nSimPts * sizeof(Float64);
        BLAS.blascopy!(nSimPts, pppw, 1, ppw, 1);
        BLAS.blascopy!(nSimPts, pptw, 1, ptw, 1);
      end
    end

    # create null distribution samples
    if test == "eLRT" || test == "eRLRT"

      # set effective sample size
      if test == "eLRT"
        ne = n;
      else
        ne = n - rankX;
      end

      # find number of non-zero points
      if test == "eLRT"
        #nzIdx = sum(W .* evalAdjV, 1) ./ totalSumWConst' .> sum(evalV) / ne;
        threshold = BLAS.asum(length(evalV), evalV, 1) / ne;
      else
        #nzIdx = sum(W .* evalAdjV, 1) ./ totalSumWConst' .> sum(evalAdjV) / ne;
        threshold = BLAS.asum(length(evalAdjV), evalAdjV, 1) / ne;
      end

      if pvalueComputing == "MonteCarlo"
        counter = 0;
        for i = 1 : nSimPts
          pW = pointer(WPreSim) + (i - 1) * windowSize * sizeof(Float64);
          if BLAS.dot(rankAdjV, pW, 1, evalAdjV, 1) / totalSumWConst[i] > threshold
            counter += 1;
            partialSumW[counter] = partialSumWConst[i];
            totalSumW[counter] = totalSumWConst[i];
            pWupdated = pointer(W) + (counter - 1) * windowSize * sizeof(Float64);
            BLAS.blascopy!(rankAdjV, pW, 1, pWupdated, 1);
          end
        end
        #lambda = zeros(1, counter);
        #fill!(lambda, 0.0);
        fill!(lambda, 1e-5);
      else
        counter = 0;
        for i = 1 : nSimPts
          pW = pointer(WPreSim) + (i - 1) * windowSize * sizeof(Float64);
          if BLAS.dot(rankAdjV, pW, 1, evalAdjV, 1) / totalSumWConst[i] > threshold
            counter += 1;
            if counter <= nPtsChi2
              partialSumW[counter] = partialSumWConst[i];
              totalSumW[counter] = totalSumWConst[i];
              pWupdated = pointer(W) + (counter - 1) * windowSize * sizeof(Float64);
              BLAS.blascopy!(rankAdjV, pW, 1, pWupdated, 1);
            end
          end
        end
        patzero = (nSimPts - counter) / nSimPts;
        counter = nPtsChi2;
        #fill!(lambda, 0.0);
        fill!(lambda, 1e-5);
      end

      # Newton-Raphson iteration

      if test == "eRLRT"

        for iNewton = 1 : nNewtonIter
          #BLAS.gemm!('N', 'N', 1.0, reshape(evalAdjV, rankAdjV, 1), lambda,
          #           0.0, tmpmat0);

          for j = 1 : counter
            tmpsum1 = 0.0;
            tmpsum2 = 0.0;
            tmpsum3 = 0.0;
            tmpsum4 = 0.0;
            tmpsum5 = 0.0;
            for i = 1 : rankAdjV
              #tmpmat0[i, j] += 1;
              tmpmat0[i, j] = evalAdjV[i] * lambda[j] + 1;
              tmpmat1[i, j] = W[i, j] / tmpmat0[i, j];
              tmpsum1 += tmpmat1[i, j];
              #tmpmat1[i, j] = W[i, j] * tmpmat0[i, j];
              tmpmat4[i, j] =  evalAdjV[i] / tmpmat0[i, j];
              tmpsum4 += tmpmat4[i, j];
              #tmpmat2[i, j] = tmpmat1[i, j] / tmpmat0[i, j] * evalAdjV[i];
              tmpmat2[i, j] = tmpmat1[i, j] * tmpmat4[i, j];
              tmpsum2 += tmpmat2[i, j];
              #tmpmat3[i, j] = tmpmat2[i, j] / tmpmat0[i, j] * evalAdjV[i];
              tmpmat3[i, j] = tmpmat2[i, j] * tmpmat4[i, j];
              tmpsum3 += tmpmat3[i, j];
              #tmpmat4[i, j] = evalAdjV[i] / tmpmat0[i, j];
              #tmpmat4[i, j] = evalAdjV[i] * tmpmat0[i, j];
              tmpmat5[i, j] = tmpmat4[i, j] ^ 2;
              tmpsum5 += tmpmat5[i, j];
            end
            #p1 = pointer(tmpmat1) + (j - 1) * windowSize * sizeof(Float64);
            #denomvec[j] = BLAS.asum(rankAdjV, p1, 1) + partialSumW[j];
            denomvec[j] = tmpsum1 + partialSumW[j];
            #p2 = pointer(tmpmat2) + (j - 1) * windowSize * sizeof(Float64);
            #d1f[j] = BLAS.asum(rankAdjV, p2, 1) / denomvec[j];
            d1f[j] = tmpsum2 / denomvec[j];
            #p3 = pointer(tmpmat3) + (j - 1) * windowSize * sizeof(Float64);
            #d2f[j] = d1f[j] ^ 2 - 2 * BLAS.asum(rankAdjV, p3, 1) / denomvec[j];
            d2f[j] = d1f[j] ^ 2 - 2 * tmpsum3 / denomvec[j];
            #p4 = pointer(tmpmat4) + (j - 1) * windowSize * sizeof(Float64);
            #d1f[j] -= BLAS.asum(rankAdjV, p4, 1) / ne;
            d1f[j] -= tmpsum4 / ne;
            #p5 = pointer(tmpmat5) + (j - 1) * windowSize * sizeof(Float64);
            #d2f[j] += BLAS.asum(rankAdjV, p5, 1) / ne;
            d2f[j] += tmpsum5 / ne;
            if d2f[j] >= 0
              d2f[j] = -1;
            end
            lambda[j] -= d1f[j] / d2f[j];
            if lambda[j] < 0
              lambda[j] = 0;
            end
          end
        end

      else

        for iNewton = 1 : nNewtonIter
          #BLAS.gemm!('N', 'N', 1.0, reshape(evalAdjV, rankAdjV, 1), lambda,
          #           0.0, tmpmat0);

          for j = 1 : counter
            tmpsum1 = 0.0;
            tmpsum2 = 0.0;
            tmpsum3 = 0.0;
            tmpsum4 = 0.0;
            tmpsum5 = 0.0;
            for i = 1 : rankAdjV
              #tmpmat0[i, j] += 1;
              tmpmat0[i, j] = evalAdjV[i] * lambda[j] + 1;
              tmpmat1[i, j] = W[i, j] / tmpmat0[i, j];
              tmpsum1 += tmpmat1[i, j];
              tmpmat2[i, j] = tmpmat1[i, j] / tmpmat0[i, j] * evalAdjV[i];
              tmpsum2 += tmpmat2[i, j];
              tmpmat3[i, j] = tmpmat2[i, j] / tmpmat0[i, j] * evalAdjV[i];
              tmpsum3 += tmpmat3[i, j];
            end
            for i = 1 : rankV
              tmpmat4[i, j] = 1 / (1 / evalV[i] + lambda[j]);
              tmpsum4 += tmpmat4[i, j];
              tmpmat5[i, j] = tmpmat4[i, j] ^ 2;
              tmpsum5 += tmpmat5[i, j];
            end
            #p1 = pointer(tmpmat1) + (j - 1) * windowSize * sizeof(Float64);
            #denomvec[j] = BLAS.asum(rankAdjV, p1, 1) + partialSumW[j];
            denomvec[j] = tmpsum1 + partialSumW[j];
            #p2 = pointer(tmpmat2) + (j - 1) * windowSize * sizeof(Float64);
            #d1f[j] = BLAS.asum(rankAdjV, p2, 1) / denomvec[j];
            d1f[j] = tmpsum2 / denomvec[j];
            #p3 = pointer(tmpmat3) + (j - 1) * windowSize * sizeof(Float64);
            #d2f[j] = d1f[j] ^ 2 - 2 * BLAS.asum(rankAdjV, p3, 1) / denomvec[j];
            d2f[j] = d1f[j] ^ 2 - 2 * tmpsum3 / denomvec[j];
            #p4 = pointer(tmpmat4) + (j - 1) * windowSize * sizeof(Float64);
            #d1f[j] -= BLAS.asum(rankV, p4, 1) / ne;
            d1f[j] -= tmpsum4 / ne;
            #p5 = pointer(tmpmat5) + (j - 1) * windowSize * sizeof(Float64);
            #d2f[j] += BLAS.asum(rankV, p5, 1) / ne;
            d2f[j] += tmpsum5 / ne;
            if d2f[j] >= 0
              d2f[j] = -1;
            end
            lambda[j] -= d1f[j] / d2f[j];
            if lambda[j] < 0
              lambda[j] = 0;
            end
          end
        end

      end

      # collect null distribution samples
      if test == "eLRT"
        tmpmat6 = Array(Float64, rankV, counter);
        BLAS.gemm!('N', 'N', 1.0, reshape(evalV, rankV, 1), lambda, 0.0,
                   tmpmat6);
        simnull = Array(Float64, counter);
        for j = 1 : counter
          tmpprod = 1.0;
          for i = 1 : rankV
            tmpprod *= (tmpmat6[i, j] + 1);
          end
          simnull[j] = ne * log(totalSumW[j]) - ne * log(denomvec[j]) -
            log(tmpprod);
        end
      else
        simnull = Array(Float64, counter);
        for j = 1 : counter
          tmpprod = 1.0;
          for i = 1 : rankAdjV
            tmpprod *= tmpmat0[i, j];
          end
          simnull[j] = ne * log(totalSumW[j]) - ne * log(denomvec[j]) -
            log(tmpprod);
        end
      end
      if pvalueComputing == "MonteCarlo"
        if nSimPts > counter
          simnull = [simnull; zeros(nSimPts - counter)];
        end
      end
      simnull = simnull[simnull .>= 0];

    elseif test == "eScore"

      threshold = sum(evalV) / n;
      if pvalueComputing == "MonteCarlo"
        simnull = Array(Float64, nSimPts);
        for j = 1 : nSimPts
          tmpsum = 0.0;
          for i = 1 : rankAdjV
            tmpsum += WPreSim[i, j] * evalAdjV[i];
          end
          simnull[j] = max(tmpsum / totalSumWConst[j], threshold);
        end
      else
        simnull = Array(Float64, nPtsChi2);
        counter = 0;
        for j = 1 : nSimPts
          tmpsum = 0.0;
          for i = 1 : rankAdjV
            tmpsum += WPreSim[i, j] * evalAdjV[i];
          end
          tmpvalue = tmpsum / totalSumWConst[j];
          if tmpvalue > threshold
            counter += 1;
            if counter <= nPtsChi2
              simnull[counter] = tmpvalue;
            end
          end
        end
        patzero = (nSimPts - counter) / nSimPts;
      end

    end

    # compute empirical p-value
    if pvalueComputing == "MonteCarlo"
      pvalue = countnz(simnull .>= teststat) / length(simnull);
      return pvalue;
    elseif pvalueComputing == "chi2"
      ahat = var(simnull) / (2 * mean(simnull));
      bhat = 2 * (mean(simnull) ^ 2) / var(simnull);
      if mean(simnull) == 0 && var(simnull) == 0
        pvalue = 1.0;
        println("Invalid p-value occur! Starting from row ", offset+1);
      elseif length(simnull) == 0
        pvalue = 1.0;
        println("All simulated statistics are negative! Starting from row ", offset+1);
      else
        #println("mean of simnull = ", mean(simnull), ", var = ", var(simnull));
        #println("length of simnull = ", length(simnull), ", sum = ", sum(simnull));
        pvalue = (1 - patzero) * (1 - cdf(Chisq(bhat), teststat / ahat));
        if teststat == 0; pvalue = pvalue + patzero; end;
      end
      return pvalue;
    end

  elseif device == "GPU"

    error("vctestnullsim:notSupportedMode\n", "The GPU mode is currently not supported");

  end

end
