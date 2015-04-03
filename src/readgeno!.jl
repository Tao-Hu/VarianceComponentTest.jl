function readgeno!(geno::Matrix{Float64}, curSNPSize::Int, nPer::Int,
                   SNPSize::Int, idxRead::Int, bin2geno::Matrix{Float64},
                   rawdata::Matrix{Int8}, offset::Int)

  offset = (idxRead - 1) * SNPSize;
  nrow = iceil(nPer / 4);

  for j = 1 : curSNPSize
    cumPer = 0;

    for i = 1 : nrow - 1
      curvalue = rawdata[i, offset + j];
      if curvalue >= 0
        curvalue += 1;
      else
        curvalue += 257;
      end
      for k = 1 : 4
        cumPer = cumPer + 1;
        geno[cumPer, j] = bin2geno[k, curvalue];
      end
    end

    curvalue = rawdata[nrow, offset + j];
    if curvalue >= 0
      curvalue += 1;
    else
      curvalue += 257;
    end
    for k = 1 : 4
      if cumPer < nPer
        cumPer = cumPer + 1;
        geno[cumPer, j] = bin2geno[k, curvalue];
      else
        break;
      end
    end

  end

end
