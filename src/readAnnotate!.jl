function readAnnotate!(annotationFile::String, snpID::Vector{String},
                       nSNP::Int, grpInfo::Vector{Int64},
                       offsetSize::Vector{Int64})

  annotatedata = readdlm(annotationFile, ',');
  annotatedata = convert(Array{String, 2}, annotatedata);
  annotatesnp = annotatedata[:, 2];
  nAnno = size(annotatedata, 1);
  nGrp = 0;
  annoIdx = 1;
  pregene = "";
  curgene = "";

  for i = 1 : nSNP

    while annoIdx <= nAnno
      if annotatesnp[annoIdx] == snpID[i]
        curgene = annotatedata[annoIdx, 1];
        annoIdx += 1;
        break;
      else
        annoIdx += 1;
      end
    end

    if curgene == pregene
      grpInfo[nGrp] += 1;
    else
      nGrp += 1;
      grpInfo[nGrp] += 1;
      if nGrp > 1
        offsetSize[nGrp] = offsetSize[nGrp - 1] + grpInfo[nGrp - 1];
      end
      pregene = curgene;
    end

  end

end
