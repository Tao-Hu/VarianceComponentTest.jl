function kinshipcoef(father, mother, start, finish)
#
# Author: Kenneth Lange
#
# This function computes the kinship matrix of a pedigree.
# Children must appear after parents. The pedigree occurs
# in a block and is defined implicitly by giving fathers
# and mothers. Person i is a pedigree founder if and only if
# father[i] = 0.
#
   ped_size = finish - start + 1
   kinship = zeros(ped_size, ped_size)
   q = start - 1
   for i = 1 : ped_size
      if father[i + q] == 0
         kinship[i, i] = 0.5
      else
         j = father[i + q] - q
         k = mother[i + q] - q
         kinship[i, i] = 0.5 + 0.5 * kinship[j, k]
         for m = 1 : i - 1
            kinship[i, m] = 0.5 * (kinship[j, m] + kinship[k, m])
            kinship[m, i] = kinship[i, m]
         end
      end
   end
   return kinship
end
