% function [i,j,k] = index_to_ijk(indexcount, N1, N2, N3)
%     k = min(floor(indexcount/(N1*N2)) + 1, N3);
%     j = floor(rem(indexcount, (N1*N2)) / N1) + 1;
%     i = rem(rem(indexcount, (N1*N2)), N1);
%     if rem(rem(indexcount, (N1*N2)), N1) == 0
%         i = N1;
%         j = j - 1;
%     end
%     if rem(indexcount, (N1*N2)) == 0
%         j = N2;
%         k = k - 1;
%     end
% end


function [i,j,k] = index_to_ijk(idx, N1, N2, N3)
    k = floor((idx-1)/(N1*N2)) + 1;
    rem1 = idx - (k-1)*N1*N2;
    j = floor((rem1-1)/N1) + 1;
    i = rem1 - (j-1)*N1;
end