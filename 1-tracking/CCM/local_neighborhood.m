function [is, ie, js, je, ks, ke] = local_neighborhood(i,j,k,N1,N2,N3,halfWidth)
    is = max(1, i-halfWidth);
    ie = min(i+halfWidth, N1);
    js = max(1, j-halfWidth);
    je = min(j+halfWidth, N2);
    ks = max(1, k-halfWidth);
    ke = min(k+halfWidth, N3);
end
