function [Fx,Fy,Fz] = grad3nu(F,X,Y,Z)
Fx = g(F,X,1); Fy = g(F,Y,2); Fz = g(F,Z,3);

    function dF = g(F,c,dim)
        c=c(:); N=numel(c); sz=size(F);
        p=[dim,1:dim-1,dim+1:ndims(F)];
        A=permute(F,p); A=reshape(A,N,[]);
        dA=zeros(size(A),'like',A);

        dc=diff(c);
        dA(1,:)   = (A(2,:)-A(1,:))/dc(1);          % 1st-order forward
        dA(N,:)   = (A(N,:)-A(N-1,:))/dc(end);      % 1st-order backward
        den = (c(3:end)-c(1:end-2));                % nonuniform central
        dA(2:N-1,:) = (A(3:end,:)-A(1:end-2,:)) ./ den;

        dF = ipermute(reshape(dA,[N,sz(p(2:end))]),p);
    end
end
