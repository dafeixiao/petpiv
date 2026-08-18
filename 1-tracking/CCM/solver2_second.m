function [Sol12,count] = solver2_second(point,X,Y,Z,U,V,W,R,minN,T, UNCUp, UNCVp, UNCWp,cutoff,zloc)

Xg = point(1);
Yg = point(2);
Zg = point(3);
R1 = R(1);
R2 = R(2);
R3 = R(3);
valid=((Xg-X).^2./(R1^2) + (Yg-Y).^2./(R2^2)+ (Zg-Z).^2./(R3^2)<1);
X=X(valid);
Y=Y(valid);
Z=Z(valid);
U=U(valid);
V=V(valid);
W=W(valid);
T=T(valid);

% dafei 
UNCU = UNCUp(valid);
UNCV = UNCVp(valid);
UNCW = UNCWp(valid);
%


count=length(X);
if count>40
    R=(Xg-X).^2./(R1^2) + (Yg-Y).^2./(R2^2)+ (Zg-Z).^2./(R3^2);
    [B,I] = sort(R);
    valid=I(1:40);
    X=X(valid);
    Y=Y(valid);
    Z=Z(valid);
    U=U(valid);
    V=V(valid);
    W=W(valid);
    T=T(valid);
    UNCU = UNCUp(valid);
    UNCV = UNCVp(valid);
    UNCW = UNCWp(valid);
end

UNCU(UNCU==0) = 2*max(UNCU);
UNCV(UNCV==0) = 2*max(UNCV);
UNCW(UNCW==0) = 2*max(UNCW);


count=length(X);
Sol12=0;
flag=0;
clear R B
if  count>=minN
    MAT_LHS=zeros(count*3,1);
    A=zeros(count*3,33);
    for i=1:count
        dx=X(i)-Xg;
        dy=Y(i)-Yg;
        dz=Z(i)-Zg;
        temp=[dx,dy,dz];
        temp2=[0.5*dx^2,0.5*dy^2,0.5*dz^2,dx*dy,dx*dz,dy*dz];
        MAT_LHS((i-1)*3+1:i*3)=[U(i);V(i);W(i)];
        A((i-1)*3+1,:)=[1,0,0,temp,0,0,0,0,0,0,temp2,0,0,0,0,0,0,0,0,0,0,0,0,T(i),0,0];
        A((i-1)*3+2,:)=[0,1,0,0,0,0,temp,0,0,0,0,0,0,0,0,0,temp2,0,0,0,0,0,0,0,T(i),0];
        A((i-1)*3+3,:)=[0,0,1,0,0,0,0,0,0,temp,0,0,0,0,0,0,0,0,0,0,0,0,temp2,0,0,T(i)];

        w((i-1)*3+1) = 1 / UNCU(i)^2;   % U
        w((i-1)*3+2) = 1 / UNCV(i)^2;   % V
        w((i-1)*3+3) = 1 / UNCW(i)^2;   % W
    end

    %% Normal pinv, taking all singular values
    % Sol12 = pinv(A)*MAT_LHS;


    %% Truncated SVD
    [~,S,~] = svd(A,'econ');
    sv = diag(S);

    % % Take the knee, i.e. where singular values rapidly fall off
    % % [~,knee] = max(sv(1:end-1)./sv(2:end));
    % % Sol12 = pinv(A, sv(knee+1)) * MAT_LHS;
    % 
    % Take 1% of the largest singular value as the tolerance. Basically
    % says that any direction where the contribution is 100 times lower, is
    % probably not good, i.e., that could be noise
    [~,ind] = min(abs(zloc - Zg));
    Sol12 = pinv(A, cutoff(ind)*sv(1)) * MAT_LHS;

    % Energy based tolerance, which says that take at most 99% of the
    % energy containing in the singular values. Something like above. Can
    % change the 99% to something else. Modeled after paper-
    % A Rank Revealing Randomized Singular Value Decomposition (R3SVD)
    % Algorithm for Low- rank Matrix Approximations Hao Ji, Student Member, IEEE, Wenjian Yu, Member, IEEE, and Yaohang Li, Member, IEEE
    % k = find(cumsum(sv.^2)/sum(sv.^2) >= 0.999, 1);
    % Sol12 = pinv(A, sv(k)) * MAT_LHS;
end
