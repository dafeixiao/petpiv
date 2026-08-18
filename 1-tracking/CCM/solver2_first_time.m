%% second order 4D SVD %%just 1st order
function [Sol12,count] = solver2_first_time(point,X,Y,Z,U,V,W,R,minN,T)
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

% % valid=det_outlier(U,V,W,20,0.1);
% % valid=(sqrt(U.^2+V.^2+W.^2)<600);
% valid=abs(U-mean(U))<=3*std(U) & abs(V-mean(V))<=3*std(V) & abs(W-mean(W))<=3*std(W);
% X=X(valid);
% Y=Y(valid);
% Z=Z(valid);
% U=U(valid);
% V=V(valid);
% W=W(valid);
% T=T(valid);

count=length(X);
MAT_LHS1=zeros(count,1);MAT_LHS2=zeros(count,1);MAT_LHS3=zeros(count,1);
A1=zeros(count,5);A2=zeros(count,5);A3=zeros(count,5);
for i=1:count
    temp=[X(i)-Xg,Y(i)-Yg,Z(i)-Zg];
    MAT_LHS1(i)=U(i);MAT_LHS2(i)=V(i);MAT_LHS3(i)=W(i);
    A1(i,:)=[1,temp,T(i)];A2(i,:)=[1,temp,T(i)];A3(i,:)=[1,temp,T(i)];
end
Sol1 = pinv(A1)*MAT_LHS1;Sol2 = pinv(A2)*MAT_LHS2;Sol3 = pinv(A3)*MAT_LHS3;
Res1=MAT_LHS1-A1*Sol1;Res2=MAT_LHS2-A2*Sol2;Res3=MAT_LHS3-A3*Sol3;
valid1=abs(Res1-mean(Res1))<=3*std(Res1);valid2=abs(Res2-mean(Res2))<=3*std(Res2);valid3=abs(Res3-mean(Res3))<=3*std(Res3);
valid=valid1 & valid2 & valid3;
X=X(valid);
Y=Y(valid);
Z=Z(valid);
U=U(valid);
V=V(valid);
W=W(valid);
T=T(valid);

count=length(X);
if count>25 %%max number of particles
    R=(Xg-X).^2./(R1^2) + (Yg-Y).^2./(R2^2)+ (Zg-Z).^2./(R3^2);
    [B,I] = sort(R);
    valid=I(1:25);
    X=X(valid);
    Y=Y(valid);
    Z=Z(valid);
    U=U(valid);
    V=V(valid);
    W=W(valid);
    T=T(valid);
end
count=length(X);
Sol12=zeros(15,1);
flag=0;
clear R

if  count>=minN
    MAT_LHS=zeros(count*3,1);
    A=zeros(count*3,15);
    for i=1:count
        temp=[X(i)-Xg,Y(i)-Yg,Z(i)-Zg];
        MAT_LHS((i-1)*3+1:i*3)=[U(i);V(i);W(i)];
        A((i-1)*3+1,:)=[1,0,0,temp,0,0,0,0,0,0,T(i),0,0];
        A((i-1)*3+2,:)=[0,1,0,0,0,0,temp,0,0,0,0,T(i),0];
        A((i-1)*3+3,:)=[0,0,1,0,0,0,0,0,0,temp,0,0,T(i)];
    end
    Sol12 = pinv(A)*MAT_LHS;
    
%     Res=(MAT_LHS-A*Sol12).^2;
%     Res_2=zeros(length(Res),1);
%     for i=1:length(Res_2)/3
%        Res_2((i-1)*3+1:i*3)=sum(Res((i-1)*3+1:i*3));
%     end
%     wi=var(Res_2)./Res_2;
%     Wi=diag(wi);
%     Sol12=inv(A'*Wi*A)*A'*Wi*MAT_LHS;

%     f=[zeros(15,1);ones(count*3,1);ones(count*3,1)];
%     A1=[];
%     b1=[];
%     Aeq=[A,eye(count*3),-eye(count*3)];
%     beq=MAT_LHS;
%     l=[-1000000*ones(15,1);zeros(count*3,1);zeros(count*3,1)];
%     u=[];
%     options=optimoptions('linprog','Display','none');
%     theta=linprog(f,A1,b1,Aeq,beq,l,u,options);
%     Sol12=theta(1:15);


end
end
% end