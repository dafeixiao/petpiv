function [Sol3,count] = solver3_acc_2(point,X,Y,Z,ax,ay,az,UNCAx,UNCAy,UNCAz,R,minN,xprior,S,up,weight,visc,T,k,itr)
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
ax=ax(valid);
ay=ay(valid);
az=az(valid);
weight=weight(valid);
T=T(valid);
UNCAx=UNCAx(valid);
UNCAy=UNCAy(valid);
UNCAz=UNCAz(valid);

% valid=(sqrt(ax.^2+ay.^2+az.^2)<150);
% X=X(valid);
% Y=Y(valid);
% Z=Z(valid);
% ax=ax(valid);
% ay=ay(valid);
% az=az(valid);
% weight=weight(valid);
% T=T(valid);
% UNCAx=UNCAx(valid);
% UNCAy=UNCAy(valid);
% UNCAz=UNCAz(valid);

% valid=det_outlier(ax,ay,az,20,0.1);
% % valid=abs(ax-mean(ax))<10*abs(mean(ax));
valid=abs(ax-mean(ax))<=3*std(ax) & abs(ay-mean(ay))<=3*std(ay) & abs(az-mean(az))<=3*std(az);
X=X(valid);
Y=Y(valid);
Z=Z(valid);
ax=ax(valid);
ay=ay(valid);
az=az(valid);
weight=weight(valid);
T=T(valid);
UNCAx=UNCAx(valid);
UNCAy=UNCAy(valid);
UNCAz=UNCAz(valid);


count=length(X);
Sol3=zeros(15,1);
flag=0;
scale=1;
if count>30 %% change to 5 for non-erroneous data
    R=(Xg-X).^2./(R1^2) + (Yg-Y).^2./(R2^2)+ (Zg-Z).^2./(R3^2);
    [B,I] = sort(R);
    valid=I(1:30);
    X=X(valid);
    Y=Y(valid);
    Z=Z(valid);
    ax=ax(valid);
    ay=ay(valid);
    az=az(valid);
    weight=weight(valid);
    T=T(valid);
    UNCAx=UNCAx(valid);
    UNCAy=UNCAy(valid);
    UNCAz=UNCAz(valid);
end



count=length(X);
if  count>=minN
    y=zeros(count*3,1);
    A=zeros(count*3,15);
    w=zeros(count*3,1);
%     yes=abs(ax-mean(ax))<10*mean(ax);
    for i=1:count
        temp=[X(i)-Xg,Y(i)-Yg,Z(i)-Zg];
        y((i-1)*3+1:i*3)=[ax(i);ay(i);az(i)];
        A((i-1)*3+1,:)=[1,0,0,temp,0,0,0,0,0,0,T(i),0,0];
        A((i-1)*3+2,:)=[0,1,0,0,0,0,temp,0,0,0,0,T(i),0];
        A((i-1)*3+3,:)=[0,0,1,0,0,0,0,0,0,temp,0,0,T(i)];
        %         if abs(U(i)-mean(U))<0.30*mean(U)
        %          w((i-1)*3+1:i*3)=[1/0.1,1/0.1,1/0.1];
%         if nnz(yes)>1 && abs(ax(i)-mean(ax))<mean(ax)
%             w((i-1)*3+1:i*3)=[1/cov(ax(yes)),1/cov(ay(yes)),1/cov(az(yes))];
% %             w(isnan(w))=10^-5;
%         elseif count>1
            w((i-1)*3+1:i*3)=[1/cov(ax),1/cov(ay),1/cov(az)];
%         else
%             w((i-1)*3+1:i*3)=[1/2,1/2,1/2].*10^-6;
%         end
%         if weight(i)~=0
%             w((i-1)*3+1:i*3)=w((i-1)*3+1:i*3).*weight(i);
%         end
%         w(isinf(w))=10^-5;


%         w((i-1)*3+1:i*3)=[UNCAx(i),UNCAy(i),UNCAz(i)];

    end

%         xprior(4:12)=xprior(4:12)./scale;
%         %         xprior(13:30)=xprior(13:30)./scale^2;
%         s(4:12)=s(4:12).*scale;
%         visc=visc./scale^3;
%         w=w*scale;
%         %         s(13:30)=s(13:30).*scale^2;
%     end
%     
%      div=sqrt((xprior(11)-xprior(9))^2+(xprior(6)-xprior(10))^2+(xprior(7)-xprior(5))^2);
%     S=S/(div+1);
%     if nnz(xprior(:))<4
%         S=S./3;
%     end
    S=S.*k;
%     Sp1=zeros(15,15);
% Sp1(1:3,1:3)=[0 0 0;0 1 -1;0 1 -1].*norm(S(1:3,1:3))./100;
% Sp2=zeros(15,15);
% Sp2(1:3,1:3)=[1 0 -1;0 0 0;1 0 -1].*norm(S(1:3,1:3))./100;
% Sp3=zeros(15,15);
% Sp3(1:3,1:3)=[-1 1 0;-1 1 0;0 0 0].*norm(S(1:3,1:3))./100;
% up1=zeros(15,1);
% up1(1:3)=[0;up(1);up(2)];
% up2=zeros(15,1);
% up2(1:3)=[up(3);0;up(4)];
% up3=zeros(15,1);
% up3(1:3)=[up(5);up(6);0];
 Wi=diag(w);
%    w=[ax;ay;az];
%    if det(cov(w'))~=0
%        M=inv(cov(w'));
%        for j=1:count
%            Wi(3*j-2:3*j,3*j-2:3*j)=M;
%        end
%    end
    b=zeros(15,1);
    c=zeros(15,1);
    d=zeros(15,1);
    b(11)=1;
    b(9)=-1;
    c(6)=1;
    c(10)=-1;
    d(7)=1;
    d(5)=-1;
    B=[b,c,d];
    if nnz(up)==0
        I=pinv(S+transpose(A)*Wi*A);
        O(1,:) = [transpose(b)*(I*b), transpose(b)*(I*c), transpose(b)*(I*d)];
        O(2,:) = [transpose(c)*(I*b), transpose(c)*(I*c), transpose(c)*(I*d)];
        O(3,:) = [transpose(d)*(I*b), transpose(d)*(I*c), transpose(d)*(I*d)];
%         O=transpose(B)*I*B;
        
%         figure(1)
%         scatter(Y,ay);
%         hold on

%         figure(2)
%         hold on
%         plot([1:12],xprior(1:12),'-k');
%         hold on
%         xlim([0 13]);
%         xticks([1:12]);
%         xticklabels({'ax','ay','az','daxdx','daxdy','daxdz','daydx','daydy','daydz','dazdx','dazdy','dazdz'});
        
        for it=1:itr
            P=[transpose(b)*(I*(transpose(A)*Wi*y+S*xprior))-visc(1);transpose(c)*(I*(transpose(A)*Wi*y+S*xprior))-visc(2);transpose(d)*(I*(transpose(A)*Wi*y+S*xprior))-visc(3)];
%             P=transpose(B)*I*(transpose(A)*Wi*y+S*xprior)-transpose(visc);
            Mu = linsolve(O,P);
%             xest=I*(transpose(A)*Wi*y+S*xprior-B*Mu);
            xest=I*(transpose(A)*Wi*y+S*xprior-Mu(1)*b-Mu(2)*c-Mu(3)*d);

% J=transpose(y-A*xest)*Wi*(y-A*xest)+transpose(xest-xprior)*(S/k)*(xest-xprior)+transpose(Mu)*(transpose(B)*xest-transpose(visc));       
            
            xprior=xest;            

%             figure(3)
%             if mod(it,2)==1
%             plot(it,J,'ok','markersize',5);
%             else
%             plot(it,J,'^k','markersize',5);
%             end
%             hold on
%             xlabel('iteration');
%             ylabel('J');
%             figure(2)
%             plot([1:12],xprior(1:12));
%             hold on
%             xlim([0 13]);
%             xticks([1:12]);
%             xticklabels({'ax','ay','az','daxdx','daxdy','daxdz','daydx','daydy','daydz','dazdx','dazdy','dazdz'});

        end
    else
        I=pinv(S+Sp1+Sp2+Sp3+transpose(A)*Wi*A);
        O(1,:) = [transpose(b)*(I*b), transpose(b)*(I*c), transpose(b)*(I*d)];
        O(2,:) = [transpose(c)*(I*b), transpose(c)*(I*c), transpose(c)*(I*d)];
        O(3,:) = [transpose(d)*(I*b), transpose(d)*(I*c), transpose(d)*(I*d)];
        for i=1:2
            P=[transpose(b)*(I*(transpose(A)*Wi*y+S*xprior-Sp1*up1-Sp2*up2-Sp3*up3))-visc(1);transpose(c)*(I*(transpose(A)*Wi*y+S*xprior-Sp1*up1-Sp2*up2-Sp3*up3))-visc(2);transpose(d)*(I*(transpose(A)*Wi*y+S*xprior-Sp1*up1-Sp2*up2-Sp3*up3))-visc(3)];
            Mu = linsolve(O,P);
            xprior=I*(transpose(A)*Wi*y+S*xprior-Sp1*up1-Sp2*up2-Sp3*up3-Mu(1)*b-Mu(2)*c-Mu(3)*d);
        end
    end
    Sol3=xest;
    if flag==1
        Sol3(4:12)=Sol3(4:12).*scale;
        %     Sol3(13:30)=Sol3(13:30).*scale^2;
    end
%     count=[cov(ax),cov(ay),cov(az)].^0.5;
end