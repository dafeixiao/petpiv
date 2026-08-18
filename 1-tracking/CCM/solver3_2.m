function [Sol3,count] = solver3_2(point,X,Y,Z,U,V,W,UNCU,UNCV,UNCW,R,minN,xprior,S,up,weight,T,k,itr)
global delta
global scale
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
weight=weight(valid);
T=T(valid);
UNCU=UNCU(valid);
UNCV=UNCV(valid);
UNCW=UNCW(valid);
count=max(size(X));
Sol3=0;
% newR=R1;
scale=1;
flag=0;
count=length(X);
if count>30
    R=(Xg-X).^2./(R1^2) + (Yg-Y).^2./(R2^2)+ (Zg-Z).^2./(R3^2);
    [B,I] = sort(R);
    valid=I(1:30);
    X=X(valid);
    Y=Y(valid);
    Z=Z(valid);
    U=U(valid);
    V=V(valid);
    W=W(valid);
    T=T(valid);
    weight=weight(valid);
    UNCU=UNCU(valid);
    UNCV=UNCV(valid);
    UNCW=UNCW(valid);
    count=30;
end
 deltanu=1.0006e-3;

if  count>=minN &&nnz(T)>0
    y=zeros(count*3,1);
    A=zeros(count*3,33);
    for i=1:count
        dx=X(i)-Xg;
        dy=Y(i)-Yg;
        dz=Z(i)-Zg;
        temp=[dx,dy,dz];
        temp2=[0.5*dx^2,0.5*dy^2,0.5*dz^2,dx*dy,dx*dz,dy*dz];
        y((i-1)*3+1:i*3)=[U(i);V(i);W(i)];
                 A((i-1)*3+1,:)=[1,0,0,temp,0,0,0,0,0,0,temp2,0,0,0,0,0,0,0,0,0,0,0,0,T(i),0,0];
                A((i-1)*3+2,:)=[0,1,0,0,0,0,temp,0,0,0,0,0,0,0,0,0,temp2,0,0,0,0,0,0,0,T(i),0];
                A((i-1)*3+3,:)=[0,0,1,0,0,0,0,0,0,temp,0,0,0,0,0,0,0,0,0,0,0,0,temp2,0,0,T(i)];
 
        %                 if abs(U(i)-mean(U))<0.30*mean(U)
        %                  w((i-1)*3+1:i*3)=[1/0.1,1/0.1,1/0.1];
        %         if nnz(yes)>1 && abs(U(i)-mean(U))<mean(U)
        %             w((i-1)*3+1:i*3)=[1/cov(U(yes)),1/cov(V(yes)),1/cov(W(yes))];
        %         elseif count>1
%         w((i-1)*3+1:i*3)=[1/cov(U),1/cov(V),1/cov(W)];
%         %         else
%         %             w((i-1)*3+1:i*3)=[1/0.05,1/0.01,1/0.02].^2;
%         %         end
%                 if weight(i)~=0
%                     w((i-1)*3+1:i*3)=w((i-1)*3+1:i*3).*weight(i);
%                 end
      
          w((i-1)*3+1:i*3)=[UNCU(i),UNCV(i),UNCW(i)];
                
    end
%         xprior(4:12)=xprior(4:12)./scale;
%                 xprior(13:30)=xprior(13:30)./scale^2;
%         s(4:12)=s(4:12).*scale;
%         w=w.*scale;
%         %         s(13:30)=s(13:30).*scale^2;
%     end
   Wi=diag(w);
%    w=[U;V;W];
%    if det(cov(w'))~=0
%        M=inv(cov(w'));
%        for j=1:count
% Wi(3*j-2:3*j,3*j-2:3*j)=M;
% end
% end


 
%     div=((xprior(4)+xprior(8)+xprior(12))^2);
%     if div~=0
%         div=div/(xprior(4)^2+xprior(8)^2+xprior(12)^2);
%     end
%     S=S./(div+1);
%     if nnz(xprior(:))<4
%         S=S./3;
%     end
    S=S.*k; % scalar k
    Sp=zeros(33,33);
%     Sp(1:3,1:3)=ones(3,3).*norm(S(1:3,1:3))./100;
    b=zeros(33,1);
    b(4)=1;
    b(8)=1;
    b(12)=1;
    c=zeros(33,1);
%     c(1)=1;

    % mu=transpose(b)*(I*(transpose(A)*Wi*y-S*xprior))./D;
    %  J0= transpose(y-A*xprior)*Wi*(y-A*xprior)+(xprior(4)+xprior(8)+xprior(12)).*mu;
    %  J=J0;
    % while J>J0*10^-2
    if nnz(up)==0
        I=pinv(S+transpose(A)*Wi*A);
        D=transpose(b)*(I*b);
                    
%         O(1,:) = [transpose(b)*(I*b), transpose(b)*(I*c)];
%         O(2,:) = [transpose(c)*(I*b), transpose(c)*(I*c)];
    
%         figure(3)
%         hold on
% %         J=transpose(y-A*xprior)*Wi*(y-A*xprior);
% %         J=transpose(b)*xprior;
% %         J=transpose(y-A*xprior)*Wi*(y-A*xprior)+transpose(b)*xprior;
% %         J=transpose(y-A*xprior)*Wi*(y-A*xprior)+transpose(xprior)*(S/k)*xprior+transpose(b)*xprior;
%         plot(0,J,'sk','markersize',5);
%         figure(2)
%         hold on
%         plot([1:12],xprior(1:12),'-k');
%         hold on
%         xlim([0 13]);
%         xticks([1:12]);
%         xticklabels({'u','v','w','dudx','dudy','dudz','dvdx','dvdy','dvdz','dwdx','dwdy','dwdz'});
        for it=1:itr

%             P=[transpose(b)*(I*(transpose(A)*Wi*y-S*xprior));transpose(c)*(I*(transpose(A)*Wi*y-S*xprior))-xprior(1)];
%             Mu = linsolve(O,P);
%             xprior=I*(transpose(A)*Wi*y-S*xprior-Mu(1)*b-Mu(2)*c);
            
        
            %     xest=I*(transpose(A)*Wi*y-S*xprior-transpose(b)*(I*(transpose(A)*Wi*y-S*xprior)).*b./D);
            
            mu=transpose(b)*(I*(transpose(A)*Wi*y+S*xprior))./D;
            xest=I*(transpose(A)*Wi*y+S*xprior-mu.*b);
% %              J= transpose(y-A*xest)*Wi*(y-A*xest)+transpose(xest-xprior)*S*(xest-xprior)+(xest(4)+xest(8)+xest(12)).*mu
%             %  S=S/1.005;
%             %   I=pinv(S+transpose(A)*Wi*A);
%             %   D=transpose(b)*(I*b);
% %             xprior=(1.5*xest+0.5*xprior)./2;

% J=transpose(y-A*xest)*Wi*(y-A*xest)+transpose(xest-xprior)*(S/k)*(xest-xprior)+transpose(b)*xest;
%  J=transpose(y-A*xest)*Wi*(y-A*xest);
%  J=transpose(b)*xest;
% J=transpose(xest-xprior)*S*(xest-xprior);

            xprior=xest;           
% xest=xprior;
            

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
%             xticklabels({'u','v','w','dudx','dudy','dudz','dvdx','dvdy','dvdz','dwdx','dwdy','dwdz'});
            % i=i+1;
        end
        % i
        
        %     count=;
        %      if cond(A)>1*10^6
        %           Sol3=zeros(30,1);
        %           Sol3(1:2)=[U(In(1)),V(In(1))];
        % %             Sol30(1:2)=[mean(U(In(1:3))),mean(V(In(1:3)))];
        %       end
    else
%             S(1:3,1:3)=zeros(3,3);
        I=pinv(S+Sp+transpose(A)*Wi*A);
        D=transpose(b)*(I*b);
        for it=1:itr
            mu=transpose(b)*(I*(transpose(A)*Wi*y+S*xprior-Sp*up))./D;
            xest=I*(transpose(A)*Wi*y+S*xprior-Sp*up-mu.*b);
            xprior=xest;
        end
        %     b=transpose(b);
        %     fun = @(x)transpose(y-A*x)*Wi*(y-A*x)+transpose(x-xprior)*S*(x-xprior);
        %     options = optimoptions(@fmincon,'ConstraintTolerance',0.0001,'TolX',1e-5)
        %     %     lb=[-0.35;-0.35;-1.5;-1.5;-1.5;-1.5;-12;-12;-12;-12;-12;-12];
        %     %     ub=[0.35;0.35;1.5;1.5;1.5;1.5;12;12;-12;-12;-12;-12];
        %     % x = fmincon(fun,xprior,[],[],b,0,lb,ub,[],options);
        %     [x,fval,exitflag,output,lambda] = fmincon(fun,xprior,[],[],b,0);
        %     Sol3=x;
        %     count=lambda.eqlin;
        % end
    end
    Sol3=xest;
    if flag==1
        Sol3(4:12)=Sol3(4:12).*scale;
            Sol3(13:30)=Sol3(13:30).*scale^2;
    end
    %     J=term1+term2+term3; 
 
end
%p=rxx-rxy*inv(ryy)*rxy';