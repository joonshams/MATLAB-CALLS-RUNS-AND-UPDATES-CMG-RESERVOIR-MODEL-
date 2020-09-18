function [Xa2] = modelspaceENKF(xf,xo,Hxf,Pm,sgx,sgz,ens)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%time=3;
%ens=40;
%xf=Xf(:,:,k,:);xo=Xo(:,:,k);Hxf=HXf(:,:,k,:);xa=Xa(:,:,k-1,:);sgx=sgx;sgz=sgz;ens=ens;
f=xf(:,:,1,:);
%a=xa(:,:,1,:);
clear xf
xf=permute(f,[2 1 4 3]);
%xa=permute(a,[2 1 4 3]);
% for i=1:1600
%     for j=1:ens
%    if i~=156 && i~=1516
%     xf(i,3,j)=0;
%    end
%    end
% end

for en=1:ens
temp=[xf(:,5,en);xf(:,6,en);xf(:,2,en);xf(:,1,en);xf(:,3,en);xf(:,4,en)];
%temp2=[xa(:,5,en);xa(:,6,en);xa(:,2,en);xa(:,1,en);xa(:,3,en);xa(:,4,en)];
xf2(:,en)=temp;
%xa2(:,en)=temp2;
end
presloc=find(xo(:,3));

XO=[xo(:,1)',xo(:,2)',xo(presloc(:),3)']';

mxf2=mean(xf2')';
exf2=xf2-mxf2;
m=length(xf2);
n=length(XO);
H=eye(n,m);
sum1=1;
for i=3201:3200+length(presloc)
    H(i,i)=0;
    H(i,3200+presloc(sum1))=1;
    sum1=sum1+1;
end
A=exf2'*H';
K=1/(ens-1)*exf2*A*inv((1/(ens-1)*(A')*A+eye(length(XO))));
Xa2=xf2+K*(XO-H*xf2);
%  ef=exf2*exf2';
%  K=1/(ens-1)*ef*H'*(H*ef*H'+eye(length(XO)))^(-1);
%  Xa2=xf2+K*(XO-H*xf2);
vp=permute(Xa2(1:1600,:),[3 1 2]);
qp=permute(Xa2(1601:3200,:),[3 1 2]);
p=permute(Xa2(3201:4800,:),[3 1 2]);
sg=permute(Xa2(4801:6400,:),[3 1 2]);
phi=permute(Xa2(6401:8000,:),[3 1 2]);
k=permute(Xa2(8001:9600,:),[3 1 2]);
qp(qp<0.01)=0.03+abs(normrnd(0,0.1));

clear Xa2
Xa2(1,:,:)=sg;
Xa2(2,:,:)=p;
Xa2(3,:,:)=phi;
Xa2(4,:,:)=k;
Xa2(5,:,:)=vp;
Xa2(6,:,:)=qp;

end

