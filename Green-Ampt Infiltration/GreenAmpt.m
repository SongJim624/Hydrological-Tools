function [ t ] = GreenAmpt( Ks,sf,H,z,thetas,theta0 )
%GREENAMPT 此处显示有关此函数的摘要
%   此处显示详细说明
L=length(z);
t=zeros(L,1);
for i=1:L
    t(i)=(thetas-theta0)/Ks*(z(i)-(sf+H)*log((z(i)+sf+H)/(sf+H)));
end
end

