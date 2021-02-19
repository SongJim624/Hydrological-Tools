function [ t ] = GreenAmpt( Ks,sf,H,z,thetas,theta0 )
%GREENAMPT �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
L=length(z);
t=zeros(L,1);
for i=1:L
    t(i)=(thetas-theta0)/Ks*(z(i)-(sf+H)*log((z(i)+sf+H)/(sf+H)));
end
end

