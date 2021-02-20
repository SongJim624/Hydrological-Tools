function zf=GreAmpt(t,H)
global Ks;
global thitas;
global thitai;
global sf;
syms t1;
syms z;
equ=t1==(thitas-thitai)/Ks*(z-(sf+H)*log((z+sf+H)/(sf+H)));
z=solve(equ,z);
z=subs(z,t1,t);
zf=double(z);