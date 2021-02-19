load t;
load z;
thetas=0.441904762;
theta0=0.054;
H=8;
fun=@(x,xdata) GreenAmpt(x(1),x(2),H,xdata,thetas,theta0);
options=optimset('MaxIter',20000);
para=lsqcurvefit(fun,[1 10],z,t,[],[10 10],options);
plot(t,z,'x');
hold on
plot(fun(para,z),z)
hold off
A=10*10*pi/4;
z=Qc/A;
i=10*A*para(1)*(para(2)+z+H)./z;
i=para(1)*(para(2)+z+H)./z