A = [0 1 0; 0 0 1; 0 0 0];
b = [0; 0; 1];

mu=1/2;
alpha=1;
gamma=4.5;

X=sdpvar(3,3);
y=sdpvar(1,3);
H=[-1-2*mu 0 0;0 -1-1*mu 0;0 0 -1];

F=[(A*X+X*A'+b*y+y'*b'+alpha*X)<=0];
F=F+[-gamma*X<=(X*H+H*X)];
F=F+[(X*H+H*X)<0];
F=F+[X>0];

solvesdp(F);
X=double(X);
y=double(y);

k=y/X;
P=inv(X);
