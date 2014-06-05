% Euler Beam eqn simulations
ccc

nEl=3;
pL = 20;
pN = 5;
[phi, D, xL, wL,gamma,lval] = GLL_Basis(pL);
d4L = D*D*D*D;
[phi, D, xN, wN,gamma,lval] = GLL_Basis(pN);
d4N = D*D*D*D;
%
kL=1;
kEps=0.0001;
kVec=zeros(pL+1+pN,1);
kVec(1:pL+1)=kL;
kVec((pL+1):(pL+1+pN))=kVec((pL+1):(pL+1+pN))+ones(pN+1,1)*kEps;

pL1=pL+1;
pN1=pN+1;
%% element construction
% elemental stiffness matrix which includes one linear and one 
% nonlinear element
Kel=zeros(pN1+pL1-1);
Kel(1:pL1,1:pL1)=d4L;
Kel(pL1:(pL1+pN),pL1:(pL1+pN))=d4N;
Wel(1:pL1)=wL;
Wel(pL1:(pL1+pN))=wN;
%% element abcissa construction
aL=0; bL=5; aN=0; bN=0.1;
xL_Mapped=(bL-aL)/2*(xL-1)+bL;
xN_Mapped=(bN-aN)/2*(xN-1)+bN;

xEl=zeros((pN1+pL),1);
xEl(1:pL1)=xL_Mapped; xEl(pL1:(pL1+pN))=xEl(pL1)+xN_Mapped;

%% Assembly Globally
N1=pL1+pN;

K=zeros(nEl*(N1-1)+1);
x=zeros(nEl*(N1-1)+1,1);
Coeff=eye(nEl*(N1-1)+1);
W=Coeff;
for k = 1:nEl
    for i = 1:N1
        rInd=(N1-1)*(k-1)+i;
        x(rInd)=xEl(i)+x((N1-1)*(k-1)+1);
        Coeff(rInd,rInd)=kVec(i);
        M(rInd,rInd)=Wel(i);
        for j = 1:N1
            cInd=(N1-1)*(k-1)+j;
            K(rInd,cInd)=Kel(i,j)+K(rInd,cInd);
        end
    end
end
% spy(K)
% figure
% spy(M)

u=sin(x/x(end));
plot(x,u)
% (M-dt/2*K)\(M+dt/2*K)