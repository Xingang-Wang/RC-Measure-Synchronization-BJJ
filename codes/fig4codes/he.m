clc,clear

load para2.mat opt_result
hyperpara_set=opt_result;
eig_rho = hyperpara_set(1);
W_in_a = hyperpara_set(2);
a = hyperpara_set(3);
reg = hyperpara_set(4);
density = hyperpara_set(5);

load sl.mat ph

rng(4)
Ydata=ph;
inSize = 5; 
outSize = 4;
indata=Ydata(:,1:inSize);
outdata=Ydata(:,1:outSize);
trainLen =9000;
testLen = 200000;
resSize =800; % size of the reservoir nodes;  

Win=(2.0*rand(resSize, inSize)-1.0)*W_in_a;
WW=sprand(resSize, resSize, density);
rhoW = abs(eigs(WW,1));

W = WW .* (eig_rho /rhoW); % the spectral radius is 0.1:
X = zeros(resSize,trainLen);
Yt = outdata(2:trainLen+1,:)';
x = (2.0*rand(resSize,1)-1.0)*0.5;
for t = 1:trainLen
    u = indata(t,:)';
    x = (1-a)*x + a*tanh( Win*u + W*x );
    X(:,t) = x;
    X(1:2:end,t) = X(1:2:end,t).^2; 
end
X(:,6000:6100)=[];
Yt(:,6000:6100)=[];
X(:,3000:3100)=[];
Yt(:,3000:3100)=[];
X(:,1:100)=[];
Yt(:,1:100)=[];

X_T = X';

Wout = Yt*X_T / (X*X_T + reg*eye(resSize));


C1=[];
eps=-0.1;
for i=1:140
xx=x;
sum1=0;
sum2=0;
A1=zeros(testLen,outSize);
discard=2000;
Y1= zeros(outSize,testLen);
u1=[0 0.2 0 0.4 eps]';
u1(inSize)=eps;
for t = 1:1 
    xx = (1-a)*xx + a*tanh( Win*u1 + W*xx ); 
    xxx=xx;
    xxx(1:2:end)=xxx(1:2:end).^2;
    y1 = Wout*xxx;
end   
for t = 1:discard 
    xx = (1-a)*xx + a*tanh( Win*u1 + W*xx ); 
    xxx=xx;
    xxx(1:2:end)=xxx(1:2:end).^2;
    y1 = Wout*xxx;
    u1(1:outSize) = y1;
    u1(inSize)=eps;
end

for t = 1:testLen 
   xx = (1-a)*xx + a*tanh( Win*u1 + W*xx );
    xxx=xx;
    xxx(1:2:end)=xxx(1:2:end).^2;
    y1 = Wout*xxx;
    Y1(:,t) = y1;
    u1(1:outSize) = y1;
    u1(inSize)=eps;
end
A1=Y1';
for j=1:testLen
    sum1=sum1+0.5*1.2*A1(j,2)^2-1.0*cos(A1(j,1))*(1.0-A1(j,2)^2)^0.5;
    sum2=sum2+0.5*1.2*A1(j,4)^2-1.0*cos(A1(j,3))*(1.0-A1(j,4)^2)^0.5;
end
sum1=sum1/testLen;
sum2=sum2/testLen;
dx=[eps sum1 sum2];
C1=[C1;dx];
eps=eps+0.001;
end

load model.txt
figure(1);
hold on;
plot(C1(:,1),C1(:,2),'r--^','MarkerSize',5);
plot(C1(:,1),C1(:,3),'r--^','MarkerSize',5);
plot(model(:,1),model(:,2),'k--o','MarkerSize',5);
plot(model(:,1),model(:,3),'k--o','MarkerSize',5);
xlabel('u_{ab}','FontName','Times New Roman','FontSize',26);
ylabel('h_{a,b}','FontName','Times New Roman','FontSize',26);
