close all 
clear all
clc;


% given parameters of the filter
A=12;
M1=20;

% x and y are initialized with zero
X=zeros(20,60);
Y=zeros(20,60);
% Variance of the unifornly distribyted white noise
variance=0.5;
V=sqrt(variance).*randn(1,size(Y,2));
V=sqrt(variance)*randn(20,60);%unifornly distribyted white noise
N=(0:59); %Creating variables with zero initialisation          
omega=2*(pi/M1);
theta=0;  

%N as 60 and number of samples as 20 for expectation calculation  
for i=1:20
    for j=1:60
        Y(i,j)=A*cos(omega*N(j)+theta);
        X(i,j)=Y(i,j)+V(i,j);
    end
end

figure
plot(N,Y)
ylabel('Value of Y');
title("Original Signal")

figure
plot(N,X);
ylabel('Value of X');
title("Noise added Signal");

%Calculating cross correlation of Y,X and auto-correlation of order 10
Ryx10=corrtn(Y,X,10);
Ry10=corrtn(Y,Y,10);
% w=inverse(R)*p for order M=10
h10=Ry10\Ryx10(:,1);

%Calculating cross correlation of Y,X and auto-correlation of order 15
Ryx15=corrtn(Y,X,15); 
Ry15=corrtn(Y,Y,15); 
% w=inverse(R)*p for order M=15
h15=Ry15\Ryx15(:,1); 

%Calculating cross correlation of Y,X and auto-correlation of order 20
Ryx20=corrtn(Y,X,20);
Ry20=corrtn(Y,Y,20); 
% w=inverse(R)*p for order M=20
h20=Ry20\Ryx20(:,1);  

 % ycap=summation of (w(i)*x(i) from i=1 to i=10,15,20)
ycap10=conv(X(1,:),h10);  
ycap15=conv(X(1,:),h15); 
ycap20=conv(X(1,:),h20);
s=size(X(1,:));  

% finding value of N for plotting N vs ycap vs y
figure
plot((1:s(1,2)),Y(1,:));
hold on
plot((1:s(1,2)),ycap10);   
legend('Y','Yest');   
hold off
xlabel('Value of N');
ylabel('Value of Y');
title('Order of FIlte: 10');


figure
plot((1:s(1,2)),Y(1,:));
hold on
plot((1:s(1,2)),ycap15);
legend('Y','Yest');
hold off
xlabel('Value of N');
ylabel('Value of Y');
title('Order of Filter: 15');


figure
plot((1:s(1,2)),Y(1,:));
hold on
plot((1:s(1,2)),ycap20);
legend('Y','Yest');
hold off
xlabel('Value of N');
ylabel('Value of Y');
title('Order of Filter: 20');

%this is for correlation function
function cq=corrtn(a,b,m)
[i,j]=size(b);temp=zeros(size(b));rab=zeros(1,m);Rab=zeros(m,m);
Xab=zeros(i,m);
for x=1:i
    for y=1:m
        for z=1:j-(y)
            if((y-1)~=0)
                temp(x,z+y-1)=b(x,z);
            else
                temp(x,:)=b(x,:);
            end
        end
        for z=1:y
            if((y-1)~=0)
                temp(x,y-1)=0;
            else
                temp(x,:)=b(x,:);
            end
        end
        Xab(x,y)=sum(a(x,:).*conj(temp(x,:)));
    end    
    rab(1,:)=rab(1,:)+Xab(x,:);
end
rab=rab/i;
for i=1:m      % Creating the R matrix from r0 to r(m-1)
    for j=1:m
        Rab(i,j)=rab(1,abs(i-j)+1);
    end
end
cq=Rab;
end

% this is function for convolution
function cq=conv(X,h)
h=h';
temp=zeros(size(X));
s1=size(h);s2=size(X);
%sum with one that spans all i in X(i)and other that spans all j in w(j)
for i=1:s2(1,2)    
    for j=1:s1(1,2)  
        if ((i+j-1)<=s2(1,2))
            temp(1,i)=temp(1,i)+(conj(h(1,j))*X(1,(i+j-1)));
        else
            temp(1,i)=temp(1,i)+(conj(h(1,j))*0);
        end
    end
end
cq=temp; % return summation of (w(j)*X(i+j-1)) as 
end  
