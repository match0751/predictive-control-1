clear;
clc;
delt=2;%采样周期
N=150; %模型长度
p=2;m=2; %系统维数-2输入2输出
t=0:delt:delt*N;
P=24;%优化时域 
M=12; %控制时域
nn=N/delt+1; %在线计算步数
rd=0;%扰动

%传递函数建立
num={[5],[3];[6],[9]};
den={[3 1 3],[1 2 5];[2 1 7],[2 3 6]};
sys=tf(num,den);%模型传递函数
sys11=tf(num{1,1},den{1,1}); % 通用传递函数模型转换为 MPC 传递函数模型 
sys12=tf(num{1,2},den{1,2});
sys21=tf(num{2,1},den{2,1});
sys22=tf(num{2,2},den{2,2});

%建立阶跃响应模型
[y11,t11,x11]=step(sys11,t);
model{1,1}=y11(1:N);
[y12,t12,x12]=step(sys12,t);
model{1,2}=y12(1:N);
[y21,t21,x21]=step(sys21,t);
model{2,1}=y21(1:N);
[y22,t22,x22]=step(sys22,t);
model{2,2}=y22(1:N);

%模型预测
%建立动态响应矩阵 A
for i=1:p
  for j=1:m  
    for i1=1:M
      for i2=1:P
        if i2<i1
           B{i,j}(i2,i1)=0; 
        else
           B{i,j}(i2,i1)=model{i,j}(i2-i1+1);%未超出时赋值为对应的阶跃响应系数
        end
      end
    end
  end
  C{i}=cat(2,B{i,1:m}); % 联结为动态响应矩阵的行数
end
A=cat(1,C{1:p});
%按列连接为动态响应矩阵 由阶跃响应系数组成的 P*M阵
for i=1:p
  aa{i}=cat(2,model{i,1:m});
end
A0=cat(1,aa{1:p});
ym=zeros(p*N,1);
y=zeros(nn,p);
du=zeros(m,1);
uu=zeros(m,1);


%滚动优化
ywt=[1,1];%Q矩阵
uwt=[300,300];%R矩阵
%得到Q矩阵
Q=[];
for i=1:1:p
    qi=ywt(1,i);
    Qi=qi*eye(P);
    Q=blkdiag(Q,Qi);
end
%计算R矩阵
R=[];
for j=1:1:m
    pj=uwt(1,j);
    Pj=pj*eye(M);
    R=blkdiag(R,Pj);
end

L=zeros(m,m*M);  

for i=1:m
      L(i,M*(i-1)+1)=1;
 end
D=L*inv(A'*Q*A+R)*A'*Q; %控制向量计算公式

sp1=1;sp2=1;%预测输出设定值
w{1}=sp1*ones(P,1);
w{2}=sp2*ones(P,1);
W=cat(1,w{1:p});
%输出期望值

%反馈校正
%建立矩阵H,S,S0  
%计算反馈校正H矩阵
alpha=[1,1];%H矩阵
H=[];
for i=1:p
    h=alpha(1,i)*ones(N,1);
    H=blkdiag(H,h); % 生成以H，h为对角线元素的对角阵
end  

ss=zeros(p,p*N); 
for i=1:m
    ss(i,N*(i-1)+1)=1;
end
  S=zeros(N);
  
for i=1:N
    for j=1:N
       if j==i+1
          S(i,j)=1;
       end
    end
end
S(N,N)=1;

for i=1:p
    s{i}=S;
end
S0=blkdiag(s{1:p});

%在线计算初始状态设置
u1_1=0;u1_2=0;
u2_1=0;u2_2=0;
y1_1=0;y1_2=0;y1_3=0;y1_4=0;
y2_1=0;y2_2=0;
y3_1=0;y3_2=0;
y4_1=0;y4_2=0;
u1=0;u2=0;
z11=c2d(sys11,0.5*delt,'zoh');
[num1,den1]=tfdata(z11,'v');
z12=c2d(sys12,0.5*delt,'zoh');
[num2,den2]=tfdata(z12,'v');
z21=c2d(sys21,0.5*delt,'zoh');
[num3,den3]=tfdata(z21,'v');
z22=c2d(sys22,0.5*delt,'zoh');
[num4,den4]=tfdata(z22,'v');

%在线计算
for k=1:nn
y111=num1(1)*u1+num1(2)*u1_1+num1(3)*u1_2-den1(2)*y1_1-den1(3)*y1_2;
y111=y111/den1(1);
y121=num2(1)*u2+num2(2)*u2_1+num2(3)*u2_2-den2(2)*y2_1-den2(3)*y2_2;
y121=y121/den2(1);
y131=num3(1)*u1+num3(2)*u1_1+num3(3)*u1_2-den3(2)*y3_1-den3(3)*y3_2;
y131=y131/den3(1);
y141=num4(1)*u2+num4(2)*u2_1+num4(3)*u2_2-den4(2)*y4_1-den4(3)*y4_2;
y141=y141/den4(1);

y1(k)=y111+y121+rd;
y2(k)=y131+y141+rd;

yn=ym+A0*du; %预测值
e=[y1(k),y2(k)]'-ss*yn; %预测误差
ycor=yn+H*e; %以对 e 加权的方式修正对未来输出的
ym=S0*ycor;%初始预测值
for i=1:p
    pp{i}=ym(N*(i-1)+1:N*(i-1)+P);
end
yp0=cat(1,pp{1:p});
du=D*(W-yp0);
uu=uu+du;
u(:,k)=uu;

%依次赋值，进行下一步运算
u1_2=u1_1;u1_1=u1;u1=u(1,k);
u2_2=u2_1;u2_1=u2;u2=u(2,k);
y1_2=y1_1;y1_1=y111;
y2_2=y2_1;y2_1=y121;
y3_2=y3_1;y3_1=y131;
y4_2=y4_1;y4_1=y141;
end


%控制量和输出量的变化曲线
y=[y1;y2];
x=0:delt:N;
subplot(211);
plot(x,y(1,:),'r');
hold on;
plot(x,y(2,:),'b');
legend('Y1','Y2' )
title('输出量 Y1、Y2 变化曲线');

subplot(212);
stairs(x,u(1,:),'r');
hold on;
stairs(x,u(2,:),'b');
legend('U1','U2' )
title('控制量 U1、U2 变化曲线');



