function cost = MY_costfunction(x, State_Initial, Np, Nc, T, Xref, Yref, PHIref, Q, R)
cost = 1;
l = 1;
%================================================
%
%================================================
X=State_Initial(1,1);
Y=State_Initial(2,1);
PHI=State_Initial(3,1);

X_predict=zeros(Np,1);
Y_predict=zeros(Np,1);
PHI_predict=zeros(Np,1);

X_error=zeros(Np+1,1);
Y_error=zeros(Np+1,1);
PHI_error=zeros(Np+1,1);

v=zeros(Np,1);
delta_f=zeros(Np,1);

for i = 1:1:Np
    if i == 1
        v(i,1)=x(1);
        delta_f(i,1)=x(2);
        X_predict(i,1)=X+T*v(i,1)*cos(PHI);
        Y_predict(i,1)=Y+T*v(i,1)*sin(PHI);
        PHI_predict(i,1)=PHI+T*v(i,1)*tan(delta_f(i,1))/l;
    else
        v(i,1)=x(3);
        delta_f(i,1)=x(4);
        X_predict(i,1)=X_predict(i-1)+T*v(i,1)*cos(PHI_predict(i-1));
        Y_predict(i,1)=Y_predict(i-1)+T*v(i,1)*sin(PHI_predict(i-1));
        PHI_predict(i,1)=PHI_predict(i-1)+T*v(i,1)*tan(delta_f(i,1))/l;
    end
    X_real=zeros(Np+1,1);
    Y_real=zeros(Np+1,1);
    X_real(1,1)=X;
    X_real(2:Np+1,1) = X_predict;
    Y_real(1,1)=Y;
    Y_real(2:Np+1,1) = Y_predict;
    X_error(i,1)=X_real(i,1)-Xref(i,1);
    Y_error(i,1)=Y_real(i,1)-Yref(i,1);
    PHI_error(i,1)=PHI_predict(i,1)-PHIref(i,1);
end

cost = cost+Y_error'*R*Y_error + X_error'*Q*X_error;

end