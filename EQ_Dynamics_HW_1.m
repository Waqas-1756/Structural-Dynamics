clear all
n=2;
uo=[1;2];
uo_dot=[0;0];
zeta=[0.01;0.02];
Ic1=723;
Ic2=425;
l1=15*12;
l2=12*12;
E=29e06;
k1=((12*E*Ic1)/(l1^3))*2;
k2=((12*E*Ic2)/(l2^3))*2;
K=[k1+k2 -k2;-k2 k2];
m1=(3000*20)/390;
m2=(1500*20)/390;
m=[m1 m2];
M=diag(m);
t=0:0.01:3.5;
Ag=0.5*390;
[phi,w]=eig(K,M);
phin=(phi)./(phi(1,:));    
wn=sqrt(diag(w));
Mmod=phin'*M*phin;
Kmod=phin'*K*phin;
fo=M*Ag;
foo=diag(fo);
IC=[.5/wn(1,:) .5*wn(1,:);.5/wn(2,:) .5*wn(2,:)];
Cons=inv(IC)*zeta;
C=(Cons(1,:)*Mmod)+(Cons(2,:)*Kmod);
Cmod=diag(C);
for r=1:n
    phinn=phin(:,r);
    wnr=wn(r,:);
end

for r=1:n
    qro_dot=(phinn(:,r)'*M*uo_dot)/(phinn(:,r)'*M*phinn(:,r));
    qro=(phinn(:,r)'*M*uo)/(phinn(:,r)'*M*phinn(:,r));
    wnr=wn(r,:);
    qr(r,:)=(qro_dot/wnr)*sin(wnr*t)+qro*cos(wnr*t);
end
f=phinn*qr;

for r=1:n;
    wnr=wn(r,:);
    zetar=zeta(r,:);
    wdr=wnr*(sqrt(1-zetar.^2));
    qro_dot=(phinn(:,r)'*M*uo_dot)/(phinn(:,r)'*M*phinn(:,r));
    qro=(phinn(:,r)'*M*uo)/(phinn(:,r)'*M*phinn(:,r));
    qrr(r,:)=(exp(-zetar*wnr*t)).*(((qro_dot+(zetar*wnr*qro))/wdr)*sin(wdr*t)+(qro*cos(wdr*t)));
end
f2=phinn*qrr;


for r=1:n % Forced Damped System
    Kmodr=phin(:,r)'*K*phin(:,r);
    Mmodr=phin(:,r)'*M*phin(:,r);
    wnr=wn(r,:);
    w_bar=wnr;
    T_bar=(2*pi)/w_bar;
    zetar=zeta(r,:);
    Cmodr=Cmod(r,:);
    wdr=wnr*(sqrt(1-zetar^2));
    fmodr=phin(:,r)'*foo;
    qro_dot=(phin(:,r)'*M*uo_dot)/(phin(:,r)'*M*phin(:,r));
    qro=(phin(:,r)'*M*uo)/(phin(:,r)'*M*phin(:,r));
    ustr=fmodr/Kmodr;
    qsst=(ustr*wnr^2*sin(t*w_bar)*(-w_bar^2+Kmodr/Mmodr))/((-w_bar^2+Kmodr/Mmodr)^2 + 4*w_bar^2*wnr^2*zetar^2) - (Cmodr*ustr*w_bar*wnr^2*cos(t*w_bar))/(Mmodr*((- w_bar^2 + Kmodr/Mmodr)^2 + 4*w_bar^2*wnr^2*zetar^2));
    Ar=exp(-t*wnr*zetar);
    Br=(cos(t*wdr)*(qro+(Cmodr*ustr*w_bar*wnr^2)/(Mmodr*((-w_bar^2+Kmodr/Mmodr)^2 + 4*w_bar^2*wnr^2*zetar^2))) + (sin(t*wdr)*(qro_dot + wnr*zetar*(qro + (Cmodr*ustr*w_bar*wnr^2)/(Mmodr*((- w_bar^2 + Kmodr/Mmodr)^2 + 4*w_bar^2*wnr^2*zetar^2))) - (ustr*w_bar*wnr^2*(- w_bar^2 + Kmodr/Mmodr))/((- w_bar^2 + Kmodr/Mmodr)^2 + 4*w_bar^2*wnr^2*zetar^2)))/wdr);
    qtrans=Ar.*Br;
    q=qtrans+qsst;
end
f4=phinn*q;
