%[text] ## Equation electriques du PMSM
%[text] $\\frac{d}{dt} i\_{abc} = -r/L i\_{abc} +p \\psi\_f /L e\_{abc}(\\theta) \\omega\_{meca} + u\_{abc}$
%[text] $\\frac{d}{dt} \\omega\_\_{meca} = p \\psi\_f/J e\_{abc}(\\theta)^T i\_{abc} -B/J \\omega\_{meca} +T\_m$
%[text] $e\_a=\\cos(p \\theta\_{meca})$ 
%[text] $e\_b = e\_a(t-T/3)$
bemfm_a=-ScalarPhasorArray([0 1/2/1i],'isreal',true); %basic motiv of back emf, un per unit
bemfm_abc=bemfm_a.PhaseShift([0;-2*pi/3;2*pi/3]);
figure(1)
plot(bemfm_abc,title='back emf') %[output:2af1bf65]

%[text] parametre de la machine
model.r=0.082; %Ohm
model.L=1.5e-3; %H
model.J=0.01650; %kg m²
model.p=4; %pole pair
model.psi_f=0.176; %Wb
model.Bf=0.02; %frottement

model.v_in=600;

bemf_abc=bemfm_abc*model.psi_f*model.p;
%[text] construction de la matrice
A=[-model.r/model.L*eye(3)    ,-1/model.L*bemf_abc;
    1/model.J*bemf_abc'       ,-model.Bf/model.J];


figure(2)
plot(A,1,0:0.01:2,title='A(t)') %[output:7912476c] %[output:54145e5e]
C3=(eye(3)-ones((3))/3) %[output:911fa605]

B=([model.v_in*eye(3)/model.L;zeros(1,3)]);
B=PhasorArray(B);

model.A=A;
model.B=B;
model.state_name={'i_{a}','i_{b}','i_{c}','\omega_{m}'};


COS = ScalarPhasorArray([0 1/2],'isreal', true) %[output:667910d9]
COS3 = COS.PhaseShift([0 -(2*pi/3) (2*pi/3)]);
SIN3 = COS3.PhaseShift(-pi/2);
Abc2Dq0 = (sqrt(2/3))*[COS3;-SIN3;(sqrt(2)/2*ones(1,3))] %[output:968e98ab]

figure %[output:18075dc5]
plot(Abc2Dq0,2*pi) %[output:18075dc5]

figure %[output:39d73856]
plot(Abc2Dq0*Abc2Dq0.',2*pi) %[output:39d73856]
ytickformat('%f') %[output:39d73856]

value(reduce((Abc2Dq0)*bemf_abc)) %[output:70d065db]
%%
%[text] ## I. régulation de la vitesse moteur
%[text] ### Ia. etat augmenté contenant un intégrateur sur la vitesse mecanique

C=PhasorArray([0 0 0 1]);
C=[[0 0 0 1];[Abc2Dq0{1,:} 0]];
model.C=C;
model.output_name={'\omega_{m}','i_{d}'};

model.int_name=strcat("\int ",model.output_name);

Atilde=[[A;C],zeros(size(A,1)+size(C,1),size(C,1))];
Btilde=[B;zeros(size(C,1),size(B,2))];

nx=size(Atilde,1);
nu=size(Btilde,2);

Q=PhasorArray(diag([1,1,1,1,1000,30]))*1e4;
R=PhasorArray(eye(nu))*1e2;
N=PhasorArray(zeros(nu,nx));

% figure(3)
% plot([Atilde,Btilde],title='(A | B)')

% figure(4)
% stem([Atilde,Btilde])

pause(0.2)
%%
%[text] ### Ib. Synthèse du probleme en LMI harmoniques
hP=5; %nombre de phaseur dans P
hlmi=15; %ordre de troncature de la lmi
skip_non_poly_LMI = true %[output:0c5099ec]

%P une sdpvar dans la classe Phasor Array
P=PhasorArray.ndsdpvar(nx,nx,hP);
PT=P.TB(hlmi);

Y=PhasorArray.ndsdpvar(nu,nx,hP,"PhasorType",'full');
YT=Y.TB(hlmi);

P2=PhasorArray.ndsdpvar(nx,nx,hP);
P2T=P2.TB(hlmi);

X=PhasorArray.ndsdpvar(nu,nu,hP);
XT=X.TB(hlmi);

[AHpJ,AJHm,AHp,AHm]=Atilde.TBHankel(hlmi);
[~,AtJHm,AtHp,~]=TBHankel(Atilde.',hlmi);

[YHpJ,~,~,YHm]=Y.TBHankel(hlmi);
[YtHpJ,YtJHm,YtHp,YtHm]=TBHankel(Y.',hlmi);

[PHpJ,PJHm,PHp,PHm]=P.TBHankel(hlmi);
[P2HpJ,P2JHm,P2Hp,P2Hm]=P2.TBHankel(hlmi);
[BHpJ,BJHm,BHp,BHm]=Btilde.TBHankel(hlmi);

ATB=Atilde.TB(hlmi);
BTB=Btilde.TB(hlmi);

QT=Q.TB(hlmi);

RT=R.TB(hlmi);

NT=N.TB(hlmi);

RPM_bornes=[-5000 5000];
f_meca_bornes=RPM_bornes/60;
f_elec_bornes=model.p*f_meca_bornes;
T_elec_bornes=1./f_elec_bornes;

%[text] on voit le parametre $\\omega$ dans $\\mathcal{N$ comme un polytope, fonction des bornes
%[text] de vitesses max du moteur -\> on recherche un correcteur valable quelque
%[text] soit la fréquence
%[text] $\\left\[ \\matrix{A^TP+PA+Q & PB+N^T \\cr \n\* & R}\\right\]\>=0$
if ~skip_non_poly_LMI
    F0=[PT>=0];
    F=[F0];
    for T=T_elec_bornes
        F11=[ATB'*PT+PT*ATB+...
            AtHp*PHm+PHp*AHm+...
            AtJHm*PHpJ+PJHm*AHpJ+...
            NTB(Atilde,hlmi,T)*PT-PT*NTB(Atilde,hlmi,T)+QT];
        F12=[PT*BTB+PHp*BHm+PJHm*BHpJ+NT'];
        F21=F12';
        F22=R.TB(hlmi);
        F1=[F11,F12;F21,F22];
        F=[F;F1>=0];
    end


    F
    obj=-trace(PT);
    optimize(F,obj)
    PP=sdpval(P); %on extrait du PhasorArray P contenant des sdpvar, un phasorArray valant la valeur de P après optim
    K1=R\(Btilde.'*PP+N);

    clear F0 F
end
%[text] Mais ce n'est pas une LMI polytopique normalement.... 
%[text] Donc l'utilisation de ce K1 devrait échouer.
%[text] ### LMI LQR Polytopique, K commun
%[text] E. Feron, V. Balakrishan, S. Boyd, L. El. Ghaoui, "Numerical Methods for H\_2 related problems"
%[text] P. Xia, H. Shi,  H. Wen, Q. Bu, Y. Hu and Y. Yang, "Robust LMI-LQR Control for  Dual-Active-Bridge DC–DC Converters With High Parameter Uncertainties,"  in *IEEE Transactions on Transportation Electrification*, vol. 6, no. 1, pp. 131-145, March 2020, doi: 10.1109/TTE.2020.2975313.
%[text] C. Olalla, R. Leyva, A. El Aroudi and I. Queinnec, "Robust LQR Control for PWM Converters: An LMI Approach," in *IEEE Transactions on Industrial Electronics*, vol. 56, no. 7, pp. 2548-2558, July 2009, doi: 10.1109/TIE.2009.2017556.[![](text:image:792d)](https://udl.alma.exlibrisgroup.com/view/action/uresolver.do?operation=resolveService&package_service_id=61344629040005596&institutionId=5596&customerId=5595)
%[text] $\\text{minimize } \\quad \\gamma$
%[text] $\\text{subject to :}$
%[text] $P=P^T\>0$
%[text] $-Tr(QP) - Tr(X) + \\gamma \>0 $
%[text] $-A(\\omega\_i)P - P A(\\omega\_i)^T +BY +Y^TB^T -I \>0\n$
%[text] $\\left\[\\matrix{X & R^{1/2}Y \\cr Y^TR^{1/2} & P}\\right\] \>0\n$
%[text] $\\omega\_i$ the vertex of polytopic parameters
%[text] Then 
%[text] $K/P$ is the optimal state feedback gain
G=[P2T>=0];
for T=T_elec_bornes
    G11=[(ATB-NTB(Atilde,hlmi,T))*P2T+AHp*P2Hm+AJHm*P2HpJ];
    G12=[BTB*YT+BHp*YHm+BJHm*YHpJ];
    G1=-G11-G11'+G12+G12'-eye(size(P2T));
    G=[G;G1>=0];
end

sdpvar gam
Rb=reduce(R^(1/2));
G2=[X Rb*Y;(Y.')*Rb P2]; %[output:635343b0]
G=[G;... %[output:group:9c6b34b6] %[output:2850b3b5] %[output:163145a2]
    G2.TB(hlmi)>=0;... %[output:2850b3b5] %[output:163145a2]
    phas(-trace(Q*P2)-trace(X)+gam,0)>=0]; %[output:group:9c6b34b6] %[output:2850b3b5] %[output:163145a2]
%[text] Optimisation de la LMI et Calcul du correcteur
G %[output:4f20a75e]
obj2=gam;
optimize(G,obj2) %[output:2007d292] %[output:5b17b169]
PP=sdpval(P2); %on extrait du PhasorArray P contenant des sdpvar, un phasorArray valant la valeur de P après optim
K2=sdpval(Y/PP).trunc(20);

%[text] 
%%

clear G G11 G12 G1 G2 AHm AHp AHpJ AJHm ATB AtHp AtJHm P2T PT QT YT BHm BHp BHpJ BJHm BTB NT XP RT XT P2Hm P2Hp P2HpJ P2JHm PHm PHp PHpJ PJHm

figure(5)
% plot(K1)
% hold on
plot(K2,title="K(t)") %[output:7eb77596]

figure(6)
stem(K2) %[output:7eb77596]

%boucle fermée
% ABF=Atilde-Btilde*K1;
ABF = Atilde-Btilde*K2; %[output:0f0272c1]
figure(7)
% plot(ABF)
% hold on
plot(ABF)
figure(8)
% plot(ABF.HmqNEig(20,T),'*')
% hold on
plot(HmqNEig(ABF,20,T),'o')
%%
%[text] ## Simulation du systeme
%storage for simulation 
clear StoreSimu

StoreSimu.t = 0;
StoreSimu.y = 0;
StoreSimu.x0 = 0;
StoreSimu.Fs = 0;
StoreSimu.fode = 0;
StoreSimu.w_sc = 0;
StoreSimu.Bw = 0;
StoreSimu.K = 0;
StoreSimu.outmes = 0;
StoreSimu.u = 0;
StoreSimu.win = 0;
StoreSimu = StoreSimu([],1)

%%
clear StoreSimu_i
polyBase = @(tt) (3*tt.^2+-2*tt.^3).*(tt>0).*(tt<1) + 1.*(tt>=1);
polyBase2 = @(tt) (4*tt.^3+-3*tt.^4).*(tt>0).*(tt<1) + 1.*(tt>=1);

%input
Bw=PhasorArray([zeros(3,3);[1/model.J 0 0];[0 ,-1 0; 0 0 -1]]);
w_target.value = 2500*2*pi/60;
w_target.tIn = 1;
w_target.Dt = .1;
torqueIn0.value = -1;
torqueIn0.tIn = 3;
torqueIn0.Dt = 0.1;

torqueIn1.value = 7/8/10;
torqueIn1.tIn = 3;
torqueIn1.Dt = 0.1;
torqueIn1.dephase = 0.43;




w = @(y,tt) PhasorArray(...
    cat(3,...
    [torqueIn0.value*polyBase2((tt-torqueIn0.tIn)/torqueIn0.Dt); ...
    min(max(0,w_target.value*polyBase2((tt-w_target.tIn)/w_target.Dt)),w_target.value); ...
    0],...
    [polyBase2((tt-torqueIn1.tIn)/torqueIn1.Dt)*exp(1i*torqueIn1.dephase)*y(4)*torqueIn1.value; ...
    0; ...
    0]) ...
,'isreal',true);

w_sc = @(y,tt) cat(3, ...
    [-2*imag(polyBase2((tt-torqueIn1.tIn)/torqueIn1.Dt)*exp(1i*torqueIn1.dephase)*y(4)*torqueIn1.value) ...
    ;0;0], ...
    [torqueIn0.value*polyBase2((tt-torqueIn0.tIn)/torqueIn0.Dt); ...
    min(max(0,w_target.value*polyBase2((tt-w_target.tIn)/w_target.Dt)),w_target.value); ...
    0], ...
    [2*real(polyBase2((tt-torqueIn1.tIn)/torqueIn1.Dt)*exp(1i*torqueIn1.dephase)*y(4)*torqueIn1.value) ...
    ;0;0]...
    );

K=K2.trunc(2);
fod = @(t,y) fode(t,y,ABF,model.p,Bw,w);

% ABF=Atilde-Btilde*K;
phase_error=2*pi/100*0;


Atilde=[[A;C.PhaseShift(phase_error)],zeros(size(A,1)+size(C,1),size(C,1))];
Btilde=[B;zeros(size(C,1),size(B,2))];

ABF=Atilde-Btilde*K.PhaseShift(phase_error);

ABFsc = ABF.SinCosForm();
Asc = A.SinCosForm();
Bsc = B.SinCosForm();
Csc = SinCosForm(C.PhaseShift(phase_error));
Ksc = SinCosForm(K.PhaseShift(phase_error));
Bwsc = Bw.SinCosForm();
% w_sc = @(y,tt) w(y,tt).SinCosForm();

fod = @(t,y) fodeCosSin(t,y,ABFsc,model.p,Bwsc,w_sc);
fodd = @(t,y,Z) fodeCosSinWithDelay(t,y,Z,Asc,Bsc,Csc,Ksc,model.p,Bwsc,w_sc);
fodbis = @(t,y) fodeCosSinWithDelay(t,y,[y y],Asc,Bsc,Csc,Ksc,model.p,Bwsc,w_sc);
% fodeCosSinWithDelay(t,y,Z,Asincos,Bsincos,Csincos,Ksincos,p,Bwsincos,wsincos)

delays=[0.001 0.001];

history=[zeros(nx+1,1)];

Fs=40e3;
tspan = [0:1/Fs:8];

x0=[zeros(nx+1,1)];
x0(4)=0*2*pi/60;
full_state=[model.state_name model.int_name 'harmonic Phase'];

dFdy = @(t,y) jacobianCosSin(t,y,ABFsc,model.p);
opts = odeset('Jacobian',dFdy);

t1=tic();

[t,y]=ode15s(fod,tspan,x0,opts);
% [t,y]=ode15s(fod,tspan,x0);
% [t,y]=ode15s(fodbis,tspan,x0);
% [t,y]=ode15s(fodbis,tspan,x0,opts);
% sol = dde23(fodd,delays,history,tspan)

toc(t1)


theta_sim=y(:,end)';

Csct=evalSC(Csc,theta_sim);
z=squeeze(pagemtimes(Csct,permute(y(:,1:4),[2,3,1])));

Csct_real=evalSC(SinCosForm(C),theta_sim);
z_real=squeeze(pagemtimes(Csct_real,permute(y(:,1:4),[2,3,1])));


Ksct=evalSC(Ksc,theta_sim);
u=squeeze(pagemtimes(Ksct,permute(y(:,1:end-1),[2,3,1])));


Win=zeros(3,numel(t));
for ti = 1:numel(t)
    Win(:,ti)=evalSC((w_sc(y(ti,1:4),t(ti))),y(ti,end));
end
 
toc(t1)

StoreSimu_i.t = t;
StoreSimu_i.y = y;
StoreSimu_i.x0 = x0;
StoreSimu_i.Fs = Fs;
StoreSimu_i.fode = fod;
StoreSimu_i.w_sc = w_sc;
StoreSimu_i.Bw = Bw;
StoreSimu_i.K = K;
StoreSimu_i.outmes = z;
StoreSimu_i.u = u;
StoreSimu_i.win = Win;
toc(t1)

StoreSimu(end+1)=StoreSimu_i

%%

f1=figure;
clf(f1)
T1=tiledlayout(f1,size(Atilde,1)+1,1,"TileSpacing","compact","Padding","compact");

f2=figure;
clf(f2)
T2=tiledlayout(f2,size(C,1),1,"TileSpacing","compact","Padding","compact");

f3=figure;
clf(f3)
T3=tiledlayout(f3,size(K2,1),1,"TileSpacing","compact","Padding","compact");

    for ii=1:size(y,2)
        nexttile(T1,ii);
        hold on
        plot(t,y(:,ii))
        grid on
        title(full_state{ii})
    end
    sgtitle('x(t)')
    linkaxes(T1.Children,'x')
    pause(1)

    for ii=1:size(z,1)
        nexttile(T2,ii)
        hold on
        plot(t,z(ii,:))
        plot(t,z_real(ii,:))
        grid on
        title(model.output_name{ii})
    end
    sgtitle('C(t) x')
    linkaxes(T2.Children,'x')
    pause(1)

    for ii=1:size(u,1)
        nexttile(T3,ii)
        hold on
        plot(t,u(ii,:))
        grid on
        %         title(model.output_name{ii})
    end
    sgtitle('Kx')
    pause(1)
    linkaxes(T3.Children,'x')


f4=figure(12);
clf
T4=tiledlayout(f4,size(Win,1),1,"TileSpacing","compact","Padding","compact");

    for ii=1:size(Win,1)
        nexttile(T4,ii)
        hold on
        plot(t,Win(ii,:))
        grid on
        %         title(model.output_name{ii})
    end


theta_hmq = y(:,end);
omeg_mot = y(:,4);
omeg_hmq = model.p*y(:,4);
domeg_mot=gradient(omeg_mot,1/Fs);
domeg_hmq=gradient(omeg_hmq,1/Fs);

figure(13)
clf
T=tiledlayout(4,1);
nexttile   
plot(t,y(:,4))
hold on
plot(t,z(1,:),'r--')
hold off
title('\omega_m(t)')

grid on
nexttile   
plot(t,domeg_hmq)
title('D_t \omega_h(t)')
grid on
nexttile   
plot(t,domeg_hmq./omeg_hmq)
title('D_t( \omega_h(t) )/ \omega_h')
grid on
nexttile   
plot(t,domeg_hmq./omeg_hmq./omeg_hmq)
title('D_t( \omega_h(t) )/ \omega_h^2')
grid on

linkaxes(T.Children,'x')


figure(131)
TT=tiledlayout(3,1);
nexttile
plot(theta_hmq,omeg_mot)
grid minor
title('\omega_m')
nexttile
plot(theta_hmq,omeg_hmq)
grid minor
title('\omega_{elec}')
nexttile
plot(theta_hmq,gradient(omeg_hmq,theta_hmq))    
grid minor
title('d \omega_h / d \theta_h')

xlabel('phase (\theta_{elec})')
linkaxes(TT.Children,'x')
xlim([0 4*pi])
%%

figure(14)
I=find(model.p*y(:,4)>0,1)-1;

angularsft(theta_hmq(I:end),t(I:end),omeg_hmq(I:end), ...
    {y(I:end,1),y(I:end,4),domeg_mot(I:end),domeg_mot(I:end)./y(I:end,4)}, ...
    0:5,{'i_{a}','\omega_m','D_t \omega_m','D_t\omega_m / \omega_m'},...
    [1 1 0 1 1])

% figure
% 
% angularsft(y(I:end,end),t(I:end),model.p*y(I:end,4),{y(I:end,1),y(I:end,4),domeg(I:end),domeg(I:end)./y(I:end,4)},0:5,{'i_{a}','\omega_m','D_t \omega_m','D_t\omega_m / \omega_m'},[1 1 0 1 1])

angularsft(y(I:end,end),t(I:end),model.p*y(I:end,4),...
    {y(I:end,4)},0:5,{'\omega_m'},[1 1 0 1 1])
%%

angularsft(y(I:end,end),t(I:end),model.p*y(I:end,4),...
    {domeg_mot(I:end)./y(I:end,4)},0:5,{'D_t\omega_m / \omega_m'},...
    [1 1 0 1 1])

%%
I=1
angularsft(y(I:end,end),t(I:end),model.p*y(I:end,4),...
    {y(I:end,1),y(I:end,4),domeg_mot(I:end),domeg_mot(I:end)./y(I:end,4)},1:5,{'i_{a}','\omega_m','D_t \omega_m','D_t\omega_m / \omega_m'},...
    [1 1 0 1 1])
%%
%[text] ### Verif de la stab
%[text] $\\mathcal{A-N(\\omega\_0)-BK}$ se balade dans un polytope entre les deux vitesse moyennes maximale du moteur. K est indep de $\\omega\_0$. 
%[text] Pour être exact :
%[text] $\\mathcal{A-N(\\omega\_0)-BK} = \\mathcal{A-\\omega\_0 \\times N(1)-BK} =  \\mathcal{A}\_0 + \\omega\_0  \\mathcal{A}\_1$
%[text] C'est à dire que la matrice harmonique A dépend de manière affine du parametre, et donc se balade dans un polytope de dim 1. 
%[text] Si on trouve une fonction de lyapunov valable sur tout ce polytope, alors on a prouvé la stabilité. 
%
% hP=20
% hlmi=10
% P2=PhasorArray.ndsdpvar(size(ABF,1),size(ABF,1),hP)
% PT=P2.TB(hlmi)
%
% ABFT=ABF.TB(hlmi)
% [AHpJ,AJHm,AHp,AHm]=ABF.TBHankel(hlmi);
% [AtHpJ,AtJHm,AtHp,AtHm]=TBHankel(ABF',hlmi);
% [PHpJ,PJHm,PHp,PHm]=P2.TBHankel(hlmi);
% F=[PT>=0]
%
%
% RPM_bornes=[-5000 5000]
% f_meca_bornes=RPM_bornes/60;
% f_elec_bornes=p*f_meca_bornes;
% T_elec_bornes=1./f_elec_bornes;
%
%
% for T=T_elec_bornes
%     base_lmi=(ABFT-NTB(ABF,hlmi,T))'*PT+PT*(ABFT-NTB(ABF,hlmi,T))
%     cor_lmi=AtHp*PHm+PHp*AHm+AtJHm*PHpJ+PJHm*AHpJ
%     F11=base_lmi+cor_lmi
%     F=[F;F11<=0]
% end
%
% optimize(F)
%
%
% P2=P2.sdpval
% figure
% plot(P2)
% hold on
% plot(PP)
%
% figure
% stem(P2)

%%
%[text] ## Correcteur gain schedulé
%[text] Tran, G. Q. B., Pham, T.-P., Sename, O., Gáspár, P. (2023) "Design of an LMI-based Polytopic LQR Cruise Controller for an Autonomous Vehicle towards Riding Comfort", Periodica Polytechnica Transportation Engineering, 51(1), pp. 1–7. [https://doi.org/10.3311/PPtr.20075](https://doi.org/10.3311/PPtr.20075)
%[text] 
hP=5;
hlmi=10

BTB=Btilde.TB(hlmi);
Cz=[PhasorArray(Q^(1/2));zeros(size(R,1),size(Q,2))];
Dz=[zeros(size(Q,1),size(R,2));PhasorArray(R^(1/2))];
Bw=PhasorArray.eye(Atilde.size(1));
BwT=Bw.TB(hlmi);

P3=PhasorArray.ndsdpvar(Atilde.size(1),Atilde.size(1),hP);

ATB=Atilde.spTB(hlmi);
P3T=P3.TB(hlmi);

[AHpJ,AJHm,AHp,AHm]=Atilde.spTBHankel(hlmi);
[AtHpJ,AtJHm,AtHp,AtHm]=spTBHankel(Atilde.',hlmi);

CzT=Cz.spTB(hlmi);
[CzHpJ,CzJHm,CzHp,CzHm]=Cz.spTBHankel(hlmi);
[CztHpJ,CztJHm,CztHp,CztHm]=spTBHankel(Cz.',hlmi);

DzT=Dz.spTB(hlmi);
[DzHpJ,DzJHm,DzHp,DzHm]=Dz.spTBHankel(hlmi);
[DztHpJ,DztJHm,DztHp,DztHm]=spTBHankel(Dz.',hlmi);

[P3HpJ,P3JHm,P3Hp,P3Hm]=P3.spTBHankel(hlmi);

[BHpJ,BJHm,BHp,BHm]=Btilde.spTBHankel(hlmi);

gam = sdpvar;
G=[gam>=0,P3T>=0]
Nv=1; %number of parameter
clear cY cW YT WT YHpJ YJHm YHp YHm
for ii = 1:2^Nv
    T=T_elec_bornes(ii);
    cY{ii} = PhasorArray.ndsdpvar(Btilde.size(2),Btilde.size(1),hP,"PhasorType",'full');
    cW{ii} = PhasorArray.ndsdpvar(Btilde.size(2)+Atilde.size(1),Btilde.size(2)+Atilde.size(1),hP,"PhasorType",'symmetric') ;
    YT{ii} = cY{ii}.spTB(hlmi);
    WT{ii} = cW{ii}.spTB(hlmi);
    [YHpJ{ii},YJHm{ii},YHp{ii},YHm{ii}]=cY{ii}.spTBHankel(hlmi);

    G11=[(ATB-NTB(Atilde,hlmi,T))*P3T+AHp*P3Hm+AJHm*P3HpJ];
    G12=[BTB*YT{ii}+BHp*YHm{ii}+BJHm*YHpJ{ii}];
    G1=G11+G11'+G12+G12'+BwT*BwT';
    G=[G;G1<=0];

    G2=[P3T,(CzT*P3T+CzHp*P3Hm+CzJHm*P3HpJ+DzT*YT{ii}+DzHp*YHm{ii}+DzJHm*YHpJ{ii})';...
        (CzT*P3T+CzHp*P3Hm+CzJHm*P3HpJ+DzT*YT{ii}+DzHp*YHm{ii}+DzJHm*YHpJ{ii}), WT{ii}];
    G=[G;G2>=0];
    G=[G;trace(WT{ii})<=gam];
end
obj = gam;
G
optimize(G,obj)
%%
figure(14)
stem(P3)
figure(15)
plot(sdpval(P3))

figure(16)
for ii=1:2
    hold on
    cK{ii}=trunc(sdpval(cY{ii})/sdpval(P3),10);
    plot(cK{ii})
end
%%
figure(17)
clf

for ii=1:2
    hold on
    plot(-cK{ii},linetype='-.')
end
% hold on
% plot(K1)
hold on
plot(K2)

bound_sch=2*pi./T_elec_bornes;

Bw=PhasorArray([zeros(3,3);[1/model.J 0 0];[0 ,-1 0; 0 0 -1]]);
fod = @(t,y) fode_QLPV(t,y,Atilde,Btilde,model.p,Bw,w,cK{1}.trunc(2),cK{2}.trunc(2),bound_sch);
Fs=20e3;


[t,y]=ode15s(fod,[0:1/Fs:20],[zeros(nx+1,1)]);

for ii=1:size(y,2)
    nexttile(T1,ii);
    hold on
    plot(t,y(:,ii))
    grid on
    title(full_state{ii})
end
sgtitle('x(t)')

Ctime=arrayfun(@(ang) real(evalp(C,ang)),y(:,end),'UniformOutput',false);
Ctime=cat(3,Ctime{:});
z=squeeze(pagemtimes(Ctime,permute(y(:,1:4),[2,3,1])));
for ii=1:size(z,1)
    nexttile(T2,ii)
    hold on
    plot(t,z(ii,:))
    grid on
    title(model.output_name{ii})
end
sgtitle('C(t) x')

alph= @(om) (om-bound_sch(1))/(bound_sch(2)-bound_sch(1))

KK1=value(trunc(cK{1},2));
KK2=value(trunc(cK{2},2));
Ktime=arrayfun(@(om,ang) real(evalp(PhasorArray(KK1*(1-alph(om))+alph(om)*KK2),ang)),y(:,4)*model.p,y(:,end),'UniformOutput',false);
Ktime=cat(3,Ktime{:});
u=squeeze(pagemtimes(Ktime,permute(y(:,1:end-1),[2,3,1])));
for ii=1:size(u,1)
    nexttile(T3,ii)
    hold on
    plot(t,u(ii,:))
    grid on
    %         title(model.output_name{ii})
end
sgtitle('Kx')

%%
figure(18)
XX=angularsft(y(:,end),t,y(:,4)*model.p,y(:,1:4),0:4)
figure(19)
ZZ=angularsft(y(:,end),t,y(:,4)*model.p,z,0:4)
figure(20)
SS=angularsft(y(:,end),t,y(:,4)*model.p,u,-4:4)
%%

SSp=PhasorArray(permute([SS{1}(:,end) SS{2}(:,end) SS{3}(:,end)],[2 3 1]))
figure(21)
plot(SSp)

%%
%[text] ## Ajout de precompensateur
Pr1=-(C*(1/(A-B*K2{:,1:4}))*B*[1 0;0 1;-1 -1]);
figure(22)
plot(Pr1);
figure(23)
stem(Pr1);
Pr1=Pr1.trunc(50);
Pr0=1/Pr1;
Pr0=Pr0.trunc(5);

S0=Pr0*[2500*2*pi/60;0]
figure(24)
plot([1 0;0 1;-1 -1]*S0);
hold on
plot(SSp);

figure(25)
stem([1 0;0 1;-1 -1]*S0,'explosed',true,'display','both')
figure
stem(SSp,'explosed',true,'display','both')


figure(26)
lsim((A-B*K2{:,1:4}),0:1/5000:0.2,[],1/(50),[],B*[1 0;0 1;-1 -1]*S0,checkreal=true);


figure(27)
plot(Pr0)

Bw2=Bw;
Bw2{1:4,2:3}=B*[1 0;0 1;-1 -1]*Pr0;
N=[1 0;0 1;-1 -1]*Pr0;
figure(28)
plot(Bw)
hold on
plot(Bw2)

F10=figure(29);
T10=tiledlayout(A.size(1)+C.size(1)+1,1);

F11=figure(30);
T11=tiledlayout(C.size(1),1);

F12=figure(31);
T12=tiledlayout(B.size(2),1);
x0=zeros(7,1)
%%

simandplot_static(model,K2,N*0,Bw,w,Fs,10,x0,T10,T11,T12);
simandplot_static(model,K2,N  ,Bw,w,Fs,10,x0,T10,T11,T12);

%%
P=sqrt(2/3)*[d(bemfm_abc)';bemfm_abc'; [1 1 1]*sqrt(2)/2] %abc vers dq0
PP=[P{1:2,:} [0 ;0];0 0 0 1]
Pdq=P{1:2,:};

Pinv=trunc(1/P,1)
Pinv.Value

PPinv=[Pinv{:,1:2} [0;0;0];[0 0] 1]
figure(32)
plot(Pinv)
%%

K_x=trunc(K2{:,1:4},5)
K2_z=trunc(K2{:,5:6},5)

figure(33)
plot(Pdq*K_x*PPinv)
hold on
plot(trunc(Pdq*K_x*PPinv,0))

figure(60)
clf
plot(Pdq*K2_z)

%%
KK=cK{1};


K_x=trunc(KK{:,1:4},5)
K_z=trunc(KK{:,5:6},5)

figure(34)
plot(Pdq*K_x*PPinv)
hold on
plot(trunc(Pdq*K_x*PPinv,0))

figure(60)
hold on
plot(trunc(Pdq*K_z,1))

KK=cK{2};


K_x=trunc(KK{:,1:4},5)
K_z=trunc(KK{:,5:6},5)

figure(25)
plot(Pdq*K_x*PPinv)
hold on
plot(trunc(Pdq*K_x*PPinv,0))

figure(60)
hold on
plot(trunc(Pdq*K_z,1))

%%
function [t,y,z,s] = simandplot_static(model,K,N,Bw,w,Fs,tmax,x0,T1,T2,T3)
    
    Atilde=[[model.A;model.C],zeros(size(model.A,1)+size(model.C,1),size(model.C,1))];
    Btilde=[model.B;zeros(size(model.C,1),size(model.B,2))];
    full_state=[model.state_name model.int_name 'harmonic Phase'];
    
    Bw{1:4,2:3}=model.B*N;

    ABF=Atilde-Btilde*K;
    fod = @(t,y) fode(t,y,ABF,model.p,Bw,w);
    [t,y]=ode15s(fod,[0:1/Fs:tmax],x0);

    for ii=1:size(y,2)
        nexttile(T1,ii);
        hold on
        plot(t,y(:,ii))
        grid on
        title(full_state{ii})
    end
    sgtitle('x(t)')

    Ctime=arrayfun(@(ang) real(evalp(model.C,ang)),y(:,end),'UniformOutput',false);
    Ctime=cat(3,Ctime{:});
    z=squeeze(pagemtimes(Ctime,permute(y(:,1:4),[2,3,1])));
    for ii=1:size(z,1)
        nexttile(T2,ii)
        hold on
        plot(t,z(ii,:))
        grid on
        title(model.output_name{ii})
    end
    sgtitle('C(t) x')

    Ktime=arrayfun(@(ang) real(evalp(K,ang)),y(:,end),'UniformOutput',false);
    Ktime=cat(3,Ktime{:});
    sk=squeeze(pagemtimes(Ktime,permute(y(:,1:end-1),[2,3,1])));
    Nrtime=arrayfun(@(ang) real(evalp(N,ang)),y(:,end),'UniformOutput',false);
    Nrtime=cat(3,Nrtime{:});
    rtime=arrayfun(@(y,t) evalp(w(y,t),y),y(:,end),t,'UniformOutput',false);
    rtime=cat(3,rtime{:});
    sr=squeeze(pagemtimes(Nrtime,rtime(2:3,1,:)));
    s=sk+sr;
    for ii=1:size(s,1)
        nexttile(T3,ii)
        hold on
        plot(t,s(ii,:))
        grid on
        %         title(model.output_name{ii})
    end
    sgtitle('Kx')
end

function dy =fode(t,y,A,p,Bw,w)


dy=(evalp(A,y(end),"forceReal",true))*y(1:end-1)+(evalp(Bw,y(end),"forceReal",true))*(evalp(w(y,t),y(end),"forceReal",true));%-[0;0;0;0;200]-[0;0;0;50/0.05;0]*(t>0.5);
dy(end+1)=p*y(4); %On integre la fréquence electrique = p * w_meca, pour obtenir l'angle electrique
if ~iscolumn(dy)
    dy=dy';
end
end

function dy =fodeCosSin(t,y,Asincos,p,Bwsincos,wsincos)
winsincos = wsincos(y,t);
hA=(size(Asincos,3)-1)/2;
hB=(size(Bwsincos,3)-1)/2;
hw=(size(winsincos,3)-1)/2;

h=max([hA,hB,hw]);

theta= y(end);
eit = [sin((h:-1:1)'*theta);cos((0:h)'*theta)];


At=tensorprod(Asincos,double(eit((1:(2*hA+1))+h-hA)),3,1); %est un 3D array dont Mt(:,:,k) est M(t(k))
Bt=tensorprod(Bwsincos,double(eit((1:(2*hB+1))+h-hB)),3,1); %est un 3D array dont Mt(:,:,k) est M(t(k))
wt=tensorprod(winsincos,double(eit((1:(2*hw+1))+h-hw)),3,1); %est un 3D array dont Mt(:,:,k) est M(t(k))

dy = At*y(1:end-1) + Bt*wt;
% dy=(evalp(A,y(end),"forceReal",true))*y(1:end-1)+(evalp(Bw,y(end),"forceReal",true))*(evalp(w(y,t),y(end),"forceReal",true));%-[0;0;0;0;200]-[0;0;0;50/0.05;0]*(t>0.5);
dy(end+1)=p*y(4); %On integre la fréquence electrique = p * w_meca, pour obtenir l'angle electrique
if ~iscolumn(dy)
    dy=dy';
end
end


function dy =fodeCosSinWithDelay(t,y,Z,Asincos,Bsincos,Csincos,Ksincos,p,Bwsincos,wsincos)
winsincos = wsincos(y,t);

hA  = (size(Asincos,3)-1)/2;
hB  = (size(Bsincos,3)-1)/2;
hC  = (size(Csincos,3)-1)/2;
hK  = (size(Ksincos,3)-1)/2;
hBw = (size(Bwsincos,3)-1)/2;
hw  = (size(winsincos,3)-1)/2;

h=max([hA,hBw,hw,hB,hK,hC]);

theta_reel= y(end);
eit_real = [sin((h:-1:1)'*theta_reel);cos((0:h)'*theta_reel)];

theta_delay= Z(end,2);
eit_delay = [sin((h:-1:1)'*theta_delay);cos((0:h)'*theta_delay)];


At  = tensorprod(Asincos   ,double(eit_real((1:(2*hA+1))+h-hA))   ,3,1); %est un 3D array dont Mt(:,:,k) est M(t(k))
Bwt = tensorprod(Bwsincos  ,double(eit_real((1:(2*hBw+1))+h-hBw)) ,3,1); %est un 3D array dont Mt(:,:,k) est M(t(k))
Bt  = tensorprod(Bsincos  ,double(eit_real((1:(2*hB+1))+h-hB))  ,3,1); %est un 3D array dont Mt(:,:,k) est M(t(k))
wt  = tensorprod(winsincos ,double(eit_real((1:(2*hw+1))+h-hw))   ,3,1); %est un 3D array dont Mt(:,:,k) est M(t(k))

Kt  = tensorprod(Ksincos   ,double(eit_delay((1:(2*hK+1))+h-hK))   ,3,1);
Ct  = tensorprod(Csincos   ,double(eit_delay((1:(2*hC+1))+h-hC))   ,3,1);

nx = size(At,1);
nc = size(Ct,1);
% nu = size(Bt,1);

dy = [At;zeros(nc,nx)]*y(1:nx) + [Bt*Kt;[Ct , zeros(nc)]]*Z(1:nx+nc,1)+ Bwt*wt;
% dy=(evalp(A,y(end),"forceReal",true))*y(1:end-1)+(evalp(Bw,y(end),"forceReal",true))*(evalp(w(y,t),y(end),"forceReal",true));%-[0;0;0;0;200]-[0;0;0;50/0.05;0]*(t>0.5);
dy(end+1)=p*y(4); %On integre la fréquence electrique = p * w_meca, pour obtenir l'angle electrique

if ~iscolumn(dy)
    dy=dy';
end
end

function J = jacobianCosSin(t,y,Asincos,p)

h  = (size(Asincos,3)-1)/2;

theta_reel= y(end);
eit_real = [sin((h:-1:1)'*theta_reel);cos((0:h)'*theta_reel)];

At  = tensorprod(Asincos   ,double(eit_real)   ,3,1); %est un 3D array dont Mt(:,:,k) est M(t(k))

%d/dtheta A(theta)*y
% si A(theta) = sum a_sk sin(k th) + a_ck cos(kth)
% alors d/dth A(th) = sum  k a_sk cos(k th) - k a_ck sin(kth)
% donc on flip et on multiplie par -k:k;
AsincosDer=flip(Asincos,3);
K=permute((-h):(h),[1 3 2]);
AsincosDer=bsxfun(@times,AsincosDer,K);


AtDth  = tensorprod(AsincosDer ,double(eit_real)   ,3,1); %est un 3D array dont Mt(:,:,k) est M(t(k))
AtDth=AtDth*y(1:end-1);

J = [At AtDth;
    zeros(1, numel(y))];

J(end,4)=p;

end



function dy =fode_QLPV(t,y,A,B,p,Bw,w,K1,K2,bound_sch)
At=real(evalp(A,y(end)));
Bwt=real(evalp(Bw,y(end)));
Bt=real(evalp(B,y(end)));
wt= real(evalp(w(y,t),y(end)));

scheduling_var=p*wt(2);
scheduling_var=p*y(4);
alpha=max(-1,min(1,(+scheduling_var-bound_sch(1))/(bound_sch(2)-bound_sch(1))));

K=alpha*K2.value+(1-alpha)*K1.value;

Kt=real(evalp(PhasorArray(K),y(end)));
dy=(At+Bt*Kt)*y(1:end-1)+Bwt*wt;

% dy=real(evalp(A,y(end)))*y(1:end-1)+real(evalp(Bw,y(end)))*real(evalp(w(y,t),y(end)));%-[0;0;0;0;200]-[0;0;0;50/0.05;0]*(t>0.5);
dy(end+1)=p*y(4); %On integre la fréquence electrique = p * w_meca, pour obtenir l'angle electrique
if ~iscolumn(dy)
    dy=dy';
end
end

function Atheta  =  evalSC(ASC,theta)

h  = (size(ASC,3)-1)/2;
eit_real = [sin((h:-1:1)'*theta);cos((0:h)'*theta)];

Atheta  = tensorprod(ASC ,double(eit_real)   ,3,1);
end




%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright","rightPanelPercent":40}
%---
%[text:image:792d]
%   data: {"align":"baseline","height":15,"src":"data:image\/svg+xml;base64,PHN2ZyBpZD0iQ2FscXVlXzEiIGRhdGEtbmFtZT0iQ2FscXVlIDEiIHhtbG5zPSJodHRwOi8vd3d3LnczLm9yZy8yMDAwL3N2ZyIgdmlld0JveD0iMCAwIDUwIDI4LjE4Ij48cmVjdCB3aWR0aD0iNTAiIGhlaWdodD0iMjguMTgiIHJ4PSIxNC4wOSIgc3R5bGU9ImZpbGw6IzIxY2UyMSIvPjxwYXRoIGQ9Ik0xOC45MywxNi40NWEuMzIuMzIsMCwwLDAtLjU0LjIydjJoLTRhNS4wOCw1LjA4LDAsMCwxLDAtMTAuMTZoNS4yMlY1LjQ0SDE0LjQ0YTguMTQsOC4xNCwwLDAsMCwwLDE2LjI4aDR2MmEuMzIuMzIsMCwwLDAsLjU0LjIzbDMuNTEtMy41MWEuMzQuMzQsMCwwLDAsMC0uNDZaIiBzdHlsZT0iZmlsbDojZmZmIi8+PHBhdGggZD0iTTQyLjUyLDIxLjA3bC00Ljg0LTUuNjlhMCwwLDAsMCwxLDAsMCw1LjMsNS4zLDAsMCwwLDIuNTktNC45LDUuNDMsNS40MywwLDAsMC01LjQ4LTVIMjguNDhhMCwwLDAsMCwwLDAsMFYyMS42OWEwLDAsMCwwLDAsMCwwaDMuMDVsMCwwVjguNTVzMCwwLDAsMEgzNWEyLjI0LDIuMjQsMCwwLDEsMi4yMiwyLjUzLDIuMywyLjMsMCwwLDEtMi4zMiwySDMyLjRhLjM2LjM2LDAsMCwwLS4yOC41OWw2LjYyLDcuNzlhLjkzLjkzLDAsMCwwLC42OC4zMWgyLjgxQS4zOS4zOSwwLDAsMCw0Mi41MiwyMS4wN1oiIHN0eWxlPSJmaWxsOiNmZmYiLz48Y2lyY2xlIGN4PSIyNC4yIiBjeT0iOC4yIiByPSIyLjc2IiBzdHlsZT0iZmlsbDojZmZmIi8+PC9zdmc+","width":27}
%---
%[output:2af1bf65]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAdwAAAEfCAYAAADr33fvAAAAAXNSR0IArs4c6QAAIABJREFUeF7tnX+MFseZ58v2ZDVD7IQZg2IwILAhjqVEJPbaGWM4y1I2P5ecTkILwylreVHCrSJHtwsyPxyt16szw1ghJ7GyIv5gLW9WDOzpdLnFu3u2bOUcgzOy48TcZUMcsMGAwTqbGYJ9DOclntPTL8879RbV\/VS\/XV3969sSmmGqu+qpT1XXU\/XUU09fNTU1NaVwgQAIgAAIgAAI5ErgKijcXPkicxAAARAAARCICEDhoiOAAAiAAAiAQAACULgBIKMIEAABEAABEIDCRR8AARAAARAAgQAEoHADQEYRIAACIAACIACFiz4AAiAAAiAAAgEIQOEGgIwiQAAEQAAEQAAKF30ABEAABEAABAIQgMINABlF1IvAkSNH1P3336++973vqcHBwdwrF7q8NBXat2+f2rJlS\/TI+vXr1aZNm9I8jntBoFEEoHAb1dyorA8CoRVg6PJcGY2Pj6t169apNWvWqNWrV7s+hvtAoLEEoHAb2\/SoeLcEQivA0OW5cimrXK7y4z4QCE0ACjc0cZRXeQKsaGh1t3v3bnX69OlYkyqvAg8dOtSut830OjY2ptauXWu9x6bY2JS7cuVKtX37dtXX12flapY\/d+5c9cQTT6glS5ZE93P617\/+dfXqq6+q\/fv3R3\/nfHfu3Kl27doV\/W3p0qVRfQcGBpRuSqY0M9\/KNzIqAAI5EIDCzQEqsqw3AVaApGj37NkT7eOy4qLfeR\/TZnLlZx944IG2GZaVV1xepsLl+4eHhxNNufwcKU+WiZ7967\/+67bS1RUyl6\/Xj8uYnJxUmzdvjhqWFTxWuPXu56idfwJQuP6ZIseaE7ApTaoyrVI3btzYVmak3Pbu3dteFTKWkZGRaFVMiosUmW0flPKi+2hFefbs2baT1ptvvhk5KUnKlsqi5ykfXpXS30zFyeXrEwWbcqVnzfpA4da8o6N63glA4XpHigzrTiBO0SQ5EemrRt1ke+rUKdHjmZ+9+eab1QsvvODkDWxbcXO76IqT\/mYqfCjcuvdg1K8oAlC4RZFHuZUlQArwwQcfVI899lh7L5QqYypcc\/+U9271FS7t7dLeLZtzbVB0Zb1ixYpI6Sbdr8ui7x3refN+LBRuZbshBK8gASjcCjYaRC6WgOsK12bSZVMvm5TTKFw690uK0txLtdFwPbJjuw8r3GL7F0qvLwEo3Pq2LWqWEwGXPdx58+ZZFSMruAULFiTu4epK\/frrr+8wO8eVb1ZXX0nrXszSHjIUbk4dB9k2ngAUbuO7AACkJcAKj57jIzY2j2BSbHTMRj+GQ3+jYzb6cR7TS9lUeLZ9XlveZj1sMvHxI9MjWg9eAYWbtkfgfhBwIwCF68YJd4FAm0DcOVyb5zArWH6Y7qFLP5pD\/097DpdXyvSs7oVsNpPrOVwoXHRwEMifABRu\/oxRAgiAAAiAAAgoKFx0AhAAARAAARAIQAAKNwBkFAECIAACIAACULjoAyAAAiAAAiAQgAAUbgDIKAIEQAAEQAAEoHDRB0AABEAABEAgAAEo3ACQUQQIgAAIgAAIQOGiD4AACIAACIBAAAJQuAEgowgQAAEQAAEQgMJFHwABEAABEACBAASgcANARhEgAAIgAAIgEEThcjD0oaEhNTg4COogAAIgAAIg0DgCuStcPXi69NHsxtFHhUEABEAABBpDIFeFSyvbnTt3qlWrVqmNGzeqTZs2YYXbmK6FioIACIAACOgEclW4XBCvcqFw0flAAASqTODE+EVF\/+g6OdH6yf8\/MT6ZqmoLBvra9y8Y6I1+n9\/f+rl88cxUeeHmahAoVOHedNNN1aAEKUEABBpD4IMFd6tLs26J6vvhjFnqwxnXRz9t19UX3m3\/+eoLZ1Mxony5jLgHOX\/Km36nn72\/\/u+pysHNyQSee+45tWjRoiCYCle4b7zxRmxFjx07lghCSqeMpXukdJoUJMnoowwfeYSQU2LlUo+6yOmDhY88JJ4+yvCRRxnkpHpc8\/E50YqUVqcHjk6ok7RinZhetVIf5tXmgv5eNX+gV9FKlP5Gq88br5kQB2eJl0v6W7\/rj2TUV882WVlekvXuxf3RWHr3zTNFOSUZXN5llzyytrtLGdI9Uroko09NXGqF67Oi3eYVsjG6lZGeg5xZ6F35LHjWg+foy2+3FeuB18+1K0UKlBXq8suKihTqH3\/xNnGC7ZdMd7kt\/PSd6u9+9Iw6eLlOB49ORBmZdSTlSxOGTV9a2F1BGZ+qwnsUUkYoXKFDhWyMLH0bcmahB4Xrl154nrQSJOUT\/Tw60VY8vFolxUOKlZRqpGwv75maktbhPYpjwavhoTvmtFfCee8VV4FnSBmDKNy4l1mqqGQKkNJ9mEUkGX2U4SOPEHL64F0XOX2w8JGHxNNHGT7y8C0nKxU2C\/PKjhSprlx1heJSD99y2sY+SQ4p3cWiZebBzl40KTEnJGyObnGbdtaS5JDSu5HT5OVShnSPlC61uc\/JaOEKV9+wJjB08Qb2Cy+8oObNm9f+f9p0uv\/UqVNqxYoVUb7m8z7SuTFIZlv+Lun6s3r99fz4927T6bmsPKXny8KTGcW1R1V4p+Fp1tXs73HpcSz0SaCt75mTRLNv+E4nGf72xZPq\/1w9O9p3JQU792M9qqenJ1KwGwd7O951s3yTpS1df1e7TWee+k9bW2RJjxvL4saGpPv\/\/jdXRdUeebo19rKpnfaDl804GTt28ngSN7a69r+4d9X1eR\/9t1EKV3JI8jm7QF4gAALVIMDOQrQHyysy3TxMZtG8zaHVIOVPSmJOvE0FzCboovaB\/dXQnhMUbt6EkT8IgEBpCXx79HB74NdNxEN33FBamesoGCtg0wT96nfvqlV1oXAvN6dke5fS2SyRdMZKykNK91GGjzxCyBmiDLDoHMsk5lnTffDOkgfvxf7NT46pV9662DZpDt05J3Jw8rmnKLHKUg+91aRysqb7kFOSwVbGgaPnIse00ZfPRM5pZNL\/xl3zL++Z2wN1SOVkTffBAgoXCveKGWQZOqYkg4\/O7yOPEHKGKKMqLLqRc+Tp4x2m4qWf6FFf\/ux8lbSKlZhnTe+mHralXlY5pOd9yJm1DFK4P3jmsPqXs62jSGSJINOzaXaWysma7oMFFG6tDBaoDAiAgGmeZFMx9mKr3zdspudNX1p0WQmXfxugUQo3yUs5ycuUZzb0k03GuL\/Tyxs8msMjyQvZ9q6YK4M8nj946Ig68PZH2itZ0wSpr27yKB\/5T0fqC8WXInmZjm40qfqjT05FMwt9e68s7dMohVv20I6SycOHScNHHiHkDFEGWFR7D5f3ZEdfOtM2NdKxnXvnK7Vq+a2xy7kQfStEGei\/000cZ3bWz\/tKbSKl++ANhYs9XOzhJgQTl15CKd3HSxqijKrISQPrj8aOqNd+2xOtbthc\/PjQtIKVeEnpVWFRFTlD8NZZ2LYUou2EG\/5V3b10SaETMSjc6m9roAYgUHsCNIhSwARdyVL4RBzfqX3TZ6ogK1\/2do5zuMpUSIqHoXBTwMKtIAAC4QjYVip1O5cZjiZKIgLksc6RropwtmqUwkVox5YjATs1UAe0OYFlTad8pdCMWdNJxqyhNHkIyhIqkxkmOY1l5Sk974N3Gp5xTjFxLHRzn60uZvr+w+9H52TpJ69Ibv7oZHtfNkRoR71vmaZbk1U36Xrfsz3vkl610I5xDqfcf5NCN1Kb+wztSH1L35b4yuJe9a3Pz+wYH21OVy7912w7vX0bpXDhNNXqCtKeStb0qpRRFTml9vBRDx95ZJGTV7M\/\/OnJdsziuGM8UjlZ04tmwQO2VI+qyBmiHt2ySGtyluoipUPhwtICAiBQCAGbl7EtoEEhwqHQxhFg5UsmZ7askI9A3OcVuwEEhdsNNTwDAiDQNQFzYKOjGwhK0TVOPOiZQNpVb5rioXAv05JMAVJ6tyYNvbFClFEVOcFiumfUgYUtjnE34fnQfzuHd6lvZE2vCu+85DQnh7TX+6dfvDV21SvxhsKFwr1igiZ1Gik9r85vCirJIaVXRc4Q9ciLhb5aoDIojvE3710S+7m7EHWtSxl5tZnv9ywE77xZ6D4Gp89fihSuzWNeqmujFC5COyI0JUJzLorGU9PTOM3\/zXvNwY7Sycv4xydV+9wsrQyG13zOWrbtefpbGUPz2dhB\/k5HTJf+UeX2pTCiTx1+X+16afpDChROMskLmxVxoxQuPkCfZrcB94JAOgLmapb2ZvUIUOlyw90gUG4C5l4vneuVnKxqpXDHxsbU2rVro1YaHh5Wq1evbreYVFHJFCCl523S4IpIckjpVZEzRD3AonNAk5gnpfOH3CXvzixl+HoHJBmq0i+qImcI3kWzMCOhxQVokfSQzynGVVNTU63POORwjY+Pqw0bNqitW7dGuW\/btk3t2LFDDQwMRP+XKip1Cim96AZ3HYyqImcI3mDRvcK1Hen5k9uuVd\/56mcS326pXaV0tFn3bWZrmLrwLku\/4K9Wmd\/qZfaSHvKpGnNVuLS6HRkZUbt371Z9fX1q8+bNamhoSA0ODjopXJ8VRV4gUFcCONJT15ZFvUIQqJXCHR0dVdu3b4+4kcJdtmxZ26w896t\/pk7\/038OwRRlgEDtCBw4ek5RAHj+eAACVNSuiVGhAAQao3AH\/vzHEc7eX\/+D+r0TB9WP9++L\/s+eZWWI7SvFBub+kCX2L2Ipd3rpch+I89LNml4V3lR\/s\/\/R8YcTH1yn6Huzb7zzvvrUzN+p\/\/iHn1U3XjPR8e6wOS+OlW7uY876u5cmne5FLOVjbf5sEtZNw3pfzpputm0e\/+c2DRlL2db\/8ui\/pGDNK5TzbuEm5b995ueYpQeYxaGIahMIEeKu2oQgPQh0R6A2K9w0TlMYULrrLHiqvgRsTlB8zKG+tUbNQCAsgdooXMKmHwvas2dP22GK0uIqSt9H1D9OjEEmbAdEacUSwNnZYvmj9GYRqJXCTWo6qaLszq1\/KUJ37Yb7\/DRdsKg+C5sTlB4tx\/YuSe0upfP+nx5ByiwnRB51KcMHT7AI+y5Lesjn9CPXPVxJUKqoS2jHaz4+J\/LEJMVLF614KWKO6SgS52SD0IHZQwfyQEI\/wdMfT3KCOvD2R9TBoxORE9TtN\/a24xqn6c\/mvebAj\/RppyYel+KcmpDeItAUPo1SuGm9w3AUQprGIL0KBMwoODjSU4VWg4x1JACFe7lVk0wrNGD94JnDUbBqXvXaYmZK5pms6T5MSD7ykOpRlTKqImc3vE0nqLkf61EPfW1JFOs17pLKyZrug7ePPKR6VKWMqsgZgndVWEDhOihcvTHJyYpMcgden\/5SBO\/1Sh0ra7qPTuUjD6keVSmjKnKm4W1bzfKWSNLeaVVYVEXONG3W7QQILDrJScyzpvvgDYXbpd1CP1rEq14a2JYvntlljngMBLojwH1RnwjC2747lngKBPIkAIXrgW7SqtdD9sgCBK4gQEq2pWg7wy1i0ofOAgLlJdAohevipZzVK\/bvf3NV28NZ3zdL4wXKpgv6mVUePO\/Py7cM7cFe9D\/86cloVLnnllnRN2dD9i94IcMLWX8XTFMr+kd8\/2iUwk3yUvZt37dFs8p6ztHHHoKPPCRWVSmjKnIy768\/\/osO3wHdcU9qEym9KiyqImcI3mDRuZKVmGdN98EbCjdn64MeyYd+lz7QnbM4yL4iBNjL+MDRifYXeshcTEd64CdQkUaEmCBgEIDCDdglTEcrVr5xHysOKBqKKgkB2\/dmly\/uTzzOUxLRIQYIgIBAAAr3MqDQ5gab8v3K4l41vOZziU0WWk6bMJIMPkwvIcooi5wcVlSP6a0Hp2gSixB1rUsZZem\/Ek8p3Uc9fOQRQs5GKdwkp6kiv4dLEa3+8WfH2oE1aOVLyvdbn595hdMUdSx8D7cVCk76frCUzpOJLDyl9oj7Hi4p2acOv6\/++WjL25gc7FZ+6tr2hEt3POHfuSweXPT\/h+y\/cU4xcSz0wdBWlzTpdC++h4vv4ep933wfuu2fIfpvoxRu2tCORdhHSPkefP1c29MZZuciWsF\/mXyMh9rWDJyCozz+eSNHECgjASjcMraKJhOd8eUPKcDhquSNZYhnOj5RMrXhq9+9q1oVgbQgAAJeCEDhesGYfyZxDldYHeXPPk0Jtg+5L+jvVXcv7o++OgUP4zQ0cS8I1IsAFO7l9pQ2zKV0fR8qrotIeUjpehl6dCteOUVHRm74V3X30iWJvVQqJ2t6aBbd8vYhJ+3FnvjgOkXHd06OX2yfkyXlSgEpfJQhtYePMnzkATmneyJYgIVtXILCrfgEile+VA3d9EwD\/oKBvug4CZkxcfkhYNuL1Sc8WMX64YxcQKCOBBqlcEOEdixDKEVa\/T77yzPqlbcuRn2WFO7ST\/Sov\/zCLISKXNRdqMcte3+hfvbWxTZT8iqmD7g\/+c072qtY+qUM7c+r1bzkifMCjau7uXrG8wgNqffNJvWPRinckKEdbbMzycwkpXdj9jO\/JEN5kLLo6emJohbZVmSSHFJ6N3KavEKUYcpJrOgiT2LdRMyTlharziAUIeQMUUZV26ws71kZ+i9YHGtPeMvKonYKd3JyUm3evFkNDQ2pwcHBNnepotKgJqX7GLAkGX2UQUrlB88cVjP7+9vHU1ihsHPPuYkJ9bXfXxTr4BNCTh+8k+Rk0\/DPf3NSvfbbnvb+q85i\/kCvuuXjl9R3vvqZWOtW3nK6tLnLPZBzugl9sJDykNJJGuk9cslDuidrug85JRl89d+sPEPIKckYO9B0kXDV1NTUVBfPOT8yPj6u1q1bpw4dOqT27NmTSuE6F5LjjSEbw6wGmaFPjE92KB5TEdP\/aZX3x1+8TVXhTPPCT9+pfvKTn0QrVrrMVaupXNnJKccmtmZdZLunqSvkTENLvhc8ZUZp7qgCz5Ay5qpwaWW7c+dOtWrVKrVx40a1adMmKNw0vdW4l4+3RD+PTrQU1mXFxbfS3jCtimklSA5arMDIUSvExavUkxOtaE08YaCyT1z+m01WMgvP7+8tzRGdkC9hlnaBnFnoXfkseDaPZ8g2z1XhctPxKtemcP02L3K7+Kl\/qz6ccb36cMasCMalWbfEQrn6wrvttKsvnE0Fj8qgi8sxH+a8KV\/+fcbP\/yZVGbgZBEAABEIQCGUdLFThhgCJMloE2PmIftLqU\/9b6\/fJ1KhoBc3Hm2h1SheCSKTGiAdAAAQaQsCrwj1y5Ii6\/\/771enTp9XKlSvV9u3bVV9fn4pb4TaEMaoJAiAAAiAAAsqrwo3jCYWLngYCIAACINB0AlC4Te8BqD8IgAAIgEAQAkEUrlmTsbExtXbt2ujPw8PDavXq1UEqm1QInxXev3+\/Wrp0qdq9e7caGBiwPkLyj46Otk3mIYV3kVPnO3fuXPXEE0+oJUuSYzn7roOLnPoWhMTct3ycn4ucfG+RlhoXOXWeJPP69eujkwGhLhcZSZaRkRG1a9euSCzzqGAIWSU59aOMLE8R75EkJ8mmy1qEjCSDi5xleNfj+lZcnIg8+mJwhUsdZMOGDWrr1q1RfbZt26Z27NgRq9zyqLQtz3379qnjx49HAxQNCAsXLrROBOi+LVu2dOxRh5KRypHk1PmSkqX79+7dmziByEN+SU6zkxNzukIqCBeeOhtWFEUoCYknyVnkRNCVpV4PGoSffPJJ9dBDD0W+HqEuF5a6LHQ\/XaEXBi5y6u8N3f\/iiy8GXwhIcvK7vmzZsohh0f1Ub9ukOBF59MfgCpdgUyehFSS9ZLYIVHlUNClP1w5Bsr\/55ptRVkV0bFc59brSoBZ6UtONnEmTnLz6Qxo5qe2feuop9d57710RMS0v+ThfVzn1gS9vmcz8XWSkex599FF13333Bbe4pGXJ9xfx\/lDZLjzNSU4R7e8ip7kIKIqprc8mxYnI4x0qROGyOZYqRAqXZz55VNAlT3O1pU8KbGblomaSaeU0X0gXFj7uSSMnm5pmz54dfBXuKicNGI888oh68MEH1WOPPVaYwuXQqLb+qZv1qA1DmxddWLLCJfnISsA\/9XCvPvqfy+Q6iaX+PE0E77nnno6APXnLqCtcFznZ6lbE9pxru+vjPPVfCoRUxFaXre1CbhVB4WqzSZfOzUqsyBWuq5xFmW5cXkKz40uTnDwGOVc5edClfeYiLDKucuqMQvN0kZEHtjVr1rRNi2ztivOX8N3uLnJymTzRevjhh4NvebnIaa4ui1gIuMhJPM0jo+fPn4+2FUP7ljRS4VbVpMyNVUTHTmNmKmply3xczExmx6cXkleQoV5CFzltDjShV2YuchbN00VGc3Aua5szS5q0PP\/888H9ClzfdZuptozvkK1vht7mSpq41XqFW3WnqSJXuKYijdv3LMrJQ+\/ULo4U+n5eUZMYSU69TqbC8L36SspPktPcHy2CpySj2X9Dr8L1CbOLg2RR5mRXOW2TnNAWA5cxyXxvinKQjHu\/aq1wqdL6sZUiPD5t4PU9MD1Klk15FTGYmatHOr5kk\/O2225rR\/viZ4o4cuPCswxHBVzkNNmbn5kMoXhd5NR5ht7D1VdlcX2TPFT1ehQhYxo5y+LglcSzbMeC4sbOuCiEId4dqYzaK1wJANJBAARAAARAoG4EgjtN1Q0g6gMCIAACIAACLgSgcF0o4R4QAAEQAAEQyEgACjcjQDwOAiAAAiAAAi4EoHBdKOEeEAABEAABEMhIAAo3I0A8DgIgAAIgAAIuBKBwXSjhHhAAARAAARDISAAKNyNAPA4CIAACIAACLgSgcF0o4R4QAAEQAAEQyEggiMItMiReRj54HARAAARAAAS8EMhd4Yb+wK8XKsgEBEAABEAABDwTyFXh0so29Ad+PfNBdiAAAiAAAiDghUCuCpclDBkc2kZl5Onj7T+fGJ90BrdgoK\/j3gUDvWrojhucn8eNIAACIAACIMAEClW4N910U+qWuDTrFvXhjFnRc+bvcZldfeHdjqSrL5x1LvfDGddH93KZtgc5f8q3593Xolvob\/SP\/+9cIG4EARAAARAIRuC5555TixYtClJe4Qr3jTfeuKKiJ8YvqoOvn1PvvPOOeu23Perk+EV14PVzV6w2Z\/cp9cm5M6O\/Pz50qxXYsWPHEmFK6TQpMGUk+eiinycnLqr\/9foZNbO\/X9Hq2Sbr3I\/1qJ6eHnX3zTMVrZrp5\/LFLbn5kuSQ0m1ymkCkPLKmU3lSHnWRU6qnCwsfeUg8fZThI48yyOlSj6bI6cJCukdKp3cgK0+XMqR7pHRJRp+auFCFu\/DTd6q\/+9EzkdI6cHTiCmVFJtwF\/b1q\/kCvWr64X83vp5+disonDFtevhqDzNoHj050TByofnQN3TEn+mlTxK718yWna3nd3gc5uyVnfw48wdMvAb+5VaF\/hpSxUIU78Oc\/jlrXVKxl2ifNozFoZcyreKq\/TRGTEt70pYXOvT8POZ0LT3Ej5EwBy+FW8HSAlOIW8EwBy+HWKvAMKWMQhRvXLlJFJVOAlO7DrCfJ6KMMzuOaj89Roy+\/HeHSlTBNSL6yuDcyW9NkhFfGOtcQcvrgXRc5fbDwkYfE00cZPvIog5wu9WiKnC4spHukdBqfsvJ0KUO6R0qXZHSYVzjfUrjC1TesCQxdvIH9wgsvqHnz5rX\/nzad7j916pRasWJFlK\/5vI90Jk0y2\/J3Sdef1et\/8NARdea9S+r1\/9unnv3lGfXKW629Y1K4ZH7+2FUX1bc+PzPiw2Xrz6flVQbeLry4jr556\/lJPKV0kjEkT5OF2d\/j0uP6jj6RtNXVnGiadfWdbr6rZv4+0vW+Z8vfJZ156j9tbZElPW4sS3r3494VHmtt\/YPaNG7s5P6dlG57n3TnJB\/pcf2XLIg\/GjsS+QDxIuaVBxZ2+POwIm6UwrU5TTlPFxp6o7kf3DquNCfTHnBDUaLaIAACNSDAW3TkC0QOt\/R\/XphI\/j9QuDXoAKGqQB2LZnC6CXr5zTPV3Yv7U+0Bh5IX5YAACICADwKsZEdfOhM5o7IvEI19aRxQoXAvt4Zke5fSTbOQrZGlPKR0H2X4yIPlpNUvHU9iMwqvfskBS6pL1nQf9fCRh1SPqpQBOTvfWKlds6b74O0jD6keVSkjDznNBQaNb0s\/0aO+\/Nn5sUGJJJ5QuFC4V8wNpE5jS+fOOfJ0a2+czgN\/4675sbO\/bsowBS1DHpIMeQwEZZ3MhWBRFZ5gMd1Lq8KC\/FimPjpbjb585opFBK9ipbpI6VC4PuwNyKODACtf6ri8v4F9X3QSEACBshHgY5NJStanzI1SuEleymk963B\/p5d3HA8+fvTDn55Up89fivY+SPkuv+Ffo1VwnOci+LrxLYKfzctUn9kjvbPtzFU5+BTPh8aiPb+61F7J0lj00NeWRBa53\/32TKRjTS\/nJC9r1\/7fKIWb5KUsmQKk9KqYuoqWk\/Z9zZWv7bxvCN5Fs+CZs1RXKd1HPXzkATmrZ0YN0WYhynDtv3\/\/m6tix58QckLhYg\/Xyx5u2v1V2i858PZHFO\/56g5XLi+Pyz0hXqC6lOGDZwgWVZETLMox+bA5PsVtb4VoMyhcnwZ65NUVAdPhatOXFl02PePzhF0BxUMg0GACtiM8ZfEhgcJtcMcsY9X1o0bmqreM8kImEACBchDQJ+4ciCLyFwn8EZokGo1SuAjt2HIEYKcN3TFAd+TImk75Zg01+F+f\/on62cR16p+Ptj6+QE4NKz91rRpe87moP5OMWUNp8ouRJVQmM0xy8srKU3reB+80POOcfuJY6GZgW13SpNvqapqZs4Z+NFmY+ftI1\/ueLX+XdIR2XBSNDRQHgMPRkqL9b\/++ZRkLFdoxTf9tlMKF01TrNZb2KrKm+y7DnLnSrJUcrcibMOljzlI9fMsZN7OV5Mia7qMePvKQ6uGjDB95hJAzRBlNZnHg6Ln2eVm2hN380Um1arn9W+U+WPnIAwo3bpTE30tHIO58b5pPC5auUhDzlXsfAAAgAElEQVQIBEDAiYD+\/tMDrWAUra+aVeWCwq1KS0HODgI2R6u4zwkCHQiAQHUJ2FazVX3XoXAv90PJBCSl+zA3hCijKnKmYaF\/0ShNPOc6sujWrF0VFlWRM03\/rXubdcPCXM1KMYyr0i+gcKFwr3jfpRdESi+q85t7vV9Z3Kv+9Iu3RkeMuh3UpLpK6UWxMOsLOaeJgEV5WdA7TGfzyRGK3ttXv3tXJGxd2qxRChehHac995K8armD088iQgf6KH\/L3l+o\/b9+vx1OkpTvtz4\/s7L1KVN7xXkp+wh9Z2t7c8BF+cWHRtTHBh\/ts\/Of\/rd69til6NN3t9\/Yq1beeq36zlc\/054Z6Aq3yu3fKIWLD9BXdx+nW8ltHs5wsuqWJp4DAX8E4k4fJFmk\/JVeTE61UrhjY2Nq7dq1Ecnh4WG1evXqNlWpopLJQkp3MXtIeUjpPsrwkUcIOX2XYXo4L795phq6c44anDWZ6WiRbzltw0CIMqrSL6oiJ9rMbra2RYGiyHJZ38Oq9AtJD\/mcBlw1NTU15TNDPa\/x8XG1YcMGtXXr1ujP27ZtUzt27FADAwPR\/6WKSi+IlF6VBq+KnHny1veJ+Lu9cateSQ4pvSq8IWfnyCS1a9Z0H7x95CHVw1cZ\/NUw+nAJXXSkR48CJckhpfuSM+lcv48yJD3kUz\/mqnBpdTsyMqJ2796t+vr61ObNm9XQ0JAaHBx0Urg+K4q8qkEAR4uq0U6QsroETCcoUrJN3tKplcIdHR1V27dvj3onKdxly5a1zcohK1rd16O5kts+G9jkgaG5PQE1z0rANBvTmdmqBajIyiDu+ZB6KPcVrqRwdQjksUwXmxCyxv4lk0fW2L7S8yx\/lti\/VYmlLLWHD942nmz6oqMJevzmOC9h7kNJXsScpvc3\/f6s6ZSvxEtKT8Mzzks0joVuirPVNU26ra6mqQ+xlKe9mNkUG+flK6Uz27i+q6eb78437prfDsGa9Dy36YoVK6JX0ta\/qE2T0m3vUxliKZOCNa9Qzru5K1yYlPOalzUvX9PJihw7qhrdpnmthxqHIhD3KTxYh+wtUJsVblanqVAdFOVUjwD2oarXZpA4XwK2uMZl+xRevgS6y702Cpeqrx8L2rNnT9thitJCVrS7psBTZSeAM71lbyHIlzcB2+QTlh936iH1UK4mZanKUkUlt3MpnfcesnwuLkQZVZGz7CxMJys2Odv6oVSXrOk+2tRHHlI9fJThI48QcoYoIxSLg4eOqBMfXKdGXzoTRYKiM+z\/8O3Wd6l9yFCWPEK0maSHJD2WJr1whYvQjtUN1ZjklMQvLP0MHYrS5ihC+1dlldcHnzinKYR2XNRWQHpfNBVKlfjRxJIdCMkJis7P3njNRFRP0ykJ7S+3f6MUbijvsDSzENxbDwK2b\/Vy4PV61BC1aAIBOEHl28pQuCUyi4Qwafgw34SQM0QZebEw97nowwnDa6ZNcOYrLdVVSvdRDx95QM7plq0aC9M\/QY8EJdUla7qPvucjD6kePsqAwoXCvWJKJ3U8Kd1HxwxRRt5y8iD2w5+ebH+1yLbXK9VVSvdRDx95QM5qKVzqnz8aO9L+So\/+LWl9UJDaNWu6j77nIw+pHj7KgMLN14KA3EEgIgDvTnSEshBIWs2WRca6ygGFW9eWRb1KScC219v0+LKlbKiaCYW92XI0aKMULryUw3vx+vCKZVMO\/axTfqfPX1IH3v5I2wt05aeube\/1ltnLuUpetjzM6ubCpshPHvStCd4ZNfry24pMxrQ3e+98FX3kHV7GLa9i01ScZ\/9olMJN8lKW7PdSug\/7fogyqiJnk1jQOUdSvDQw0gBp7qU1iUWIutaljLh3WbeiXLp0Sd00+9qOc7Np9merMl5URU4o3HJYGiAFCLT3emk1Qmcf6YpzZAEuENAJxJmMEQWqXP0ECrdc7QFpQKBNwBZKsqWAbwAlEIisIQdfP9eOAMUmY8Q0Lm\/ngMK12PBtzVV3M1QaMxNYTNMKyeLrj\/8iCq3Hq17aj1u+eOYVe1Dov8c69idNHiHbLGnol+SISycLyIGjE9G+LH1C8p5bZsV+b7bbMnyOB5IMPszBPvIIIWejFG6S05T0vVApnRpL+p5t1nR+CfA93NY3M8vAk9okqT18f3+YQ0k++8sz6pW3pvd7P3bhhPryZ+fHOpX57L9xTiVxLPTBkJ\/le7lP82AnpdP9Wb93Kz1v9i1zMPeRrr\/Ltvz1dD4v+9pve9pKds51PepP\/s0iNThrst3\/bA6FOtdu0lk2vb18O\/Rxm9bxe7hm\/26UwkVox\/KaWiBZegJscj54dKJj5Quzc3qWZXxCX8lSmy7o71VDd87BlkIZG8tRJihcR1C4DQTKTCBO+ZLM+Bh4mVtuWra4Pdnli\/uhZKvRhKKUULgiItwAAtUiYCpfkt6271utWtVTWttEqbU3DyVbxxaHwr3cqtKGuZRu7sPYOouUh5TuowwfeYSQM0QZTWBBAzr9+8efHVO7XjoXdUtWvubqV2KeNd0Hbx95SPXIswxexbLTE7eH\/rEAHjuKlFMfvyQ5sqb74O0jD6kePsqAwq3jNAp1AoEYAlj9FtM1vj16WJ0cv9jea8cqtph2KLrURilchHasV2hEennqFOoxdH0otOTUR2dHZznZ65lXXPRJQbr4s4K6Z2qclzI+QL4osibQF3jIg3zig562gqWIT\/iAu\/yBdl5F6u+CubKscv9rlMJFaMfW\/E4ynWRNr0oZVZFTag8f9dDzGHn6uGLPZ93kuXGwM\/6uuVoILWfcakWSQ0p35cmximnCYvKa3afUJ+fOVElBKCQ5pHRXOfWYyUW0WYh6VIVF7RTu5OSk2rx5sxoaGlKDg4Pt\/iVVVOoUUrqPBpdk9FGGjzxCyOmDd13k9MGimzxotUYXHU\/RFYquhBcM9EUB8SkARzdl5KEApHZPKydzIMVKF+2\/Pv\/au9E3jnUWjw\/d2q6OSxm+5bRNQCQ5pHTKM6ucLmVI90jpVZFTYhk3iezm71dNTU1NdfOg6zPj4+Nq3bp16tChQ2rPnj2pFK5rGXneF7IxstQDcmahd+WzVeL5P3\/2q3ZIQZsSptqRAiZFXNR5YB88+Qws77tyq3H4RH2i0W1v8CFnt2WneQ5ypqGVfG9IlrkqXFrZ7ty5U61atUpt3LhRbdq0CQrXXz\/pyClkp8lSBciZhV66icGBo+eivWC6TEXMq8DoZ3+vmk9BHAb6cjsfbLY7r1Dp58mJ1qq95b09GTky0UUhM\/WLA02QrHREZ34\/\/WyF0fR1oX\/6ItnKpwo8Q8qYq8LlpuNVrk3h+m1e5AYCICAR+HDGLPXhjOvVpVmfim6l3\/lv9NN2XX3h3Y4\/X33hrFRMO53y50vKn\/LlsnrefS36Xf+bc6G4EQRSEAgV8bBQhZuCB24FARAokAA5bekXrUTTXGzO5mfwdaU09HBvXQh4VbhHjhxR999\/vzp9+rRauXKl2r59u+rr61NxK9y6QEQ9QAAEQAAEQEAi4FXhxhUGhSs1A9JBAARAAATqTgAKt+4tjPqBAAiAAAiUgkAQhWvWdGxsTK1duzb68\/DwsFq9enXhMPis8P79+9XSpUvV7t271cDAgFUukn90dLRtMg8pvIucOt+5c+eqJ554Qi1ZsiSkmMpFTn0LQmKel\/AucnLZRVpqXOTUeZLM69evj04GhLpcZCRZRkZG1K5duyKxzKOCIWSV5NSPMrI8RbxHkpwkmy5rETKSDC5yluFdj+tbcXEi8uiLwRUudZANGzaorVu3RvXZtm2b2rFjR6xyy6PStjz37dunjh8\/Hg1QNCAsXLjQOhGg+7Zs2dKxRx1KRipHklPnS0qW7t+7d2\/iBCIP+SU5zU5OzOkKqSBceOpsWFEUoSQkniRnkRNBV5Z6PWgQfvLJJ9VDDz0U+XqEulxY6rLQ\/XSFXhi4yKm\/N3T\/iy++GHwhIMnJ7\/qyZcsihkX3U71tk+JE5NEfgytcgk2dhFaQ9JLZIlDlUdGkPF07BMn+5ptvRlkV0bFd5dTrSoNa6ElNN3ImTXLy6g9p5KS2f+qpp9R77713RcS0vOTjfF3l1Ae+vGUy83eRke559NFH1X333Rfc4pKWJd9fxPtDZbvwNCc5RbS\/i5zmIqAoprY+mxQnIo93qBCFy+ZYqhApXJ755FFBlzzN1ZY+KbCZlYuaSaaV03whXVj4uCeNnGxqmj17dvBVuKucNGA88sgj6sEHH1SPPfZYYQqXQ6Pa+qdu1qM2DG1edGHJCpfkIysB\/9TDvfrofy6T6ySW+vM0Ebznnns6AvbkLaOucF3kZKtbEdtzru2uj\/PUfykQUhFbXba2C7lVBIWrzSZdOjcrsSJXuK5yFmW6cXkJzY4vTXLyGORc5eRBl\/aZi7DIuMqpMwrN00VGHtjWrFnTNi2ytSvOX8J3u7vIyWXyROvhhx8OvuXlIqe5uixiIeAiJ\/E0j4yeP38+2lYM7VvSSIVbVZMyN1YRHTuNmamolS3zcTEzmR2fXkheQYZ6CV3ktDnQhF6ZuchZNE8XGc3Buaxtzixp0vL8888H9ytwfddtptoyvkO2vhl6mytp4lbrFW7VnaaKXOGaijRu37MoJw+9U7s4Uuj7eUVNYiQ59TqZCsP36ispP0lOc3+0CJ6SjGb\/Db0K1yfMLg6SRZmTXeW0TXJCWwxcxiTzvSnKQTLu\/aq1wqVK68dWivD4tIHX98D0KFk25VXEYGauHun4kk3O2267rR3ti58p4siNC88yHBVwkdNkb35mMoTidZFT5xl6D1dflcX1TfJQ1etRhIxp5CyLg1cSz7IdC4obO+OiEIZ4d6Qyaq9wJQBIBwEQAAEQAIG6EQjuNFU3gKgPCIAACIAACLgQgMJ1oYR7QAAEQAAEQCAjASjcjADxOAiAAAiAAAi4EIDCdaGEe0AABEAABEAgIwEo3IwA8TgIgAAIgAAIuBCAwnWhhHtAAARAAARAICMBKNyMAPE4CIAACIAACLgQgMJ1oYR7QAAEQAAEQCAjASjcjADxOAiAAAiAAAi4EAiicIuMQesCAfeAAAiAAAiAQN4Ecle4eqzPssRNzhsq8gcBEAABEAABk0CuCpdWtjt37lSrVq2KPji8adOm4B9yRpODAAiAAAiAQBkI5KpwuYIhv8ZQBqiQoXoEDhw9p05OXGwLfmK89fuJ8UnnyiwY6Gvfu+lLC52fw40gAALNIFCowr3pppuaQRm1DE7gwxmzojI\/nHG9ot8vzbrl8v9ntf8WJ9TVF95tJ1194ayz7FQWX1y+7WHKn\/Lln3RPz7u\/Vj3vvuZcFm4EARDwQ+C5555TixYt8pOZkEvhCveNN96IFfHYsWOJIKR0yli6R0qnSUGSjD7K8JFHCDklVi718C0nr0QPvn5OHTg6oU6OX1RvvPO+On3+UrtfLRjojX5f0N+r5g\/0KlqJ9v3uPTV79mw1v7+VtnzxzI5+KMmZhcXI08ejsp795Rn1ziStoqdX1pGcJGN\/r7p7cb86NzGhhtd8LvYdyVNOLjRLXTmPMsjpUo+myOnCQrpHSqe2z8rTpQzpHildktGnJi61wvVZ0W7zCtkY3cro0rGz5O3z2Sw8STHRPzL9knIlJcvKipXq3TfPVI8P3ZpZ5Cxydls41YXrRKZsmjwceP1cx8SBFTHVkyYJRcjZTf0gZzfU4p8BT388Q7KEwhXaLWRjZOlCdZVz9OW32ytXVj6kXEnhtFam\/dEq1VyhZmFZtgkMTzRaynhSERO+yCy97DOLo9UwK+Gsdc\/j+br2zzxYueQJni6U3O4JyTKIwo2rtlRRyRQgpVO50j1SuiSjjzJ85BFCTomVSz3i5NSVCplaX3mrZWbVzausUCQ5pPQscnJf9lGGjzzmfvXP1L\/7xrc6lDAxG7pjjiLHLR9l+MhD6p8+ypDykNJdJloueUj3ZE33Iackg8s74pJH1nZ3KUO6R0qXZHRT2253Fa5w9Q1rAkMXb2C\/8MILat68ee3\/p02n+0+dOqVWrFgR5Ws+7yOdMZPMtvxd0vVn9frr+fHv3abTc1l5Ss+n4UkK9qVfHVOvnLqo\/uWsaptO536sR91+Y6\/68mfnq4\/8v4nod+4PJoO6807Dk+6lfesTH1wXmdkPHp1oMyUF\/JXFver2eb1q1fJpczsPRLa+pQ+4Urqtb5kDttl30qabLMznfaTr76otf5d0fpf1n+bYo3M3+7bL83FjWdzY0M393KZxY6dLum380p2TfKTHjZ1p+m+jFK7kkOQ2b8BdVSDAe5S0\/8pmUTYP+9h3rQKDkDLqVoPRl8907HfTCnjojhsiCwIuEGgyASjcJrd+jequD\/i84mIFm9fea43wea+K7pQ18nTLmsTm5zLv\/3oHgQxBQCMAhXsZhmR7l9JNs5Ctl0l5SOk+yvCRRwg5Xco4eOhIZNY0V7Eh9xRd5JTuyZruo0195JFUD1LAZGkw98y5rfh9kVjkLaerHJKcUrqPevjII4ScIcqoCgsoXCjcK+YG0gsipefZ+U1TMTs6Dd055woP4iLl1KFKcmRN98HbRx5SPbiMaz4+p23mN1e\/f\/TJKTEwgFRO1vSQLJKCIEj1qIqcIepRFRZQuDB4lJ4AK9nRl85EzjnYiy19k6USkEJd0jEkU\/kiZGUqjLi5AgQapXCTvJR1LzaeLdFPm2cf0u1e2D55kbmYjus8e+xSpGR1j2JywEF7dXrZh+RhlmWuLrKk0+TqB88cVrteagXh4H1ffeWbJX\/bu+tTfuR\/5dgAvtNHRhulcMse2rHpphebuTgpmpPES0qvihkqRD3KyIL3ffWV76vfvUs87yvxktLLyCJu8SbVJWs6WHSSz8oTCvcyz6wg0TG775j6wMp7sl9Y1KO+89XPJBqJmtJmUj199D0feeQpp2l23vSlRbHRriQ5pPSys9BfCqkuWdPBovtxzTZ4QeFWwO5fRxFZyepHeEyP1TrWG3XKRgD9Jhs\/PF0sASjcYvk3rnQ6GqI7P5GSxbnMxnUDLxVm5cuBNvQwk14KQCYg4JlAoxQuQju2nMDY6YT6Uly4N+5n3aTTs3p4Pd0Rhj2MP937ThRSMc4pzWdoRzaL2erLf4tzOsqaXgRvW3198oxzWopjpZslbX0vTbrZt7ifkpPdgbc\/0vZ0Xn\/nTPW131\/U\/tCEblqVQj\/SvXqYVtOs6iPd9n7FhRW1la+zRmjHTketbvtniP7bKIULp6nWa551X0d6nsrgoBTmapaPekh5SOk+6uEjjxByhiijKixc5Nyy9xftmNm2VW8IniHKcGEhySGlV6WMqsgJhevZZND07Gg1S16lZDrm1SyZjX1\/0q7pnFF\/mYDp5UzHySjMJ\/3EBQJFEIDCLYJ6zcq0BaYgL1IMbDVr6ApXx5wIwkGvwo1ZYdGhcC83nmRakdKrYtLwKad5nEd3gJJ4ZU33UQ8feUj1qEoZTZHTDKwRt+qV2jVrug\/ePvKQ6lGVMqoiJxQuFO4V88Wkl5AGrB+NHWlHgCKzsW01K73IWdN9vGA+8pDqUZUymirnt0cPt7c\/9FWv1K5Z033w9pGHVI+qlFEVORulcBHasftQleQE9dTh96OQe6Rkl36iR\/3lF2Yh9OWiRdGEJc7LuY6hQeO8QOPqag6GZXveXPWSh\/O3Pj+z4yMKumIqm\/xl5wv5ENqxwrsAYUW3nZtFUPmwbYDSwhFAfw\/Huokl1WqFOzY2ptauXRu14\/DwsFq9enW7TaWKSqYVKb0qJg0XOflMIwUUoIsCU+iexk1iIdVVSnfhLeUhpfsow0cedZJTP9fLR4tov5d+98GqLHmEaLMQZfjgGUJOSQ\/5nIRcNTU1NeUzQz2v8fFxtWHDBrV169boz9u2bVM7duxQAwMD0f+likqwpfSqNHiSnK6enE1gwX1LqquUXod+USYWoXnaolmRz8LgrMnE7\/aG6BehWcSN3VJdpXQf9fCRRwg5JT3kUz\/mqnBpdTsyMqJ2796t+vr61ObNm9XQ0JAaHBx0Urg+K1qlvMwjPctpNXvnHBzpqVIjQtYgBGxe+dheCYK+NoXUSuGOjo6q7du3R41DCnfZsmVts3LIilahdyQd6amC\/JARBIoiYAbUYC99NjcXJRfKLT+BkHoo9xWupHD15iCPZbrYs1KKNSulkzlCj79q81rNms7y6\/F54+KvxnnNjr3bp\/7HqyfV\/sPvRx91\/8Zd8xXN0vX7+XedT5p0ek7ilTXdB28fPJlRkpdyVp7S8z54p+EZ56Ubx0I399nqkibdVlf9eR\/pJgszfz195Onjyvx4wh99ckqMxaz3PVv+Lul6DGXdHGq+q1liLbNscWOBj3RusxUrVkTVtvUvGi+S0m3jE4+NZp4h+y8pWPNKCjHsc8qQu8KFSdneXPo+FN1hOkH5bGTkBQJNJABzcxNbPX2da7PCzeo0lR5d+Z+wDQK6p2X5awAJQaBaBEwnK\/hEVKv98pa2NgqXQOnHgvbs2dN2mKK0kBXNu9GS8kdc4yLpo2wQmCbg6vUPZs0hEFIP5WpSlppMqqjkEi6l8z6Bvm9gyiTlIaUnlWGuZpPMxlI5WdOLZsHcpXpURc4Q9agKi6rIabYZ7fXSV7ToYier3\/32TOFHi0L0rRBlVKVfSHpI0mNp0gtXuHUL7Xj6\/CV14oPrFH9zlpygfvmXdseDJKce7qz0s46hCFE\/v6En45xOqhra0Rys86wfKd4f\/vSkoneXvJpf\/e5d0Rga5\/RkmzzmKZ\/tXYF8\/tqnUQo3lHdYmllIt\/eydyQ9DyeobiniORAohgD8K4rhXnSpULiXW0Aye0jpeZs0eG+Wj\/TYws25yOByj1RXKb0qZVRFzhC8q8KiKnK6ttk1H58Tfa3IPFrER\/WStqjqxiLLdlxVWEDhllzhmkd66Cs937x3iVq+eKZ1sub6omfp3HUpoyovaQjeVWFRFTm7aTNz1fuVxb3qT794azt+s+2Fl8rJml4V3lWREwq3aBuDpXxezR44OtHxrU4c6SlhY0EkEPBMwBa\/Wf9Wr+fikF1AAlC4AWFLRSFAhUQI6SDQLAIIqFGv9m6Uwi2rlzLt4fzNT46pV966GJmPaDZL4eHogtdw8z7wzuaxsrY\/vGRbx3vM0IF5emnrn8wkJayHZTXNqWif8O3j6mXeKIWb5KUccq+DXpjWzPVM22RMnsb3zldq1fJbE6d0IeWME0SSoSr7KVWRMwTvqrCoipx5tpn58QSepNu+XCTJIaVXhXdV5ITCDWixME3GC\/p78Sm8gPxRFAjUjQC++lWtFoXCzbm9+IU4eHRCHXj9XNtkDAeonMEjexBoEAGMM9VobCjcy+0kmVakdN2kYYtnTCbjjYO9hYdyq4rpJQ3vupu+waKzhSUeWdOr8o7EyWnzcuaIVua7IrGqOgu9vlJdpXQfLBqlcJOcprJ+n5WcGp795Rn1yrnropUsOTXcfmNv+8wsNWYZvoerf0uXOpDN0YOdLrpNp+ey8pSeLwtPZmQ6quj\/z8pTet4H7zQ845xy4ljoA5WtLmnSbXU1B0Kz76RNN1mYz\/tIZ0Wgv39x37a2la+zTvreLTtbcRxn2u+l871\/eOu16u6lS9ohJZOcfpL6NssWN1a4pHOb1vF7uNzOzLdRCtd3aMe4lSx5GccFpqiG4QNSggAI1I1AGmerutW9LPWBwk3ZEra9EsQyTgkRt4MACBRKIC64Bo1lWCzk1zRQuAJbM+oT3U5mGSjZ\/DolcgYBEAhHwLaIoBMU\/\/Dtz4UToiElQeFebmi2sVPno4uCUeiexbP7lPrCp+dEirbbOMbSpryUzvshWeIg+8gjhJwhygCLzlFOYp413QdvH3lI9ahKGXnIaSpfKoODbMSNfRJPKd1HPXzkEUJOKFyloiAUB18\/pzh2sb6KXb64X9ERHlwgAAIg0CQCHKCHxkb+mhHVf\/nNM7H67bIjNErhkpcyfQ6LOtLPf3MyCqW4\/\/D7ETr2Kv7yZ+er+f296sZrJqK\/I7QiQivy7Bn9ofU+xHkp5xnaMMmL1lzdQL58QhvypwTPTUyoXS+di8ZH2l6jizyfv\/b7iyLrH\/jH82+Mwr1h1X9Sd977tejIjr6CfXyoFUpRMidI6T7yCFFGVeQEi+kpNFh0LickHlnTq\/KOFC0nm591BcxjK53UoIvM0LR4yfu7vkWz4B4q9b3aKdzJyUm1efNmNTQ0pAYHB9tv6qz\/8F\/U6ns+rchETCtYcx9WAiWl+2hwl8aQ5JDSqyKnj3qE4BlCTh9l+MhD4umjDB95lEFOl3rUVc6Rp4+3\/V94AOZ4z6yE8xh\/s\/J0aTPpHildkrFLS7X1saumpqZan8DJ6RofH1fr1q1Thw4dUnv27OlQuCEr2m31qiAj1Q1ydtvC9ufAEzz9EvCbW5b+yatgkoidUHUlTN7Qdy\/uV7YPL6StRRY505ZF9\/Me98mJi5H\/z0n6KM1E68M049+\/15plSBlzVbi0st25c6datWqV2rhxo9q0aRMUbje9yOGZkJ3GQZzYWyBnFnpXPgue4OmDgO6MZVPE9DfeGyaT9IKBvqhY+htZJ\/l3vodlyqt\/0okVUqh06UpVnzjQ7yxrazVvd7TNS0Zbu+SqcLlAXuXaFK6PzoI8QAAEQAAE\/BP4cMasKNNLs25R+u\/8t6QSr77wbjv56gtnnYX7cMb17Xu5TPNhzpvypd\/55++dOOhcjn6j74iHcUIUqnC7IoOHQAAEQAAESkmAVp50cewEFvLE+KSzvLx6pgd8mLWdCw5wo1eFe+TIEXX\/\/fer06dPq5UrV6rt27ervr4+FbfCDVA\/FAECIAACIAACpSDgVeHG1QgKtxRtDSFAAARAAAQKJACFWyB8FA0CIAACINAcAkEUrolzbGxMrV27Nvrz8PCwWr16deHE+azw\/v371dKlS9Xu3bvVwMCAVS6Sf3R0tG0yDym8i5w637lz56onnnhCLVmyJKIXkYcAAAYySURBVKSYykVOfQtCYp6X8C5yctlFWmpc5NR5kszr16+PTgaEulxkJFlGRkbUrl27IrHMo4IhZJXk1I8ysjxFvEeSnCSbLmsRMpIMLnKW4V2P61txcSLy6IvBFS51kA0bNqitW7dG9dm2bZvasWNHrHLLo9K2PPft26eOHz8eDVA0ICxcuNA6EaD7tmzZ0rFHHUpGKkeSU+dLSpbu37t3b+IEIg\/5JTnNTk7M6QqpIFx46mxYURShJCSeJGeRE0FXlno9aBB+8skn1UMPPRT5eoS6XFjqstD9dIVeGLjIqb83dP+LL74YfCEgycnv+rJlyyKGRfdTvW2T4kTk0R+DK1yCTZ2EVpD0ktkiUOVR0aQ8XTsEyf7mm29GWRXRsV3l1OtKg1roSU03ciZNcvLqD2nkpLZ\/6qmn1HvvvXdFxLS85ON8XeXUB768ZTLzd5GR7nn00UfVfffdF9zikpYl31\/E+0Nlu\/A0JzlFtL+LnOYioCimtj6bFCcij3eoEIXL5liqEClcnvnkUUGXPM3Vlj4psJmVi5pJppXTfCFdWPi4J42cbGqaPXt28FW4q5w0YDzyyCPqwQcfVI899lhhCpdDo9r6p27WozYMbV50YckKl+QjKwH\/1MO9+uh\/LpPrJJb68zQRvOeeezoC9uQto65wXeRkq1sR23Ou7a6P89R\/KRBSEVtdtrYLuVUEhavNJl06NyuxIle4rnIWZbpxeQnNji9NcvIY5Fzl5EGX9pmLsMi4yqkzCs3TRUYe2NasWdM2LbK1K85fwne7u8jJZfJE6+GHHw6+5eUip7m6LGIh4CIn8TSPjJ4\/fz7aVgztW9JIhVtVkzI3VhEdO42ZqaiVLfNxMTOZHZ9eSF5BhnoJXeS0OdCEXpm5yFk0TxcZzcG5rG3OLGnS8vzzzwf3K3B9122m2jK+Q7a+GXqbK2niVusVbtWdpopc4ZqKNG7fsygnD71TuzhS6Pt5RU1iJDn1OpkKw\/fqKyk\/SU5zf7QInpKMZv8NvQrXJ8wuDpJFmZNd5bRNckJbDFzGJPO9KcpBMu79qrXCpUrrx1aK8Pi0gdf3wPQoWTblVcRgZq4e6fiSTc7bbrutHe2LnyniyI0LzzIcFXCR02RvfmYyhOJ1kVPnGXoPV1+VxfVN8lDV61GEjGnkLIuDVxLPsh0Lihs746IQhnh3pDJqr3AlAEgHARAAARAAgboRCO40VTeAqA8IgAAIgAAIuBCAwnWhhHtAAARAAARAICMBKNyMAPE4CCQR0B3b8g7u4dvJp0gvXfQqEKgjASjcOrYq6lQaAnkrWa5oXsrRtxIvTcNAEBAogAAUbgHQUWQzCHAEIKotRQGioygUo5s8yekc4u23366+\/\/3vR1Ghvve970UhTw8dOtTx0QEXT27zSFCSt7L+4QA9MpHu7ap7tRcV77gZPQS1bBoBKNymtTjqG5SAzaTMR7ceeOCBKOIS3UNHPyjUHV0cvOD6669X69atiwIvUPhDuu\/06dNXBKfXlSI9r0fD0le++nG2U6dORcfHSNGbEbT087SmMg8KD4WBQM0IQOHWrEFRnXIRiFO4ekQgXcHpgWHOnj3b\/tAHhT6MC\/quK9W44BxmkASixLLxijvuq10wK5erT0Ga6hKAwq1u20HyChCIU7h6aLskhcvfjeaq2oJFmF+JsZmH+ctctJLWL\/pmLgXnT4pQBIVbgY4GEStBAAq3Es0EIatKIKvC5S9rJX0vNslhij9g8Rd\/8Rfqr\/7qr6xfOpI+lwaFW9XeB7nLRgAKt2wtAnlqRSCLwjX3cGklu3fv3is+Y6jv4ZLpeMOGDe0vscTt4dJ9tD9MX+75gz\/4g469Yn2vlxqj6BCHteoQqEyjCUDhNrr5Ufm8CejfKjW9lHnPNM6kTF9NcomPbDo26bHKzTjaupcymZPJIYuuOG9oeCnn3UOQf5MIQOE2qbVR19oSwDnc2jYtKlYjAlC4NWpMVKXZBHzvtealxJvdSqh9kwlA4Ta59VF3EAABEACBYASgcIOhRkEgAAIgAAJNJgCF2+TWR91BAARAAASCEYDCDYYaBYEACIAACDSZABRuk1sfdQcBEAABEAhGAAo3GGoUBAIgAAIg0GQC\/x9LT+M03kNJdgAAAABJRU5ErkJggg==","height":0,"width":0}}
%---
%[output:7912476c]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAdwAAAEfCAYAAADr33fvAAAAAXNSR0IArs4c6QAAIABJREFUeF7tfX+QVsWZbmdjZWcsKjWOJMUiQX4IibXmEkXjOMKyqYrLvaaGrGQ2w+CuFkGC6HVTlRkYHKgFtvgpjrfAGMUtpEwqjCQRN1C6ZWW9uRQGNd6oU0lVshhlCr3I3TUwyaXiVMrEe9\/D7Y\/+znSf7tOnu0+fb57vH2XO6T5vP0+\/\/XS\/\/etDH3zwwQcMPyAABIAAEAACQMArAh+C4HrFF5kDASAABIAAEEgQgOCiIgABIAAEgAAQCIAABDcAyPgEEAACQAAIAAEILuoAEAACQAAIAIEACEBwA4CMTwABIAAEgAAQgOCiDgABIAAEgAAQCIAABDcAyPgEEAACQAAIAAEILuoAEAACQAAIAIEACEBwA4CMTwCBvAgcOHCA3Xvvvayjo4Nt376dNTc3S7N477332Nq1a1l7ezvr6upi9O8tW7aw22+\/nc2aNStJQ3kNDw+zvr6+vGbgfSAABBwiAMF1CCayAgIuEDhz5gxbvnw5++hHP8qOHj3K9u\/fz9ra2qRZk5geO3asJsovvvgi6+3tZfv27asJblqUXdiIPIAAEMiPAAQ3P2ZIAQS8IsBF8\/7772c7duxIxFY2On399dfZmjVr2H333VcTV5ngkrH0d8pr7969rLW11av9yBwIAAE5AhBc1AwgEBECfDRKJlEo+dChQ+yJJ56QCiUJ6KlTp2qjWx6G5sVZuXJlTagxyo2IZJgybhGA4I5b6lHwGBGgUeuyZcvYPffck8zJpv\/NbeZh5yVLliTv8Z9qhEvPSZBV4h0jFrAJCDQaAhDcRmMU5ak0AmlRTI94+eIpLsQUdhbnd7MEN+tZpUGD8UCgIghAcCtCFMxsfAT4qDU9Z8tDxeLiqay52vSiKY6cSqQbH1mUEAjEgQAENw4eYAUQSBY2LV26VImEOCdLIvzggw\/WrUamhFmjWAguKhkQKBcBCG65+OPrQCBBQBU65vDQAikSU77KGCNcVBwgUD0EILjV4wwWNyACqsVRvKh89Ltt27a6xVSYw23AyoAiNSwCENyGpRYFqxICqhAxL0N6BEz\/psMx0quUs8LGWKVcpRoBWxsRAQhuI7KKMlUKAdViqXQh0oun0qdMiaHpw4cP1x0LiX24laoSMLZBEYDgNiixKFbjI5BnERROmmr8+oASxo8ABDd+jmAhEFAiIBvlpl\/G6BYVCAjEgQAENw4eYAUQsELARExxW5AVtEgEBJwjAMF1DikyBAJAAAgAASAwFgEILmoFEAACQAAIAIEACEBwA4CMTwABIAAEgAAQgOCiDgABIAAEgAAQCIAABDcAyPgEEAACQAAIAAEILuoAEAACQAAIAIEACEBwA4CMTwABIAAEgAAQgOCiDgABIAAEgAAQCIAABDcAyI3wCX7AQnd3N6ML0vGrHgLgMB7O0ncf81ug+LnaQ0NDdWdhi+\/zd+MpDSwxRQCCa4rUOH5PbAT2798Pwa1gXQCHcZGmOv2L7j2eNm0aW7RoEVu7di2jDu7s2bNZT08P6+\/vTwqxdetWNjAwwFpbW+MqFKzRIgDB1UI0vl+gUdHu3btZZ2cn6+3tZX19fRDcilUJcBgfYVxYu7q6asbxThH3MS7KCxYsYPT+3r17WXNzc02IEWmKj1edRRBcHUJ4niCQbgwAS\/UQAIdxcCZGG8iiOXPmJGJKPz6SnTVrFuMXUyxevJgdPHiQbd++PXmHRr7t7e1MFOs4SgYrdAiUKrgzZszQ2YfnCgSee+45Nn369GD4qBprcGhPwZtvvmmf2CIlOLQATZPEhR9yYV2zZg1bv359EjrOI7gmPjj6qS8mJWn65Q\/cg+Aox99d8xV28SuPOcrNPJuQfuhMcFWLAHixZdeIUUVRFfbEiROZgjLen2dhZ17VzN\/MaqzBoRxHXR0Fh4zpMIr9uQsO+V3FO3fuZKtXr65N25iGlLNs2PHsMNvx7Ak2tbWJTb2kiT3\/xgg788Dn6ips2RhftfEou+iii2r2ka2vrb+hZqNv+1xwaNqSOhPcrCvA+EXZc+fOTcIiNA9BPyqo2EMkYOnHR24x\/FtlC9kpVoS0rb6fh6wkVBYbwTWthOP1PXBYfeZtOCRf2rRpE9uwYUOy8InmZ+lHc7c2i6aybPjM5hfqxOv5X42wRd98lc2b2cIO3X116QSQffQTBZZ3Eg7ddTWbd0WLdxttOLQ1ypngyhYBkFG0YGPLli1s5syZ7LXXXhsjuCGH87YgxZguZCXRCS46TedD+3k7XeBQjpuso91oHV8xIsjncEl8xfndlStXJiJMP\/H99E4B1cBl7oPDrG\/hdPbl2R+MGcjwZ30Lp42ptz4GOjL+Fj30Kjt5dpQ9deukMYOsrz0zkjwjIfY9sAnph04EV7UIgCoQjXzpd\/nll7PBwUFjwfUdRqh6\/iErSVaHA9MC6nl0XR0Dhwgpu+jMy+oRidknWpvYQ91XKsP2rV\/\/USJof\/jNO6VM39Ho9htLrmSXffis9Pv0\/MaZLay3rcmrfSH90IngpisNn6+94447EsFdt24do43cMsFVjY54L0sVYg71PKbwtohFyEoCwZUjoBNU3XNwCMH1JbhiKFlVD0mU6bfr5havgib7PhfTrA4B2UadgkcXT2Kd865UQqXzM93zkH5oJbh8TvbUqVN1p6FwRPgigI6ODrZ58+Y6oOhvu3btSv4WsqAuKnZMecSCXSx2xMSNqS2xYBeLHaa4xfReDNilbSCR6r5uUjK61f3o3VBzpaIt6blllZ28U+Bzvjkkh1aCmwYnaxGAKMJ5Qsq6ijLen4esJLoRLuZwqzGHq+IRixcv7IiIfR5exqHIHy04+vYLb7HDt08xWnxKgvbe6GgyivQZzRPzFuePqTy6OVrqFPCV1Xn5Mck\/ZFvqRHCpUKpFADrBRUjZfOEIYckrbshKohNcbAuSIxRTKAsc2nFUBQ7FtkA2utWVQTfK1aXP83zw5dPJNqU8235uGniBNTc1KVdV5\/m+qsMSavGuM8G1GXFiwQ0W3GQd3lHUkWJPj04T5nBt2s10GrEeyUK1Oj8gQZs9uUUZgtalz\/Oc7Ou+7s8YrY7mP5P0NComkaY9uumfSfqsdiakH0YruC4qYiPnEbKS6EZHCCkjpCxGX0JtK5E12D5Cjlkhzxj8kIeUaSsNrUymVb15+RAFzSd\/dMiFGO42CfmSPaqymabP2lYWksPSBReNdfUb61DhmEbrQIV0dF2nCRza1a4YOOQ28G02NodFiKuG7ZDQpxK3KunfHvuGOJdrk16VJiSHpQsu5v\/s5pZCVhLbxrpoqKfR04NDhJRdCAfVo6t6nkyykq3mNfUjlaCZpleVhadXzRWb5k+dAjrIg1Zgiz\/T9BBcnKWs9DddJQrZWGddfo15+GrMw4NDu45tLH6o42\/kr\/cqt\/foyiAKIm0lci1olP\/9L46yt86MFuoQ3D34C\/bjN0bqFlylQ8oylnXlD9mWljrCndS5mS1cuJBNmDAhwencuXPJf2P6t8oWbm9Zz5988kl2+vvrXXSQM\/OgLV9Zl19jS0n8W0rAoT1HujnCEI21CX8X3\/Y4+\/nG+Ykv287BkijSz2YOmNJlbSuiOWIS87aJ7yXfEBcxZc2Rp\/GXHYSRJ71sjj8Eh7yRLVVwP37bI6zt+uu9i0YjfuDFl15i\/\/6tO70XjR9iorr8OmRl9V7YwB8IhR049EdsCA51\/FE72nPrf65b+WtbYn7co2w1sG2elM70oAvdN4rOA8vyD8Ghc8FVXc9Hlxrs2bMn+d7kyZPZvn37krse6YdwZPzhSOKVH1hCnKUvvwaH4BBbu\/zea6zzwYl3fo+9+8jfKLVKF1IVn8sWT+VJLzPC9baj9FxzUfsqKbiy6\/nopiBqoLu7u1lbW9sYLrLCkXzorwpVhHru8\/QVXRgm63moSqJzdqr8+NkjkL6b1D4ndUpw6APVC3n65lDHn25a5+jRo2zKFPXJU+nn5NM\/vWdaLeybN306bCsulpKFdCn\/+fPl4XBCOf284\/G32bovzKrNNedNL7aruoGf65rjLKQsu55PnHvgo1qxAKFEwzVoMeQXCjtdOCuUHTFg7tqGUNiBQ9fMXcgvBIeh+SOBlC2eskFRtdDJJi+ehvJULcCyyVd30pZNnqo0TgRXdT3f8ePH2dKlS2vfFu931PUsioYJGj19CEcnjkwWbGBrl9y9dHUQHGJbkEljXsQHKX9dPUw\/J0GjH7\/8IG96sUwkZhs\/P5H9\/c2fVhbVJv\/0qLnItIYuJG\/Ckek7TgQ3\/TF+Pd\/27dtZc3Nz8piHl9vb21lXV1fyN4SUs1f2EUaqkHaoxpps0F1+DcGNW3DBof08u04MQvmhrQ\/aCC6lERdP6TDIek5zwuIF8zJPscnf5L5f\/q2s\/Cmf\/\/n0t9ipZ\/6bqWYWes9KcE2v56OVrXQJPf9R2Jl+fX19NcHFCTd2\/IVydJ11uvkj2y0KLufOs451K7qloEh6cGh+cQc4VHuiDx+kedIFn5xYu8A+q\/Ov8nHxOEYf\/PHFU0V8kDoEv33sVhZKh6wEN0296nq+BQsWsCNHjiQCy8PO9P98AVUsDY5OVGJ8Hgt2sdgRI0c6m2LBLhY7dHjF+DwG7HzZUPQoRd9zo6qTp0zrCb9r9+cDX6qW4KZDVnPmzGF8dCtuC5LN4eIsZZylbOogjfaer4YyL06x2JHX7hjejwE7XzYUWTzlemGTjOuiC7J4h+C2v7qmeoJrU\/mxh9N+bsmXk+XlEfPw1ZiHz+IVHFabQ11I2XZbDz+OcfV\/OpdrWxGFeOlkqUN3Xc3oIoWi23ay0pNobpn\/Ybbqlr9IqrjptqPvHv8QG3z5neSYyJBtqZOQct5Gmr8fsqC2NsaaLhbsYrEjVp50Qhdq7qgKdoBDOwR8+qBtWNjVyVI6RNIrqnXv8+fi6N0nfml7ShdchJQRUjZ1kkZ7L6SjQ3D91J4YOPRpA5\/nlN1CpEKUxIwuQODbivwgfyFXm06B2CHwiV90gostJfLqqFsqH7KS2DbWujKM9+fgMP8e0XRdLLsOxcChzoaiGOkWT6XzT49ui35flz7v0ZHpDoEOP5cdhtJHuBDc6guuKkoR6vhNQlC27SCGLUkq2+jvIR1d12kCh+ZblEROY+DQ1xwu96m\/2Po\/2J\/+6Z+yH\/bckFSjtF+Jc8QkZh1XTkgOu+Dpfc7hkj2U\/xefej+ZM77sw2cTG8WDMMTv737mZ+zRn4zUblYK7YelCy5Cyggpu+xBVimvGBrr0A1OlfgxsTUGDkPYYBK23fHscG0hkgl2Lt8xDX1TOegS+76F02qfD4Ef\/1jpghvDohGXxIfKK2Ql0Y2OwKEd6+DQDreYUsXAYSgbdAuhimwjcsGp7vuqbUSh8KMyOhNc1fV8WadShSyoC0JjyiM0dqqbn3ThrBjCuj5OueF1ocgpN+DQPIzbaByq2kvxXPqOjg7Gj8cV39+2bVvteFweoQgRKbzlO6fZ1Eua2K6bW+rCtuQDdDLVjI9NYLS4SubzIfijEfaOZ0\/U3XTEQ+D0feowrP\/LFjb3sqa6kHNIP3QmuLLr+dK3BdEhGHT6lHjSFOaOzBudsuaOxEZg\/\/79dVctZlVW3WKH8f48pKODQ3l3uWgdtOVQ1l6ShfzWtUWLFtWuNp09ezbr6elh\/f39SSG2bt3KBgYGasfm6mwoWkYxvWwF8lUbj7KLLroo2dMq+7n8vkn+6ZH4vU+8yvb8ZGRMKJnnpcPP5UDLmeDKruejXhk\/2lFmdMiCugQthrxCYUcj2927d7POzk7W29ubHNMp3m0cyo4YMHdtQyjswKFr5i7kZ8uh6jrT5cuX13yMizINUuh9Or2PLoNJ3zFua4MtKiS6U1ubktHu82+MJP+vElvbbxRJR\/O53C6y8eTZ0Uz7QuLnRHCzrucbHBxkJ0+eZENDQ0wMkfBQSBFgx3vakHOnsrOwwWHxGggOi2NYdg55OVS1l1QOPpKl+8P5rWuLFy9mBw8eTMLL9CPBTd+6FhqD0U99kf3x4kvZxa88FvrTxt8jG5t++QOj9\/NyaJSp5CUngpvOl1eUa6+9NqkoYs9MrCi2RiNdeARUghveEnzRFgFwaIuc33S8vVyzZg1bv359EjrOI7h+rUPuLhGwElzT6\/lWrVrFnn766drEf\/p6PpcFQV5uEFBxi8baDb4hcgGHIVA2\/wZfcHj48GE2efJktm\/fvkRQ+Y+m3qht3LlzJ1u9enXukLK5JXizbASsBDdttOp6vhUrVtRCJFOmTBkTCim78Pi+OQIQXHOsYn0THMbBjKq9pPURNoum4igVrDBBwIng0ofEZevi9Xzi39PX85kYiHfiQACNdRw8FLECHBZBz21aVXspzu+K7aX4fnqngFvLkJtPBJwJbh4js\/aU5cnH57uqfac+v2mStxieovfL6sSAQxO25O+AQ3Ps4IdqrKrgg2Q9OLzAYXDBFffmkhnpPWXmrujvzaw9i\/6+apazuH+P27lkyZK6jfBmOdm\/BQ7tsaOU4NAMP\/ihGqcq+CBZDw7rOQwuuHyBgGpPmZkr+ntLt2fR35ftcpbt57PLyTwVODTHyuRNcDgWJfhhds2J3Qf5yDZrD7+Jb4R8J4QfliK4tDdXtacsJMBZ36rCfFf6JK9Q2JGzg0M3aIPDbBzhh3J8quKD4ig3fWiOGw9yk0soP4TgKviK3dHLnBepirODQ3VjBA7dNNRl+WFV+KuC4IbksBTBzTqmzI0bFM8l5sY6VG9MhWIVwlmxOzo4NPNR+KF6hFuFdhR+WM9fcMGt2mR\/bGGQ9B4+s2bL7VvgsBie4NAcv1gFt2wOq+KDMQtuGRwGF1wioAp7ymJ1dOrV7tmzp67FSl\/XZd6c2b8JDu2xA4fm2MEPs6cFli5dmrwQ895ccHiBw1IE19zd8CYQAAJAAAgAgcZAAILbGDyiFEAACAABIBA5AhDcyAmCeUAACAABINAYCEBwG4NHlAIIAAEgAAQiRwCCGzlBMA8IAAEgAAQaAwEIbmPwiFIAASAABIBA5AhAcCMnCOYBASAABIBAYyAAwW0MHlEKIAAEgAAQiBwBCG7kBME8IAAEgAAQaAwEILiNwSNKAQSAABAAApEjEERwQ97GEDneMA8IlIYA\/LA06Md8WDwalR7y41nFC9s7OjqSa0ybm5vrjsMt4yjXeJCrtiXeBVesQDGf91ltGmE9EMhGAH4YVw05cOAAGx4eZnQ5ivjjl6AvWrSIrV27lnV3d7PZs2eznp4e1t\/fn7y6detWNjAwwFpbW+MqFKzRIuBVcKlHvXv3btbZ2cl6e3uTytXW1qY1Ci8AASDgDgH4oTssXeXEhbWrq6uWZfqQfy7KCxYsYFW5is8VPo2aj1fB5aDFeltEo5KKcgEBGQLwwzjqhRhtIIvmzJnD9u7dmxjHR7KzZs1iJLjHjh1jixcvZgcPHkzCy\/SjkW97ezsTxTqOksEKHQKlCu6MGTN09uG5AoHnnnuOTZ8+vXR8wKE9BaE5VAkuOCyXQy6sa9asYevXr09Cx3kE1xd\/v596I3t\/4ifZxa88Zg+Q55S\/u+Yr7E9+92vW9MsfWH8ppB+WLrhvvvmmFKgTJ05kCsp4f05OpsLOtuapFnLw\/HjDwBdy0N+z7BjvHOnK74PDLO6zBBd+KEcuBIfkdxQy3rlzJ1u9enVt6s00pKyrR7oyyJ5\/ZvMLCSBTL2liz78xwqa2NrHX1t8QTVu949lhtuPZE4ldH2tm7Kf\/a5SdeeBzVvbp8LNtT2XpohVcl4VsxLx8VBLVQg7C7\/XXX2fLli1jc+fOra2c5IIr9hDJeenHR98x\/FtlC9kpNjZpW30\/98Gha8FtRN9xWSYbDqnjs2nTJrZhw4Zk4ROJLf1ojYvNoikbG7Iw4GIrCuzdg79ggy+fZofuuprNu6LFJYRWeZGNon0nz4wy+tu8mS3s0N1X58rTNX5ZHy9dcNFYnw8L523sfVQS2UIOso0W3WzZsoXNnDmTvfbaa2ME1\/VIO5e3VPhlHxzaCi78MKwfitEkPodL4ivO765cubK2ill8P73bg+qRK\/4WPfQqe\/M\/zrE9t0xiN86ZVdc2fe2ZEXby7Ch76tZJuTvVLju9cx8cZn0Lp7O+hdNq9vH8W7\/+I7bysy1s25ILoqvrVIf0wyCCq2oEEI5Uz8HqwkCuK4lqIQc1AjTypd\/ll1\/OBgcHjQVXV4bx\/tw1h7Z9D\/hhPH5ow6GuHuXxMxolfmPJlXWjWDE9PacwsziKzJO\/rHx50lOH4BP\/L4z8UPeVtazS6Ul0afRL4eZ0lEr2fR1+Npyo0kBwFcjkqQRFK5FNet+VhM\/X3nHHHYngrlu3jg0NDUkFV9W75qN2Kp+shxvqeUzhbREL3xyaNhQQXAgu+QiJ6Y0zW+rETCZYJGhiaDlkW0k20re5mMrsI1GmH+8U6OwL6YfRCq5pYzFe3ytaSfic7KlTp5h4og3Hky\/koGebN2+ug5n+tmvXruRvRe0Yr\/zFhB04tK+FMWDnwgZahDT48jvKhVEiQmlBs0cvX0oSejGUnJU63SnIetcFfqYlKV1wXc098J6OajRl+9zl3INoW7pnVvYcbtZCDlGE84SUTSvheH0vpKPrGhz4Ydg5XJd13sUcrjgvarLQkd7no1yT9120ozS6Tc8fq9pRmm\/mo1zM4f7\/2oZQVlyhLNVCDp3gIqQsb6x1DVFMgottQXIJjCkcqRJpXT3SlYGe3\/Kd07m2\/YhzqSb5Z50ZYJKeBL77ukljwt1pwRUx4qPcyz58NnOLqQ4\/l52j0ke4cPTqOrouLGriSEUdscrpQzq6boQLP6yuH+rqkc4Pbxp4gX3+qj+rrfpNI6FKT4JGe191+Rd9vvuZn7HHXjmXq0NAZeCdgt62Jggub6wRyqpuKEsnuC57ho2Yl66hDFVmFyFJ11M56fxchCRl00qxTe3YcF6Uv6s2HmWHb5+Se+88hW1pxTAJmo7\/IvxRh2D25POLufJOvVGn4Kf3TFOWLXQbFu0I16bijac0aKzNw7hFnF03\/1PkeUwcYi+1XesRA4dFbJBtszFFgg7CoAMxVCc8meajey99yIXuffG5auW1+E4R\/PLYQu9GK7hFwxCNnj5kJUE40i7cqKuD4LD+pC8ZyjoMy34eA4c6G7IwSo8A83JAgvaVayawv7\/508pmoghH1CF4b3SU\/bBHfqxkOkKRNsKkU6DDL6+oZr3vXXDFhTjpi5OxaCquRVM2FQscVoND+GG1O022\/GUJEo1O3zozynbd3GJ9bj3lceTf3mU\/3zjfi+BSh+DRxZNY57wLB12kP6QTdMqDwtG06Er2axjBpa0mWRcnF5170K0CdfEc4chsGZ7UuZktXLiQTZgwIXnx3LlzyX9j+rfKFm5vWc+ffPJJdvr76236ObnSwA8vXISSdw4wLVjp9CEaa1\/80cpfEqK2ie8l9cn2gBjdPKnttAvfG0xbgYq0w7f\/08tJe6SaAw7BIXdYryNcfngC3fXY3Nyc3OPY3d1du4Q+ZEFztVAVeDkW7D5+2yOs7frrK4BYfCa++NJL7N+\/dad3w+CH\/iAO4Yc++KNQK922o7oBKA9iReaBs75D+d54xSXK1dN5bOQrqht6hEsVhR+UQAVNX5yMcGRc4UjV9Xx0qcGePXuSujp58mS2b9++5L5O+oHDuDiUNSjwQ3uOdOHKUIJr246mR+i8frjeR0ujZdvr8VQY88VSOg5MntM+Y9mxlbo2LI+om7zrfYSbVVGo14GfPQKuVwfKruejm4LSkQnR4qxpAR5+U4WrQj23DZW5mJLICtWFaKzp+zrBhR\/a+yCldO2HaWt0\/Omm5o4ePcqmTKnf9iMKpOy5WG9Nnq96jrF1X5iVzJOm\/YbSz59\/fo5XFtKXPRdPirJJL+7Pp\/RvN81KVlTTFiGxbA0nuDQ6Qki5mEOHCoPIrucT54\/4qDYtuNhSYsdvSMGFH9pxpEsVgkPXIWUSnh+\/MeIknMzx4Quw8t5Fq8I3z1nIOo74c9XiqRAcchu8jnBNJvtxwo28uujCJK4riep6vuPHj7OlS5fWjBTv6NT1DnVlGO\/PXXOoanjgh9UOKRfhj48qxRFfWnhc+aFKJPPmn55fzps+7Qc8PXUK6Cde7adrw0zF3PQ9r4LLw1m8wZZdnAzBjUNw01bw6\/m2b9+eLHijHw8vt7e3s66uruRvCCnbr+4MJbjww2oLbhH+0oIrWyzlStBUi6fy5p9eLJU3vUpw6e+yxVMh\/dC74GYpv27uwfccmkn+RZaj2y6H55hlpS9aSUyv56PpALqEnv8oNEm\/vr6+muAipGzav61\/ryiHdl8dmwp+aL9tKAYO8\/BHYnbJR95nGz8\/0XobUFa7ybcIUS0T28487WjH42\/X7et12Y7SMZZzL2tij6+4ruYIITksXXDRWNs1m64riep6vgULFrAjR44kAsvDzvT\/bW1tEFw76kpxdF3HF35oR6ZrP7SxIo8NWdtjbL6dTkMri+nOWtUhE7pv+L5rV3byVB78dPbrnkNwdQhF+txHJVFdzyduC5LN4eICClxAEambeDfLhx\/mNdrUBtcLm2R2Fl2Q5WOxVNrO9By2KX55eZG9X7rg4i5V80P4xTBNyEpiOzpyOfcis6Hq+cfEIfywun6oCynzbT380ni6H1ZsS0y2\/aS3FWWlTx\/HaLqth06W+vYLb9XdXETfMU1P75psO9r4r++ys7+\/iPEV1SH9sHTBRSjLrt8UspLYCq5dycZPKnBYfa5j4NDEBpcnS+lYU60G1qXTnXmsS5\/nuRhaN8EvT95Z70JwXSEZOJ+QlUQnuAgpI6QcuPpH87kY\/NDEhqJzq3kBt5krLnINX177xLliE\/zy5q96P1rBrXq40Lf9ISuJTnCxtUuOkK4OgENcz+eiIdfVI6qHdLSh6txkXT21eZ736EgKd9NCq\/QeWR4mFvcRpzGzsY\/y4PPFt\/3VNSxU4Z+lAAAgAElEQVRUpLV0wcXcUXXnjqjSZjm7rSNwh2r09LqG0kVjbJJH1hwgnxOjfGRbO0I9H+\/Hc+o6vVlRJhKWjisn1LbCpOc5Xc\/h8vz5nPEH7\/xMe7TjPUcuSjoEJnOwYl2k\/7ed4+Wj3J8PfGn8CG6onoVJw1Old6rQWPPeqaqxDvU8zx7AtDO73AOYbkxi4hB+aOf9MXCYZQO\/4s7FrUB5ETLd4kMdAtpK1Lfw\/DnHIX\/07QnP38dO\/uRfgny29BEuHN2O5xgcXTfCtSvZ+EkFDqvPdQwcZtkQciGSjE2dmBbdRlS0BlGn4NjPfsXefeRvimZllD6I4KpunNEtZzc5Ccr36Gk8jY5U1\/NlnUoVQ4NjVNMjfCk0dvBD+baRIlEMWw5Vviaead7R0cH40ari+9u2basdrco7vbKQMr9xZ9fNLUntLyssL+6tTbfptFBq\/V+2sM55VyY2ykLKRfjheWa14xPv\/B5r\/\/QVtW1CPpsK74IrViCcpXyByqLzk7aOnlWZZNfzpW8LokMw6PQp8aQpzMPHPw8PP5TX\/LL8UOZrZCG\/sWvRokW1azFnz57Nenp6WH9\/f1KIrVu3soGBgdqRq+kVwSfPjCaXy9ONQIfuupr94Tfv1MRWhkJRDEzSpxdt8ROfKJT85dkflGqfj7ZU1c56FVzqUe\/evZt1dnay3t7e5HhA3lAjHFmsH+Wjksiu56OeNT\/aUWaxDzuKIVOd1KGwgx\/6qxO2HKquwly+fHmtneSiTB3crOsVp372v7Bz89awqa1NbOolTez5N0aS\/y9j3laFNIVuuV1k48mzo9HYZ8uhTa3yKrjcINkZvFxwbYxGmvMIuJz\/zrqeb3BwkJ08eZINDQ0xMcwFDovXRJcc6qyBH+oQsnuel0OVr9HX+UiW7p7mN3YtXryYHTx4MAkv02\/t2rUsfWMX\/X30U19kf7z4UnbxK4\/ZFSRAKrKRfk2\/\/EGAr5l\/Ii+H5jnXv1mq4NoajXT+EeDOfu211ybOTrcG0TV9aWf3bwm+4AoBleC6yh\/52CHAfW3NmjVs\/fr1Seg4j+DafRWpykDAqeCqFtfA0cugNvubptfzrVq1ij399NO1xRvp6\/niKxksgh\/GVQf4YrXDhw+zyZMns3379iWCyn80bUN+tXPnTrZ69ercIeW4SgtrshBwKriqD0Fw46+Equv5VqxYUQtz0QHmGOHGzyX8MG6OVL5Ga1xsFk3FXVpYJyIAwUV9qOtpL126NPn3nDlzkjAyXT4vbklIX88H+KqDADq+8XClugpTnN8VfU18P73bI55SwRIdAkEEN21E1p4yncGhnqv2LIb6vuo7YniK3ilLAMGhfU0Ah+bYwQ\/VWFXBB8l6cHiBw+CCK+7rJDPSe8rMXdHfm1l7Fv191Sxncf8et3PJkiV1G+HNcrJ\/CxzaY0cpwaEZfvBDNU5V8EGyHhzWcxhccPkCAXHVa3d3d93+XDN39POWbs+in6\/a5yrbz2efm1lKcGiGk+lb4HAsUvDD7NoTuw\/ykW3WOQym\/hHqvRB+WIrg0r5O1Z6yUODqvlOF+a70KVC6Mrl6Ts4ODt2gCQ6zcYQfyvGpig+Ko9z0wUduPMhNLqH8EIKr4Ct2Ry9zXqQqzg4O1Y0ROHTTUJflh1XhrwqCG5LDUgQ365gyN25QPJeYG+tQvTEVilUIZ8Xu6ODQzEfhh+oRbhXaUfhhPX\/BBbdqk\/2xhUHSe\/jMmi23b4HDYniCQ3P8YhXcsjmsig\/GLLhlcBhccImAKuwpi9XRqVe7Z8+euhYrfV2XeXNm\/yY4tMcOHJpjBz\/Mnhbg++Zj3psLDi9wWIrgmrsb3gQCQAAIAAEg0BgIQHAbg0eUAggAASAABCJHAIIbOUEwDwgAASAABBoDAQhuY\/CIUgABIAAEgEDkCEBwIycI5gEBIAAEgEBjIADBbQweUQogAASAABCIHAEIbuQEwTwgAASAABBoDAQguI3BI0oBBIAAEAACkSMAwY2cIJgHBIAAEAACjYEABLcxeEQpgAAQAAJAIHIEgghuyNsYIse7suaBw8pSB8MjREA8GpXM48ezihe2d3R0JNeYNjc31x2HW8ZRrhFCWEmTvAuuWIFiPu+zkuwFMhocBgIanxk3CBw4cIANDw8zuhxF\/PFL0BctWsTWrl3Luru72ezZs1lPTw\/r7+9PXt26dSsbGBhgra2t4wavRimoV8GlUdHu3btZZ2cn6+3tTSpXW1tbo2A3LsoBDscFzShkYAS4sHZ1ddW+nD7kn4vyggULWFWu4gsMY+U+51VwORqx3hZRObZKNBgclgg+Pt1QCIgRIyrYnDlz2N69e5My8pHsrFmzGAnusWPH2OLFi9nBgweT8DL9aOTb3t7ORLFuKIAauDClCu6MGTMaGFq\/RXvuuefY9OnT\/X5EyF0luODQnoI333zTPrHDlODQHkwXfsiFdc2aNWz9+vVJ6DiP4MbG3++u+Qr7k9\/9mjX98gf2wHpM+fupN7L3J36SXfzKY8lXQvph6YKrKuyJEycyBWW8PycnC1lRsgQXHMpbB10d9cGhajEOt5A37nwxDv09yw5dGcb7cxccEmcUMt65cydbvXp1berNNKSssyEURzueHWY7nj3BprY2samXNLHn3xhhZx74HAv1fZVG8++fPDPKFn3z1eQ1bh\/Z+tvHbg3WlkYruB47OA2Rtc7JXBfSRnBd29Bo+fngULUYh7B7\/fXX2bJly9jcuXNrq1+54IojNWqg6McjKDH8W2UL2Sk26GlbfT+34ZB8adOmTWzDhg3JwicSW\/rRGhebRVM2Nvjwpdav\/ygRWPFHf5s3s4UduvtqH5\/MlednNr+QvP\/a+htq6Xgn4dBdV7N5V7Tkys\/m5dIFF45+Piyct6EI7WRZggsO4+FQthiH6hctftuyZQubOXMme+2118YIbshoiU1DFWsaWz8UIxF8DpfEV5zfXblyZW0Vs\/h+ercHidpP75lWagdp7oPDrG\/hdNa3cNqYtozsUz3z3SHi+X\/tmRF28uwoe+rWSXWdSd7hDFX\/gwiuylkQylLPwerCMLaO7rrhAofxcKhajEMNOY186Xf55ZezwcFBY8HV1cPx\/jwGP5x45\/dY14Kr2EPdV0rd2zdHV208ymZ8bIJyFEvfJ0FOj365sb7t2\/3Mz9hjr5xj31hypXQUG5JDCK5CgXxXgqL5h6wkWSINwY1HcNM88fnaO+64IxHcdevWsaGhIangqqIUPPJCecvCuqGexxTeFrGIwQ\/JhpG\/3stUYdGibY0uPQnuM39\/XTJ3K\/tRehph0k8WWtblX\/Q52bfgkxOVHZKQHEYruK5HYo2WX8hKYiu4jYa56\/IU5ZDPyZ46dYqJpxJxO\/liHHq2efPmOvPpb7t27Ur+VtQO17hUKb8YsCMbrup5UiloPvGkcHH3dZOUYiZ+m94NNVcqfpfmbsV52zQeITksXXAx\/xfP\/J+NY1JlBYdxcJi1GEcU4TwhZZs6MZ7ShGysVbhyH6SwLQnaZR8+m7zqOyrw3eMfYoMvvzNmXlS2yI5sWfTQq+y90VH26OJJdTtQfC56I5Ff+dkWtm3J+UVbsrUyITksXXCxpUTuRrowSshKYjvC1ZVhvD\/3waFqMY5OcBFSlneadCu0fXCYt8PCbSBBk4VtffkZH932tjUZb+GUrWT2ZR9hQaNbWiiVdWZBSA4huIra7bMS8J5WViXQfT9kJYHg2nWKwGH9lh0ZijqMYn8egx+KNoQUNB6qzcMRdQo+0dpUF4LOkz5PHaJv3XjFJezLsz+A4PK5I4Qj4whH5u1V8\/cRUr5wQEvsW7uyQpKhtkXY1rNY08UiuLwdJZG55CPvs42fn+g1pPzVg6fZ7MktiXDqogDpkDGFvmlOlS+y8hVSlnUIEFKO5Hi7WB069kYyhganatyJnZUYhA6dpmp3mtI+KBvluvYREjPVNhvdtyjtjTPPi7Wvn2wkHUNbipAyQsqF6jy2BcW7LciUWHBYbQ7T\/JGg0UETtHrYxfRVOuSbniu2CQmLnQKb9GLdlqXPk3\/IQYN3wRUXcaQvToajV8PRwSHmcFXi7aOx1DWmMT0P1Vjn9UFxK4xrjtIjaJv88whi3vzvHvwFe+vMaG3Pry59KA6p3noVXNqmkHVxMkJZ8YeywKE9R+nRRXr+KKSjZ412J3VuZgsXLmQTJkxIXjt37lzy35j+rbKF21vW8yeffJKd\/v5602CC1Xs2Psi3CNH5wCZzrGSYyTYiOnv42y+8xQ7fPsXofZUP3P\/iaIIFnwN2eVY2lZ3ylY3wG3oOl2+8p7sem5ubk3scu7u7a5fQx9LgWHlByYlCYQcO\/REdikNdCT5+2yOs7frrda\/huQSBF196if37t+70io2ND+aZw8xjPF\/5S2cmF\/3RKFdcPFU0P55ed9BF+jsh\/dDrCJcqCt9kT4VMX5yMkHL8IWVwaM9RGaEs1fV8dKnBnj17krZm8uTJbN++fcmdq\/SDH9pzHKKxtvVBLmh\/+M07xvtkZaIn1mMSMzpcQzzGUVfPVc\/54qk8+3h19tlsOwrBIbe7VMGlCoGfPQKqw8DtcxybUufs4LAY2q45lF3PRzcFpaNLotVZUzs8BEfv4yzlsScVhWisdT6o4u+W75xOFk9NGX2dTZmiDgEfPXrU6Dk\/D3nXzeevseP1wST9\/PnzkzTpkC61Hz+45SKmek5pKH\/T55QfnWTVOe\/CCmhd+hAcBhNc6lkjpFysUZalDlVJbMJZ7kvbmDn64FB2PZ84B8hHtWnBjWF7UhVZ9sFhGgdbH0wvHiqKr4\/tRi7zpPL++I2RzHOTy2xL6dteR7gmk\/042lHuBrowTQhHJ8vAoX24MTSHquv5jh8\/zpYuXVqraOI9qwgpX1gUpwtXltVYF\/HB9F25NmXk9Vg1N6qr51nPSSRpkd7jK65T9gdM86eyiouleIa69KHaUu+CSx\/IujgZc0f2jXnISgIO7TpFZTs6v55v+\/btyaJF+vHwcnt7O+vq6kr+hpCy2QpdWUg0lB\/a+iA\/eaqooFE4mY5JlC2W0tVz3XPdKFeXvmiHIBSHQQQ3K4yBbUH2W05CVhJwaHf8Jm+gVVseinJoej0fTenQJfT8R2Fn+vX19dUEFyFlu4BrUQ7tvlqfSteOpke5NtuEOh5\/mz3yt59OLnC3SZ+17Ue8r1a2bUcUXNVz6hDQGc20AIt+4jn1uvQhOfQaUtZVppAF1dlSteexYBeLHVXjTxfKtSmP6nq+BQsWsCNHjiQCy8PO9P9tbW0QXBughTQx1H+dDemTp\/IWWXULUd58st7XjXJ13ypy164OP9238zyH4OZBK6J3Q1YS3QgXoyO7iuGDQ9X1fOK2INkcLi4RsYti+OAwb23S2WC7mIjbUUTMTMuimn81SV90cZgOPxMbTN8pXXBxD2d17+HkozRwWH0OsXhR3mTq5g9DNtaqRl0XUqZtMV986n3Gt6Clw7JZ23roZKl9R4fZvyy\/XHmyVJFtQVQmSv\/E2xOT4xjTW47486xtQeljItMh5XGzLUin+jFUVp2NsT6PBbtY7IiVpypEB8Chfe2JATsTG2gUSL+8N\/QUGXnmRdV2JJ33ZKm0XSb45S2L6v3SR7gIR9pRGbKSVEE07FAsN1VMHCKk3LghZdvQ8ODLp9mOZ0\/k3tdq61U2x1GSSNOZyXk7EqKNIf0wWsHVhXLG+\/OQlcRWcMc7R7ryg8PzJw+JK0rTdS325zFwqLOBY6ha\/KTCmIuZy6MXZW2J+H3ZKFdlH+8QPHXrpEJ1SIefbQdClq50wcX8X\/Xn\/8Bh9TnEHG5jz+Hyox35FiFxnlM1B0tHQ9LlAiZztLqjI02PZrxp4IWEiB\/23FAjRDUHyzsES6a8a3z0o2xb0bgSXISU7fpPISuJ7QjXrmTjJ1VMHCKk3PghZfIs0y0+JGZ0DrOLW4HyerTpXG7RuVtuV0g\/LH2EC8HNWx3Pvx+ykkBw7TjSpQKHOoTifx4Dh3lt0Ilp0W1ERVmjTsHJs6OZc8cuF3Plxa9I+YIIruq2Et1ydpsTTcRQiYv0Li9GFm2j\/9edgJL1PGQlIVvB4dibTmLkUHU9X9apVKHrUpEGK7a0ttipeBLPw+7o6GD8WE7x\/W3bttWO5eSd77wRiqwL6tM37pTRjpJ982a2sEN3Xz2mnaRTr2Z8bELyrKgPhh68eBdcsQLt37+\/drqNrqCxL5Yo2z5bR7dpsMChHLWidcAHh7Lr+dK3BdEhGHT6lHjSFObhw87Dy3iiWsZve1q0aFHtSsXZs2eznp4e1t\/fn1TErVu3soGBgdpxnbp6pKqnPCTLn9OeW1qVnA4lF63ntun5\/loxPdlMP\/FOXtv8Gy6kTKOi3bt3s87OTtbb25scLcedXCe4NsIwntLonMwVFuDQFZJj8\/HBoex6Phod8aMdZaXxYYc\/1OLK2RY71TWKy5cvr7WTXJSpc+TjmlMK3T7\/xkhymfzUS5q0YdzQyD\/\/qxG26Juv1uzjttJCLpc\/Ww5tbPA+wiWjZOe3csG1MRppziMQcv4bHPqpdS45zLqeb3BwkJ08eZINDQ0xMVQJPyzOa14OVTyRJXwkS\/cW89ueFi9ezA4ePJiEl+m3du1alr7tybYUf7x4Ivv91BvZHy++lF38ymO22XhNN\/qpL3q3Ly+HtgUuVXBtjUa68AioBDe8JfiiKQK8wb722muTBptuDaJr+tINtml+eM8PApynNWvWsPXr1yeh4zyC68cq5OoDAaeCq1qYgcbaB3V+8gSHfnD1kavp9XyrVq1iTz\/9dG0BTvp6Ph+2Ic8LCPAFh4cPH2aTJ09m+\/btSwSV\/yjkT5zs3LmTrV69OndIGVhXBwGngqsqNgS3OhUCHFaXK9X1fCtWrKiFKumAAoxwy+VYxROtcbFZNFVuafD1PAhAcPOgNY7fRaepGuSrrucT\/56+nq8aJWssK1U8ifO7Ik\/i++ndHo2FTGOXJojgpiHM2lMWC9yqfadl2yeGp8iWshpPcGhfE8ChOXbwQzVWVfBBsh4cXuAwuOCKewLJjPSeMnNX9Pdm1r5Tf181y1ncv8ftXLJkSd1GeLOc7N8Ch\/bYUUpwaIYf\/FCNUxV8kKwHh\/UcBhdcvkBAXDHZ3d1dtz\/XzB39vKXbd+rnq\/a5yvbz2edmlhIcmuFk+hY4HIsU\/DC79sTug3xkm3UOg6l\/hHovhB+WIri0J1C1pywUuLrvVGHOMn2CkK5Mrp6Ts4NDN2iCw2wc4YdyfKrig+IoN33wkRsPcpNLKD+E4Cr4it3Ry5wXqYqzg0N1YwQO3TTUZflhVfirguCG5LAUwc06psyNGxTPJebGOlRvTIViFcJZsTs6ODTzUfiheoRbhXYUfljPX3DBrdpkf2xhkPQePrNmy+1b4LAYnuDQHL9YBbdsDqvigzELbhkcBhdcIqAKe8pidXTq1e7Zs6euxUpf12XenNm\/CQ7tsQOH5tjBD7OnBZYuXZq8EPPeXHB4gcNSBNfc3fAmEAACQAAIAIHGQACC2xg8ohRAAAgAASAQOQIQ3MgJgnlAAAgAASDQGAhAcBuDR5QCCAABIAAEIkcAghs5QTAPCAABIAAEGgMBCG5j8IhSAAEgAASAQOQIQHAjJwjmAQEgAASAQGMgAMFtDB5RCiAABIAAEIgcAQhu5ATBPCAABIAAEGgMBCC4jcEjSgEEgAAQAAKRI+BMcMWj\/qjM\/LhB8e9lHEEYOf5RmQcOo6LDyhhwaAVb9InQjkZPkZGBzgT3wIEDbHh4mNFh\/\/xXpQO2jdBq8JfAYfUJBofV5zBdArSjjcOpM8GlA9mnTZvGurq6auhU5Rq3xqGzWEnAYTH8YkgNDmNgwa0NaEfd4llmbk4El98GMTQ0lJRlzpw5bO\/evez48eNscHCQbd++Pfn72rVrWXt7e50ol1l4fPsCAuCw+rUBHFafQ1kJqnTZfGMy4K5UTgQ3bQ6FtY4dO8YWL17MDh48qBTcGTNmuCvJOMvpueeeY9OnT\/dWanDoDdpaxuDQP8a+v+CbQ7JfJ7hoR4uxHIJDbqEXweUhkFWrVrGHH344Ge02NzcnI9zu7m7W1taWfJ8qyptvvilF68SJE5mCMt6fZ2FXrPqdTw0OGfNdx8Chf4yrzqHoizbtKKX3jUHV8\/fth2J77ERwKZS1adMmtmHDBtba2spoHol+K1asYD09Pay\/vz\/599atW9nAwEDyjk5wXYhGOo+7B3\/BBl8+nfx53swWdujuq318xjrPPPa5riRFOFz58H9nO549kZR7amsTe239DdYY+Ei449nhKO2LiUNVxzeLj+d\/NcLue\/YEe\/6NkSh5D2Gfaw5leOsWTYWwwYdfxpJnSPycCC7vhS1dujTBkM\/hkrCKy9n3799fG91ywRWH89RToh8Plbr8d+vXf8Qmf\/Qidvj2KUn+n9n8Ajt5ZpT99J5pmd9T2ZLuOaZtzfuc2\/fzjfMTDLh9Zx74XPLvdP4+KonIlSmHH7\/tEfbxGX\/OuN33PvEq2\/OTEaay2yWnMlzS+d808AL7j\/dY0gmgZ4++NDLGvlAcx8oh1SUbP5z74DDrvm4S621rqsPVB8c2HJFPcfu+e\/xDSaeL\/F1sY\/L6aQgOZUIUSztq4nMy\/m344ziII+ii7awsvY+2VNWZcCa4Nr2VUCFlEq8bZ7awh7qvrJlJwH\/tmRF28uyodEQWMkxiY1\/ISpLF7cQ7v8fefeRv6l6hKAKN1h9dPIl1zruAeTqfEBjf8p3TY\/iNxb5YOLTxQxKzvoXT2Zdnf1DrsFIHluqyGD0KwbFsLUMo+2LgUGdDWRzIBFPWlpRtnw4\/G20bt4K76KFXpaLKSabn9EuHl0NVAlv7QlaSrAqnsoOHcflItwxHu2rjUfbI336azbuiZcznY7CPOiv\/3HuT1D6XTq7LK6\/gcjHrWzhNOj9Iz6lzS6PLUH4kljGkfTH4oc6GMjgQ+Yj9+zr8dP6T53m0I9w8hch6l3rc31hyZWajRg566K6rS2n4bO0LWUlsBJfSqDozrrjNyoeHE8WoRvr9su37yMkfs9PfXx8Cjsxv5KlLqg6i+AGTzoyvQoe2Lw92vspsOyXgY+qOyoiQsprp0gXXZu7ItKLQCGfuZU3s8RXXJQio5pbuf3E0WUyVnt\/xPXdA9i345MRkNJA1N8EXU4nzzTE4OmGqc3YSPjG07GN+j+wQ68TGf32XDf3v92vztlnzR2Jny8f8kKwOUX378Rsj7LeP3apcpe+rcZblq+NQxIXzSX6V5Ycdj7+trNs2dcB0DjBd37jfp9PTVAOfZiri5zH4YQw2hKyvrr8VEr\/SBdfntiDZ\/B0nKx3mEMNQaSeVEewiTFLEvpCVxHaEy+fJKb1sRbgLDFXzdxSxuOzDZ7Vby2gevyz7bvura6IRXBM\/lEUEsjiUiV+6LrmsA2XYF4Mf6mxwibGvtjDrTAHf9uvwcynwDSu4upCijERKI66u9VkJitoXspIUEVzCUJzTE\/Py4UhiSNE0\/7LsqxKHxJts+iMLY+Liko+8X4sw+W6sy7AvBg51EYqjR4+yKVPO786QRfpCPJ8\/\/\/zuC1kUib5f5vOQHJYuuD5Cyi++25xsAXjq1kkJyaYhaGp4V362hW1bcn5\/bpFQU1Z6Cl9z+0xDZWQLbcOglaG0WCVkJdEJro5DCvGe\/f1FySjXJpyYl8N0CFuHMYV43zozynbdfH5xldjRcl0HiEO+XqBKHH714GnW3NSUm0PTELSOYx2HPFKRl0PdlIeuHYiBwxhscDkKDJ1XSPxKF1ybDfc6QnSjx6z04ihX9x3b5y7sC1lJdIJrwmEIXGmum+ZGbQ7eKMO+KnFou7BQtuXN1m90fmuz8LGofTFwGIMNPjgNlWdI\/KIVXNNwoIwUciIa3dqEhHl4sW3ie1bpuT1Z9ruwb11Xe+Xm\/z7R2jRmL7QNRyqM06HhPHWIwp+h7Qvp6LadJhfz8DSqV20Py8ORKiRddB6+iH0hOXzvvffGHI9LmOhscIGxSz9N81i2fTr8XAp\/wwkubzjp9BubSkKjJB5etEmvE1xX9v184EuVElzChQSRRp90\/GM6VKdqTPNwQB0ZcXSb15FD2xfS0YsIbhFBIg5ocSBNhdC+XB+NbZn2heJQvAkq74l9IeZodXPEZc7R6uaIQ3FIdd+74IpHkm3btq3uaj7dZL\/NfJ84R2aTnkDhDux6\/o7nzUNftvmTMLT88\/JgguuKQ76IZuPnJxrPq5tySKOcG6+4JDn5iH7inJ9u\/o8\/p2Mg6fhPvo0s3SmQLfgw5ZDK\/ueXMvbV61tqtoV0dFsOqQN65N\/erR2JyjFJY6z7tyiKppyacMjt40eL2nIk+rxYFl0dCMEhjWx3797NOjs7WW9vL+vr6xtzRK7JtI7LkVoj5RWCQ46XV8ENfei2q4MMZOFFFxXMpX0\/+dHTQQ5NcM2hr7nS9OjWlq+Q9oVy9CIcqlZw58XXdg5Y952y7QvFIeHAR7kQXF2tyPc8JIdeBZdf8WZzrVTecCAPWbra1qNreGO3L1+VU7\/tmkMfW69cHs8Z0r5Qjl6EQzqchY8ei4T9aSRKv\/TJXzZ+JNpRtn2hONQJbtZOAYSUs7cdheTQu+AODg5mXkBvsuHexNH5Vhs+h+fCkdd9YZZ03ikdZirDPjqH98Hl85X2uRRc1xySQPItQnnmaNNl4hyrRk+2dSCUfaEc3eQCc5kfut5HK+vE2nJEdSEG+3xw+Prrr7Nly5axU6dOsY6OjqT9pPvEbUe4RTC2aetUfqpqk8q2zweHqrKWJrj8xhbd9Xi6uSH+XDZHZjNXxPPj9qmumdPN7aSf0\/zgtZc11fb45k2fLkuoSmLSWOv24aY5dDHPLs7v0TGCqjk80zlWzgf9N22f6Rywav6QnyhWNQ5V85q2fkXz7LQSnBY0mvp1lp9wAbfhWLbXWrbGQOenofxQN8LFHP9xQ8YAAAyxSURBVK79ECMkh94Fly6jV4WUs8K2eXs9RVeoynplRVc\/ik7t2r5QlaRIOFLFYdGV2pwryt8lRzzfUPbFzCEP09NBEi6iEBzbtM\/n9XORI\/r\/su0LxWGW4BKm+BVDIOtWs2I516f2Kri6xRqTOjezL33pS2PmdfKGMfKeoWqaP51ARD\/ZjTN5Ggof9oVydB2HWXZkYUSNhBjdkFVqE4zznEct61SpxCSEfTFzyMP0JudR5xHkInulRf5isS8Uh1mNvs4GEz\/Kw2EePzJta8v8vg6\/ygguGSpuR5DtHxv56721vZm2oSrXoS8xhEh507zwH37zTl0YTBdqEp\/zxlsMo+VJz98V04esJDoO84aUqRwklN3X\/Zl0G49Yzqw6YbrNyCYkTBEJbp9NeipDeo6xrJCyiR+mOdSF6U05StddWjx1\/NRIcoOU6ZGrspCx2NFyFVImW2\/\/p5fZqd++z37Yc0Otnc3KP6Qfqhr+GGxwKUqh8wqJn9cRrg44KuhVPU+OOeFHl058XuQ4P5PvFN0i5Mu+kJWkSO86K61uJbiOH19bTfh3fdsXK4e8s0Fndvv4FcU1Jvti4JBssOn05unw2Hawig4sQqQPyWHpgkuT\/UUc0NU+PF\/C4Mu+kJXEp+BSuF52ApGuoXe1p1nHu0\/7YuXQ1Z5mFbZFO7Ex2RcDhzHYoPPXmJ+HxC8awaVGV5wrNZ13UDmfaXpVRRDTk2jGZl\/ISmIruDoOKHzHL4qXfSPEHHDW3JFv+2LkMC2GOg5tn9vOwcZmXwwc6myw5Yj7ZKOn1+HnsrNQuuBSKOS7xz+UXFdHi2jyhC7IaTuunFA7is92DjgrtELPdjw7XLMvz9VtolCrtozYzg+GrCQ6wS0SziKM0tfp6eoAXRX3H++xZG7dhHNbjMkOMWxtyiFtfzl5drRmn+r7MXLI+Zh72YWzyE0w1nGWfk5zxAs+OTHpZJvkzzFM1xdfIUeaI75xZkvNvpg5jKUeuRSmkHmFxK90weX7x8iR+F2vpmD7Di2JdsRmX8hKohPcInsAbULDvsL0snL6tC82Dm3KauqrsveIR\/EyC11eMdoXA4eYwz2hXIBn0iELyWE0gstHE6ZHM8rCvKLD+giDiHPNuvx92xeyktgKrg4j\/ly1+EmWnhpd2egxb0g67YhZ6X3ZFxuHsnKacmgyNZN+J+9JUTHaFwOHOht8cpjHj2zqSIj8dfjpOoJ5npcuuGI4ks5F7fjUhOQ0Jh5mkoWqKAT97RfeSm4wUYWystJTeMj0efpEo7+74ROMVm9mhcEoBM3t87USMGQl0QmuKqRsivHGf323Npcr4ipLL4YUTfMvygHtx6aL7emOZV3oNM8pWjFxOPPu77DZky+EUHk5TTG29UNdCJt\/v8gpVTrOsp7rpjxi4FBnAwT3wghY1pbp8MsjqLp3SxfcdDjSZKtHyJBiGsBY7AtZSXSCWySkzPM2CdmHDimK5fZhXywcTvnLv2MTP39X3V3CuobD1XMTTtPnpLv6tkk+OvtCcMgvnj98+HBi8sqVK5Mr+vgPIWWElE3qMpNVVjFkKMtE5wBGHy7wUiz2hXB0E5hc2pHVmaGoweDL75QiCmKngN9lnMbGxj6X2JlwpXrHpDNRJH9dWt33dc91+Rd9nvX9EBweOHCADQ8PJyLLLzBYsmRJ7W7xEDYUxTDm9CHxi26ES8RQBZ83syW5UUb8+TpEIm9liMG+kJUkxAiXf0O2EI6LGS2qs9mzm5ffrPdd2teoHNrgLetsnTwzyhZ989XaamGbfF2lUdl3Tf9h1rXgKunxr66+PaZzt2MHmzZtGgTXEcAh\/dCZ4IrH\/xEO27ZtSyqE+Hf+N5NQCFXwyR+9KJmnpTk4aujef\/\/92r99zy3ptipw+\/gtNTT\/TL+sW2voedH5RJ7eRyVxzaENR3Sr0jv\/5\/0ER0pP87uHf3GOyW5tMs3fFeb0vfSWJNE+XZ1JP29UDnXz8KrnNP9Nnaovz\/6APfrSCNvzk5HcWwXzcpDnffL5lZ9tYV+9vqW2lbHln5czF1Mqptohnm0+a9asJJkupIz7cBvwPlwx7MErT5GD7ykPGtHS\/M3U1vPXefG7bk0rp+\/3yrTPR2Ptg0MbDvi+Z0pL3MfGuyv7GplDG96f\/9UI+69P\/ILRyDZG3mX2+eBQhR2fy+3u7mZtbW2110LaYMNr7GlC4udshEvX8IlhDgK5yNVusZNUtn0+Kgk4DMsqOAyLt4+v+eBQdgE9iW1PTw\/r7+9nfGQrRgpDjrJ94Fhmnj44VJXHieDyifyhoaHkO3PmzEnuwD1+\/DgbHBxk27dvT\/6+du1a1t7eXjf3UCbQVf+2SycDh+XUBnBYDu4uv+qSQ5ld5JubNm1iGzZsYK2trWNeIcHArxgCvjnk1jkR3HRRKTR57NgxtnjxYnbw4EGl4BaDCKl9IgAOfaIbJm9wGAZn31+hyNOePXvqPpNeD+PbBuTvBgErwZWFPJqbm2sW8VDyqlWr2MMPP5yMduk5jXDT8w9uioFc8iIADvMiFt\/74DA+TmAREMhCwEpw0xmmQx7UI6PfihUravMO9O+tW7eygYEBaVgENJWLADgsF38XXweHLlBEHkDAHwJOBJfME7eU8Dlcmm8Q\/75\/\/\/661XX+ioWcbRAAhzaoxZUGHMbFB6wBAiICzgQ3D6xZe3Pz5OPzXdUSfJ\/fNMlbd8ybSR4u3gGH9iiCQ3Ps4IdqrKrgg2Q9OLzAYXDB1e3NNXdFf2+KK3ZjG5Xrjnnzh8qFnMFhMZTBoRl+8EM1TlXwQbIeHNZzGFxwdXtzzVzR31vUG9u9ezfr7Oxkvb29yfml4iZzf1+2y1m2d9YuJ\/NU4NAcK5M3weFYlOCH2TUndh\/kI1u0pREIbtbeXJMGKsQ7vGcWs+DKjnkLgQ05Ozh0gzQ4zMYRfijHpyo+KI5y0ZYyVsoIF4118ca6zHmRqjh77I01ONT7ATiE4OprSbE3QvphKYJLIbTY9+bG7OhljYp4ta5COCv2njU4NGsk4Ydqwa1COwo\/LDmkXLXJ\/tjCILpj3syasWJvgcPi+GUd1Vcsd7PU4NAMJ9VbZfthVfiLWXDL4DD4CJcIqMLe3Fh71rEc8wYO7RtscGiOHfxQjVUVfDBmwS3DD0sRXHN3w5tAAAgAASAABBoDAQhuY\/CIUgABIAAEgEDkCEBwIycI5gEBIAAEgEBjIFBZwRUPC\/B9cADlv2DBAmcHYNDcy5EjR5JDNcbzDxxWn31wCA7zIDDe29KGENw8hOd915c4uq54ecsVw\/u+O0q8jODQH9vg0B+2oXIGh6GQLuHgCxdFo7No77333iQruoh5eHiYTZs2jV1zzTXJFYBz585lDzzwAJs8eTK7\/\/77GVWooaEhtnLlytqoUrxLVLzdSLSPNkRv2bKF3X777WzWrFlMTEN579u3L\/k7\/cQVb+Ll0OJZouJ3KK\/HH3+crVu3LrkreLz9wGH1GQeH4JAQQFtqXg8aYoTLe2gkuMuWLWP33HMP6+rqSkTw8OHDiTDSb82aNey+++5jl156KVu+fHntnGR679SpU2z79u114ieKIqVfu3Yt6+7uTkLL4qiJGp5jx44l6d9+++3EBhJ6ElgxjXhofVrMzSlrnDdl4UhwWC1+wWG1+JJZCw7DcdhwgstFlUae6VtZenp6WH9\/P\/v1r3+diDGddkV39pKw0sh4YGAg+bcsFKk6\/ov\/vb29PRF5PtoVR9zpfHn+4z2srHJ0cBiuASj6JXBYFMHy04PDcBw0nOCKwpkluEuXLq1DOR0ipodievq3LDxM4WAaxdJIWvxR+JoWWonCnqYVgrsjmQrg0QhZJwUchmsMbL6kaqzhhzZolpMGHIbDfdwKLr9AIWv+NGuxDT\/A\/x\/+4R\/YP\/7jP9ZCzSJ1qpEzRrjnESjq6OAwXEOh+hI4LJ+DohaAw6IImqcfl4KbnsOlUdQTTzxRCzFz+MQ5XAod85A0hatVc7j0Hs0PL1myhN100011c8XiXC99Q1yQZU5Z47xZxNHBYRz1ABzGwUMRK8BhEfTypa2s4PIVkrJVynzOVBWO1K045hCmFzaJZ5emVzaLq5RNVkOP91XKPGRPq83BYT6njelt+GFMbNjZAg7tcLNJVVnBtSmsTRrs4bRBLa404DAuPmysAYc2qMWVBhxWdB9u6GrkenGTr4oXGpcqfQ8cVoktua3gEBymEahaW4oRbvXrMEoABIAAEAACFUAAglsBkmAiEAACQAAIVB8BCG71OUQJgAAQAAJAoAIIQHArQBJMBAJAAAgAgeojAMGtPocoARAAAkAACFQAAQhuBUiCiUAACAABIFB9BP4vsqJuUrdP2WIAAAAASUVORK5CYII=","height":0,"width":0}}
%---
%[output:911fa605]
%   data: {"dataType":"matrix","outputData":{"columns":3,"name":"C3","rows":3,"type":"double","value":[["0.6667","-0.3333","-0.3333"],["-0.3333","0.6667","-0.3333"],["-0.3333","-0.3333","0.6667"]]}}
%---
%[output:667910d9]
%   data: {"dataType":"textualVariable","outputData":{"name":"COS","value":"1x1x3 <a href=\"matlab:helpPopup PhasorArray\">PhasorArray<\/a> of <a href=\"matlab:helpPopup double\">double<\/a> representing a 1x1 <strong>real-valued<\/strong> periodic matrix with 1 harmonics"}}
%---
%[output:968e98ab]
%   data: {"dataType":"textualVariable","outputData":{"name":"Abc2Dq0","value":"3x3x3 <a href=\"matlab:helpPopup PhasorArray\">PhasorArray<\/a> of <a href=\"matlab:helpPopup double\">double<\/a> representing a 3x3 <strong>real-valued<\/strong> periodic matrix with 1 harmonics"}}
%---
%[output:18075dc5]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAdwAAAEfCAYAAADr33fvAAAAAXNSR0IArs4c6QAAIABJREFUeF7tfQ+wVsWVZ7syU+9tMQovJqOIiCJOrIrlGMrkiVJu1nWsUIVVmSKBRzKxXGolO5qZSh4b\/ugGmY38ccWpMDGRzBKLZOs9TEyyG6ac1cSdoQAlmgTYpIIRBQYR2JU8iPNqoFI4bp37pj\/663fv7XPvPX26+\/vOV5Uyj9v3nNO\/8+s+\/ed03wvefffdd5X8BAFBQBAQBAQBQcArAhdIwPWKrwgXBAQBQUAQEAQyBCTgChEEAUFAEBAEBAEGBCTgMoAsKgQBQUAQEAQEAQm4wgFBQBAQBAQBQYABAQm4DCCLCkFAEBAEBAFBQAKucEAQEAQEAUFAEGBAgCXgnjlzRi1fvlwNDAyo\/v5+hmqJik5GQPjUyd4NUzfhVBjcu02r94A7MjKiFi9erPbt26eGhoYk4HYbw4jrK3wiBlTEKeGUkIALAa8BF0aNGzduVPPnz1dLly5Vy5Ytk4DL5dkO1CN86kCnBq6ScCqwA7pMvdeAq7HUI0gJuF3GLk\/VFT55AraLxQqnutj5jFUPGnCvvvpqxqqKqhAIHDx4kFxtUecofCKHOjqBPvgElRRORedqNoN8cSqvAsEDrquyhw4dUldddZUTfOpy0Hm7bAOjqPVi5MVsG2Bi2oe11elgo0BZ5+jyGQZfH37F4hDCvpht4+CTK+AKp9pbJ4ajwqn8Hi36gFulI5ay8SGAbXhVLG8ScKvokbLxIeCDT00DbnwoiUVVEPDFqSRnuFWAk7LxIeCDzBJw4\/Mzl0U++CQBl8t7cerxxalgAbcIZkxFMcsXPpYAY9Ybs222LzA+pmqGGF0pYVeGC2U9KGX5bosYH1PxyV7OLpKbEn7CqTEETJ9xcoplSbnIye\/79BOq\/8MfVrdcM1ktu3N6bjEh83hYUsFk+OUT6nPr\/5s68fSDlH1goSxMw+HE7sjIWbXr9dMK\/vvLI2+pU7+d0LL9ir4e9fjAdePqwmmfVh5CZ93AjPExJdkw+lLCjzPg7j7Zq3a+dkq9MXK2jffT+nrVLTMmqVuvmRQkB6ZrA+6l87+k7rzzzqwj2vn66cwpUy6aoP7k5iuyAAxEhp9OmuqEv80EMLOh2nWzSZHC82Nvn1NHfvt7avil45k\/wZcnTpxQJ5\/4OGUfWBpwn3\/++Sj4ct\/wfgUDDvhN6+tR0yb3qMm\/ey77e+LEiVknpDkPzwduukx94tp3K\/O9k\/lktn3dHm6\/\/XZUMiMV4SDgxsIpjYHvPrEup3btO6D+Zv+o2vTSWF8OvH5vr1LXTpmU\/X3RBWfVT948q946o7JBKDz\/6DU9au3CG7Pnofo4Tk4FneHao0dwAnRS6589lDkDOqGimS9VgxI5zRGw\/QajV\/AdjGAxM4TmFoxJ4NSVZ\/PO106rR549lAXSW2dMUj+4b6wjcf3WP3tYDb98POuElt15lXC+BDBuH3Prc3Elxud2+9\/74M0oM4H3u147lbUX6O+\/svC6rM\/g\/nH6OKqAq4E2HQizpCc+db3TEd20rBNTXfUACXynB0ih9kcwDccXdtB5wEBx4KZLW4MNs+PA6IWA\/Zn\/\/nMFKwWuwIuRp2cNrmN1lLKwOuuWw\/iYssPG6EsJvzJsqtajLNBWkXXhxZdlbQf6kqKBahV5Lr7b3MP4mIpTUQZcM\/B+7bn92RIFdGR5e166rE+HUJLU5ThMPTBl6nZo2Lqajc0ODt0UcFds3ZPx0zVCr+Kzb796QWuV5wd\/emMm2\/5VkefqgChl+eYdZ+cIdcHoSwk\/bPt2lYMgeddX92TF8ma0dTCBPkXLtGe7deS56qDbBcbHrn4b+zzqgKsrASN\/cAR0PNjlCiwAUq46AmbDcC37c5KZU5cecMBSsAuD6giP7XHdP7w\/W26DgSYMOOWHC4CUOHFyitJun7LM1RwYbOcNCJvo1\/JdqzxNdJjvcvo4iYAL4JgdEIz6Q6z1Uzk4ZTmw7AMJQdg9Sk4yc+kCLsISGGQgQ6fgMxhydz6xc5PLxxoHbn2x46\/zDXwMMs2660kWtp9pghunj4MHXFcG4NGjR9WcOXMyPGFZ4aEfnVTb9o+qkcc+0pbVpjPcoFxZVrMtr6i8L3muDEDbnry66OUVV5bzjh072rAzZWs8QX5eXc3n+v8\/+sMjaviVf27tLxbpN+3jzADEZJQ29T\/sr9637WTGx8fnXaL++e0ThRibWZdN+ATHK2CQM++6ieqhf3dJxm9sRmfMfDJ9kcc3bXsoPoF+Dk7ZPmrK0aby8vooGGTC9t62V0bVtruntvWxvvoo0Ln4yT1ZVvP3P3l+haesD42dU8EDbp17SvWo31xq87nGXzZ6CqE3hM4qM61O3sP9wy+9mNFBb21w+cIe8VPqpZRlB07qtsM5G9EBt04flVfvlHC27dczW1hdfOc3x53321PVVa9sHnxrVP3iobGJV8qcij7gFoGrlzZhpis\/vwjoYFsHa84O0qcuc986VB4B2AABv44f\/DKET7pPH+fVglsfH5LVNAHvihL4qkmqV\/qux\/eoI6fOejk6xOnjZAMuuE0HXdnTrUdi7Ft9n\/875zGVIlmcZPalK4Zgq\/HVQZdjbwvLD85yvnwcA385ccToMnMVQgZbbauvoMvJqegDrmtpQgfdr\/\/xpWr+reOvyrOJ5ZKny8dcjss2ja3OFqyjl5PMGF1V66CDLdwUlXeRRVV5ro4OI69K0MXIw5QBu2Moh\/GxC+MqzzH6YsCFepkV5JnLyGYmMqa+mDJ1OKWDbtEqUx29GB9X4UxZ2eQDriYGZI1iZrp1HOKDzEUyMfZhytQhs22TPbOto5eTzBhdVergOmtIgXHdAeHTO\/ere793guR8ehVMXGd6fWOC8TFV5whyMPpSwg\/bl+lgm3f7E6a+mDJ1uQJBF35UA2CMj6k4FX3AxVYUnABnFjFBFyuzW8vpGRTFOThOMlPr0qPpGJbT8rioE6m66ZwutY9dbZxbn8sejudlwZZDP0YH7CkXrTph3jfLcPo4eMB1HQsCYLAfL7hjw4vq+D+ea2WzwSiryvsc5V3HgmJ4DmSG+5CX9vcUYm+PTm3s9PPYjgVh+KCvmtv+q5Nq08cuVbfcMDPjEQc\/MPaZ7eHNdyZnl8LAoODyC0819ldRW8P62\/f7nHzSM1zKPipGDpmc2\/jMz9U3fjbaOl9uzlSL2rhvn+fJ15OCWZf3qB8Ojt3dXNc+Tk4FD7jUKfdlyw11lzCwyzBc5Xwt1+gUfMgGpLqujXP0iNGFwQ5G+N968Y1k7vDWWeR5qzuY+mLKxNJ2MD6uOsMpK4\/RlxJ+ZXWFr\/3M23LUmSCJqS+mDAWnYNur6fFQjI+pONVxAdeVUMJFBNtBlHopZZmkh8sV4PYkyvt7OcmM0eXCTieKwQUTfzb3emc7c8nTAnyXK0omwejFlKHoHCnaBMbHTqdVKIDRlxJ+RVWHfnPuxpfVbX9wSemd9VgecGGit1X0YLOOXoyPK1CmtGj0AbdORbUTKPYg6+hP8Z2ijMSmdeEkM4WuJkegmmLV9H3X6k5T+TG8T+HjKvXg1lfFNqqyrpUtKj0+5UC7rZu\/w+njjgy44Fh75OPT2anL1rM6H8k3nGRuoksf\/4G967KvUsXu65QHDBhsm\/gYI98uw62vjo1N3\/E12G5qV5X3XceFymRx+jj6gFtniUCDC0ukEExgP1KfI2siL89pIeRR68R20nX0cpIZo6uoDnaDrVPXskbNJa\/OEhuXbTY+dfRifFylo3aVxeirU48YuAI26O9Zw\/EfM+muqX0hMIFkz\/f2qlYSFbYOGB+7eIJ9HjzgujIAm17kDWcVe3vGLi0AEjSVB8CaWXNV5bmykG35vi4GBz2QJAEX8QNJ4eIQUzf8f7PR2NjlPdfv6\/fgv5wZgNBw6vIJGit8jGDKRRPaPg5g+6Op\/6nl5fFJJ1GBTyGL0+SQTz6ZftcdUMoXzUMdmnBKtxGzXdi4AGahOKXzXZZ8aJJau3Csf3T1AXn2x8SpWX91uC3py+7DQvdRwQOuK0sZO3IoKudKomoqP+X3XUlSFHXjHD3W0eVzOZ0Cv7oymiyx1dXJ8V4dHzexi1tfE1urvNspWyh2netsJXL6uOMDLjikjhOqkDfFsuZSks9vC3OSuY4uWE6H79mmvG9bxD\/KywFi4XgdHzexnVtfE1ux73ZCkpSrrsB97EdGOH3sPeDu3r1bLVq0KMNn7dq1asGCBS2sMBWl2gswl9lSv3O5KSb2TVJN5dnkN+VhfOxqPOZzKj65Oh2fmJTVl1ovZn+eWqdPedR8Al9QccrFY5+4VOFUkzuSQQ+mHpgyWFl1yrky9n32UWW+8BpwR0ZG1ODgoFq5cmVmw5o1a9SGDRtUX19f9jem8VA6Dpzg87uKVUhfVBZTX0yZMpLqm6T0rK6pPK6AS8knV2amT0y4eAJ64Oagh350svTIREp1xfQZrsBnPqfklEtvDDiXXQdKaR+lrDoB13U9bUcGXBg5rl+\/Xm3evFn19vaq5cuXq4GBAdXf348OuC4SV33eictsVTDg3tuj7CAp+dTk3F4VvGMoy+1zn3Wm5JOe3cbWR\/nCr1P3bcvwwrRzak4Fm+FCBzk8PKzWrVuX2QABd\/bs2a1lZc6KahC6eT83xKXklD6m4JNr5Ourswst17XEFto+rH5KPumAG1sfhcWiSjnXFkoVWSmV1R+1MY+G2vZTcyrqgFv3GIdeZoD\/Qlq6TmnXf9vP9d9wjOeFf7pCwef88o5N1JE3Z86cDOO8FH9bXqhjQfpzbubtW76PcQCR4UeViY4JuC4+3T38D2rSxRe3jokV8aUKn6r4H8PPKvKwfMrbtze5qdtQ2UX0O3bsUEW2af4XtUXzudkWtby857Z91HzCBlwXp6oeDaTmAEYe9HnDLx9vO\/pW1IdiOYXxWdGxHB2UODj158+cztR9ee6k7L8m131wKmjAjXW5ppOW2cocDM9CLiVRjh6bLinrGT42e9GFa4rPMUtsMdeLkk864MbaR1H6wb7kn1J27LJcq1rUnAoWcCkSEnw6s1v2c0MOLijJ3IRPXMegfPKVQnZILlDYT8knsKcJpyjq41uGDjadevQNi1\/ZViI1p4IFXD2C1MeChoaGWglT8AxTUZ\/ZbmUfOfCpt8whGL2YMnq559uvXpAtJcHVbUXnbavIM5ebiurhMwPQPMJRhU\/6eMwnrn23dfNSUz\/YS2qpyLP3c336nxoTTJ+B7YR1ubqc0u\/HjF+okxkxYqKPhsIK1zu\/Od7qB3xwqoiDXo8FuYiPqahvx+mbhuwvTfjWiwlWTcrAu3rf1vVRAp91xfjYxRPs8zxdejl92uTz13tWHTRQB40Y5Jnnc336n7qunHyKYVLQFD99ZektN8x0NiNKHlDKohzY6sEm7OfqfoCTU9EHXCdLCAqkvsyWB0HIfVvTHk4y5+niuL6SgILsIlLN1ufkEzbgsjsPoVC2UPJByrvql5NTwQOuKwMQYCvKnCzLCtajIuz7EHTPnD2bZS771OfKAKR4ro8AmJd82FjZo0Zfz0N+vGDF1j1q2yuj6olPXZ8tp1PypSq\/fJVvyhe47D3v492++FCWBe3KaIXnnHzSATeWPqoKh\/Ql\/rCFgu0D82aSGJ809Sn3+2++M1nd9dU9Ge\/hC0mcnAoecKmOjCAGfc4imGvwnEIiKOC6RYnTRM7Ro6nLlZnIiUHMulI7n8vJpxRnuLGsbMXMebBN7+dC0P30H32Q7Oiiq97RB1zOvQBzmY3y25B5I8cix2DqW1bG1\/ct69aBs4M0ddnXV2q8MfjWrWtZY4tVL3TQcze+rK5+78TsbHLsdeDkEzbgxuTbvO2xEPaF0Fm13cL53COnzqq3v\/FJCbiukYKv5+bIx+dXdHzYH+OsjrODBF1\/\/5NfqvuH92cNqZvP21blVyqrO5x8wgbcqlj7LA8DzbITCT51pygbBii7f\/xj9f+++RkW86Of4bKgYClJMYkq1qUkzg4SdC352v\/OjkHBUtG0vp4Q9ElSZypJVJx8Sing6pUtuEkOztzKD4cA9znl6ANuqKWJOza8qHp7xo6SxLbMZmNSdE9qKOx8nsMt88WUuZ9T75v9iWBnjmPjCXaJTfsLMrqh47aPyMW0HB9jwI2hnZWtUISwL4ROLN\/tcpycCh5wXRmAIe8phUy\/WZf3qB8O3pz1OXlZrlXtc2WVgh4zay8vg0+T2ewot\/\/qpNr0sUuVed6O455Sl31cGYAQKD67eYf64l3vV38293oyf9n+sP+u6n9qedR80vtaj8+7RE25aELbJSHdxCc9yIDOOOY+6tjb59S8LUezWS2ct4+1j8rrP+Hfuo1TwQNuTFnK9gxFL7PFfC1aTBnJeTM8ztEjpy7cglWapWLeUuH0cQr5HLBnCz\/JV6jf1jg5JQEX4adYL\/4O8bk9BFxtRTjJzKmrKg4plddbFGCza0uFu17cPr50\/pfUb6fdUrjMzl1\/ra9bP7fnA29OTknARXhQz3RdVyQiRJEV0cF24KbL1LI7p5PJpRbESWZOXdQ4xSYv70ae0DaGyMIHTn1g8LtZ1ntM2b96FUKSA5uzkrPfiD7gxrT5njfTDWEf3KD0t6+dVa5gG8I2vVcT4p5STMOJAZOyLiKEfUU69UDz1hmTWt8PDnUP9YUXX5bdDnRi7\/PqxNMPNu9lkRI0p3SAywu6nD4zZ7awz859R7Ldvotg5MTEtKGOXky\/gaSLs1j0AddZA8YC+kMHI499hFFruyq9rxTShiqV5yQzp64qGKRc1g66IeqiZ7bAeW4f2\/pCnnM1g63MbOmYyMkpCbgV\/aaDLpx3417K1cE2hO6KMLWKc5KZU1ddPFJ9D1Z39EyXsw524iK3j219ZTNd37jEnMzmu+4+5XNyKnjAjTnlXjvZPvqy+2SvgjOLSz40Sc29fFTNmTMnK5qXkg\/\/Du\/bz+ylGf1clzefw8j2a8\/tzy7i\/8VDczJZWmbZxd\/dlnIf+xGOIj6ZPqx6zIj6WFARn\/QsEy4TgYzYMj7ncTmP72Zd855\/\/cen1aaXTmftbO3CGzOdXMfMtK\/yOHXv906on755VsGMG2yq6jO7jdt\/2\/Lg6M99205mJpnHtYr6DJe8Iv1N+ii7\/5M+Kn+IEDzguo4F1VmTLxsNUcnTI2\/znK4PvfoGKZCt94+o6qDt9SmPc\/SI0eWzrj78XySTsh5YWbv2HVDrd446E4iw8orKAefNm5PgizchcgIA+yJOmatNpn3UHIB+5v6t+zOx9jJyU5xtWynlUcqyB2PUGGP6jTKdVZ5FH3CrVIa7rLmn4iODUXc6UK9Uz9lxkplTFzfXYtKnb6TysbVhDjDzkgK5fVymz\/f+dopbSDHxFGsLJ6ck4GK9UlKOumHYI3zuvWICSFoiOMnMqYsSoxRl6VwGWGKmGmzqdgQyi5KCuH3s0ufjPKxr0JEiX2K22eVjStsl4BKhaS79uI7rFKk0Ay10OnXlEFWJRAwnmTl1kYDTAUIgkWfn66eVvjS\/zgcjzJUcF+e5fYzVZw4W6q5GQftf\/+yhbDk95tvtOoC2bVXA+pii3tEH3JT2AsAhZsODzuOWGZMUfOavbL9qLNAezxoa3F\/7JzdfUZoBnRImnGTG6EoJO+q9qiJ5TTHRg03gMQTc73\/y0rY7mIv06vPk8BzaCVwsk\/cz7cP4mKJj1DIw+rR99oDZDrxFOEO7H37peDZwAfz+\/Qcntu4D5+IA6GnKA9NWSlnUttnyMD6m4lTwgOvKUoaKFmVO2ll1Mf0NgfdbL76hIMMQftCQbvj9CWrixIkt37167HSW7Qg\/CLTz3j8xy8bUP5O02AxCV9Yq93POrFJMlnKqfNKdhG0\/tz\/LsuKBo4\/uHkt40pyGD9vfcs1kdWTkjBodHVWnfjshCyz6uTm4xPCdk09gYx1OQZvfeeJ3stmqicMVfT0ZBvADHOAHWEDbv+0PLin8+AA3Z2PjlItzTZ9zcip4wHVlKVONLELKgeALv12vnWozAxpg0ag+pL2UujlHj5y6KDHqNFkw04P\/7Xr9dBZo3xgZG1TCrynnuX1MoS+v\/QMO8INVMFgBk184BCh8jLU++oCb0tJEGeiU9aCUBTb7lMdJZowun3Xl8j+1z1LCBONjbOeHKYfRlxJ+XBxNCROMjzFcwZRhCbhnzpxRy5cvVwMDA6q\/v79lF6aiKTlOyDyGgO89N+ET7SAppTaG6TMwHZ9dRjglnKrDm6rveA+4IyMjavHixWrfvn1qaGiocsCtWiEpHxcC1B2k8Cku\/3JbQ80nsF84xe3FuPT54FRRDb0GXBg1bty4Uc2fP18tXbpULVu2TAJuXFzzbg0lmYVP3t0VvQJKPkFlhVPRu9y7gdScKjPYa8DVivUIMi\/gekdTFARHgDoxTvgU3KVBDaDmkznLlT4qqGuDKffBqbzKBA24wdAVxUkjUBRwk66UGB8UAeFUUPi7RjlpwD1w4IC655571LFjx9S8efPUunXrVG9vb2uPxB49dg3KUtFaCAifasEmL5UgIJwSeoREgDTgFlVERo8hXdx5uoVPnefT0DUSToX2QHfol4DbHX7uqFpK59hR7oyiMsKpKNzQ8UawBFwbxd27d6tFixZl\/7x27Vq1YMGCKIDWZ\/G2bduW2bNkyZIsszrG31NPPaVeeOGF1rJ9DDYWLdf5tk341Bxh4VM7hsIp4VRzBMZLYA+4MJIcHBxUK1euzKxZs2aN2rBhg+rr6\/NRv0oyodM5fPhwFmT1iHfhwoXRDAh0ZXRgmzVrVjQB1\/TrzJkz1fr169Vtt93WdgyskjOQhYVPSKBKigmf2sERTgmnmiOQL4E94MLIETrjzZs3ZwlVeTdQ+apsVblg5\/Tp06MKuDALf\/jhh9WMGTPU3r17owm44Nft27ezrwgIn6qyur288Gk8fsIp4VQzBIrfDhJwh4eHs0ABPwi4s2fPjiqogV32jM2XA6rKhVk4\/K688kqlcYSBS+gfdFJgz5EjR7JbxcwsdZ+2ab3Cp3ooC5\/yA670UfX4BG8JpyTgVmJP0b2qlYR4KAxLf1u2bFEPPPBAFtRiCrjQyLZu3dq2csExkEoh4AqfqjeGUHwCS4VT1f2l35A+qhy7IDPcmJeUY53Z6pHjihUr2jzKNZN0NUGzk4IZN\/gYfr6TzmJf\/hM+uZiT\/zwUn3TAlT6qnt9goCR9VEQz3JgTEsC21atXq1WrVkWRxFVGebtDqtc86N4y\/Tp16lS2rQLhE40PhU\/ncRROCadoEBgvhX2Gq0eQ+liQ\/QUhXxXFyIVR7aZNm9qKxnRsyTQstg7S9ivnkSrzCIfwCcP08WWET+2YCKfq8Uj6qMiWlJu7USQIAoKAICAICALpIRBkhpseTGKxICAICAKCgCDQDAEJuM3wk7cFAUFAEBAEBAEUAhJwUTBJIUFAEBAEBAFBoBkCEnCb4SdvCwKCgCAgCAgCKAQk4KJgkkKCgCAgCAgCgkAzBCTgNsNP3hYEBAFBQBAQBFAISMBFwSSFBAFBQBAQBASBZghIwG2Gn7wtCAgCgoAgIAigEJCAi4JJCgkCgoAgIAgIAs0QYAm4sX4tpRl08nYoBIRPoZDvXL3Cqc71bUw18x5w4SLwxYsXZ5+Ti+me25icILbgERA+4bGSkjgEhFM4nKRUcwS8BlwYNW7cuFHNnz9fLV26NPtUW39\/f3OrRUJXIiB86kq3e620cMorvCLcQsBrwNW69AhSAq7wjwIB4RMFiiLDREA4JXzgQCBowL366qs56ig6AiJw8OBBcu1FnaPwiRzq6AT64BNUUjgVnavZDPLFqbwKBA+4rsoeOnRIXXXVVU7wqctB5+2yDYyi1ouRF7NtgIlpH9ZWp4ONAmWdo8tnGHx9+BWLQwj7YraNg0+ugCucam+dGI4Kp\/J7tOgDbpWOWMrGhwC24VWxvEnAraJHysaHgA8+NQ248aEkFlVBwBenkpzhVgFOysaHgA8yS8CNz89cFvngkwRcLu\/FqccXp4IF3CKYMRXFLF\/4WAKMWW\/Mttm+wPiYqhlidKWEXRkulPWglOW7LWJ8TMUnezm7SG5K+AmnxhAwfcbJKZYlZQm4tHu9VRr47pO9audrp9QbI2fb3HDLNZPVsjunjyMfdYPkJHOeriMjZxX8b9frp9WRkTPq1WOnVW9PT1bNK\/p61OMD1+VWuQrGIXIMKO2jlCUBt7gFhcB5+OUT6n\/tfUOd+u2ElmHA+2l9veqWGZPUrddMav07pX2UsjqJU8ED7vPPP99KigInwU93YJ34t9k5m6S062qTDPv8wosvy4LL8EvH1c7XT6tpfWPB5YbfH2twEydOzP4LARieww\/KfPSaHnXvhye1JahR2Hf77bejks8oZiUQcE0+3f3XL6tt+0dbdXxvr1JTLprQwgCC70\/fHBuI5GGQAv988yk2+Zx8Al7YnIqdE7v2HVB\/s39UbXpprG0D3y\/7vQnq2iljgfWiC86qn7x5Vr11RmUDUfgt+dAktXbhjbmDb4o+oKg\/r9vHUXOSk1PBA64rA5CiI+4GGdB4YDS7\/tlDWfCA0evATZe1jWDzcLDfg3f0zJcCN+4Z7jef+5l65NlD2WDi1hmT1A\/uG+tIXL\/1zx5Wwy8fzzqhZXdeRYqBS7c8xyPAyScdcFPoo+x2vPfBm1GgAu93vXaqNTj\/ysLrnH0GSnBChTg5FX3ATWlpooxjlPWwZelAC\/rNgFlFJ8yMzYCd1\/CqyNOjUE4yT5n7OXX2\/XepgZsuLRxsuOqw87XT6v6t+9GB1yVPcyLmcjHbZs+COPmEDbgh8bPbbV6gxdgHAfs\/f\/fn2YqQa6CKkYcpY\/uWq\/8MySkJuAVejpkw2jZzVJs3K6tTB1MmBC5zj7OOPM4O8pLPfEf9j6V3lI7QsXVYsXVPtiwHqwU\/+NMbW0vzNl2w8mIuF7NtITvH2AMuLB\/ft+1kRsmyGW0V\/0IAv+urezKZRbPugAu3AAAgAElEQVRdjDxMGQm4zMsAnJ0xc9W8q4PAqBsG9TIwGA+z3fuG9ztHu66KcvqYWhdgfP\/w\/my5DQYeMACRXzgEwB\/9n1qhjj3zl2xGUHOKynBYCobtI+AkDLZ1roYP+UWJhVS6Qsvh9HH0M9zQzohRPyx7QrCFRobdq6lTDx1wjpw6WzjadcnlJLMvXbpzk71dl7f9PdcDzKNHj6qTT3zcnyJLsi9ONamAzjfwMdA27dL9jGuJuUldYniX08fBA64rSxka2Jw5czK\/lGUI6mdQrizLObQ8V4adbb9dl2+\/ekE2soXO\/xPXvltYV5CzY8eOQuzM5Zw87Mzndz2+J5vl\/c+PTWjJM5\/bvtHLSfBfzgxATEZpXf\/rGf+sy3vUDwfHElKgfnXlFfG5qrymfIL3TX+VZZRS8SkPuzI+6aXTc+fOqdGnB9XhX7zE1k\/75FRVDsCg42vP7VfbXhlVv3hoTqs\/pO7zTE4B9ut3jioYdH\/\/k+dXeEye5HEmdk6F6qOCB1xXBmBKewFlvQBFPfRMy0zj961Ty4flZQg6mKVVs66co0eMriZ+yBvxN5GX57sQ8kLotANsGY8\/8NAONWHChGw1B+NjymiM0ceFn57Z6r1VLr3mShf4AKMXU6YKB3zKw\/iYilPRB1yqiqYuRwfbkcc+EqwqOuhWsYGTzBy6oPP5wy+9qKpgEMxhiSvWWJvJexw+NmHj1lfmMuBdWQKfb3fDSleT7SXf9tWVz+ljCbh1vcT8Xt\/n\/y6K86E68EPDN2+pKYKDk8xcunQg6PS9LWaKt6nLC7ZQgMvH2hhufTbm+tQAnBEPGWy1XZ0YdDl9HH3A9bmUwLUc22TpRO8dmgk7oTFxBd1OXVI2+aIDgrmnmwKfYlnGLmsTOkEKLm+BLYxQfMIGeJ\/tUS8j5wVbn3rLuHzHhhezm6qojiKFvhpVAq7h7VCkikVv3sw2Btv08nLeTDdUB4lpOJTY6T1d+7xyCkHNtJESkyaDS3jXPO6mO\/RQfAodcO09W5tXIf3258+czpaXi4JuSNuqBnBMv0G12hP9DJeqoqnJ0TOomI+i6OzlsuVlTjJz6tJ80kEXk0yWGge57bWTc\/L0c\/uYW5+usyvYcvsmTx\/sKU+b3IO+PjUGm0NzKnjAdR0LAoA66WMGrmMc+jmQGT44sOU\/3JRxBEaMJhb2TCLUcz3S1UcG7PrFdizIB5\/efGdydi766398qYIlZk6+YvkUC1+KsIFgu\/jJPdlSJQzg3vnN8Vy+c\/JJz3C5+yi4XvGhH51sHf3zwdkyjmI5ZecyhOqDmn4cgZNTwQOuHAtqH3OZ595SWa6BmS789IcCQi0BYmYjvpa6quxrl430fdkXm057AFA2owvFpxBLylVWtkJwxfYb\/A3bXvYKTyy2FfE+FKck4BZ4JBRh4JNy+\/7vuaTu77UzSkOROWTABRrpDM66F8jndWaYDqNpMA3Fda1XfzQDc3cvxseUS5cYfVT46f1rc2WrqW85OKW3VcytJSpMdP19ysP4mIpT0QdcqoqmICeFfZsiHPMaHXaGQOUbzoZTZLM926eqWyfLqXLkjdvHXPow+9excwD8iD0uGFNduHwMdZaAG4nn9fGflJNv8pZVOcnMqauMNtDxYDKXI6FeMDPs4z8YQ7h9zKWv7PgPBpcYypSt8MRgX5ENXD5OIuD6XEqIablGj\/LN+5E5lhOpl5yg0R18a+yu1xhnuBx8Sm2JjQOTPC5jznPa\/OTsHLH8bYqfvaTeVJ6NNac8nbn85bmTWsmDMfWz2pZQ217BZ7iuDMCql7kDoCl9vAC+PwmfgDtz9myW5Wrbn+LF4NCR9vb0KGh0nBmA0BnHwic92wefzr917PIG\/aPkpyujNGY+fXTzP6jnPt\/f+rSc2QmaGa\/63+G\/nHzSAdcnpyBJct6Wo61b5KCOofu8ppya9VeHlb7vvcinui2E+iBGKE4FD7iuLOWYlyIobIMLJHa9fjqKa9so6gMyzEzLTf\/x3youH3PPflx4pbrE5qpX0+d6+6TOfh+3j33qq7Ok3hR7jveL8jk4dNfR4dPHtj0ScOt4iOgdV3YmkZogYnSjm7jzEXXkpb9lsYGz4WArJElU45GqkiRlv83tY1\/6OiFJytUGwM+QsQ\/f7Y7558vHeXX2HnB3796tFi1alOleu3atWrBgQcsOTEU59x9MgHzrLTpvh9GLKQN1CV3OdTa1TiNMkU\/Q8YT4pCKWA1w80UEG7IIz23X0YvqMqrwKwSnsueOyutTBj0se2AaX4mhfF+mNoQ4+OFVUX68Bd2RkRA0ODqqVK1dm+tesWaM2bNig+vr6sr8xFY3BIT5ICskF+nL2qoE+JUwwPsZ2kKnyqcoSWwjfcum0M3Hr6KXkE\/AuBKc0H4quba2Di48+qkmQhDpAfgr0c2XX08ZQV2pOlfnCa8CFkeP69evV5s2bVW9vr1q+fLkaGBhQ\/f396ICL7YxTKtdNe3uUZE6ZT93k86K2CJ1v0eUW2PZLySfQyc2pTt23LfNf7OdzqTkVNOAODw+rdevWZTZAwJ09e3ZrWZmzotgG7btcypdb1MGG0sfQOabMp27dz61yXaGLY5R80gGXi1PdsG+b5z+9tRTrfi41p6IMuPaxCTDSPAqg\/w6dIg92UB3jOPb2uXFHAGzZtr4UjwXpZSIgMvyospQxAdfnEY6m\/Mzzv+3vqnxveoQD3jePSJRdBF\/3CEfe2Uz7+Itu72V8p+YTNuBScWrF1j1q2yujatPHLlW33DAzaxtNOaU7d6o+yuRjkU\/q9FEw2NRHH03O1uWUrrd97MhsP6E4FTTgli0pv+\/TT6gp195Q+iFj14g3lefduJQEvqEcPXIv\/\/niVuxLbJT11sfeyj5WXkUfJZ90wOXa9sq75L9K3VMuS7nKQY0DNaeCBVxMQsIHBr+b2ae\/NEMNZizyunUPj5LMGD5RzaZ98kZzAc6hxn5kogkOPo69UfIJ6sbBKfvjHk0wTfndKsmDnPWk5lSwgKtHkPpY0NDQUCthSs9+vvncz7JviXZyxh4sn3\/rxTfUE5+6Xt16zaRSLmGy9jBl7OWUMqU+5VGT2TzCkccnV8D1WdcqGBft54awz5dO13nbOnqp+YTpo5pySg+w4JvR5nJqEV\/q4FKFexx6i+pgHxWMoa4+OFWEsdcsZdcoRVe07OaZGBzShMx6VEd5DjMlTDjJjNEVE3Z5ASmEfdQ64brC+7adzD32ZralOnoxPnb1O1WeY\/S56gF72LCa8c5vjnd9wAXszcGmCzvtK5\/lMD6uwpmgM9wy5WZFO3HJtVv3bU2fc5KZUxdFA4x1ia1p3Xy2ZW4fN9HX5ArLpj6I+X29xH7rjElRbCU28XFVnIPPcM0MQLj0Hn4\/HLw5+29eBh\/8e1lGXizPgVSLn9yjjv\/judaXc+xRmp1hZ9puLwnbWKTynPOyeWg4rozSWPhh8hsue4cZ0OUXnsp4X4Xfrixl7ufffvUCNfzycfXgv5mkZl3e0zajo+A7J5\/AF004pW8Yu\/fD57+ck0KfxsEZPdjUH\/fI6++5+jhOTgUPuPb+iGvfp+qIIlT5Tvi+JQV2nKNHTl0U2GgZPmeElHa6ZOlO1Oc3nbl9XEefXtmaNrknihmcy2+hnvu4+rVOXer4uI4eeCe6gGsvs\/lcuy8DrYnevOzMJvJsOyll2aNIakw4yYzRFSt2Ib8jSoVJlW861+Udxsd1O8O89zD6bPyKjkJR4aztjFke1rY630ROuY+KLuACmLGMfOo03JjPm9WpT9N3MB1WUx36fU5dVDabclJd3bE\/SuADm1A+rsMpiissfWIYm+zQN7DV8XFdDKMMuFCZFJfZJElqPA05ycypq26DK3sv1SQq6sstyjDi9nEVfTpJyueSug\/ehZYZOomqio+bYhVtwNVBF64D00lU1EsJlPJc96Ril1gw5TBl6i7ZUWICsjjJjNEVO3Z3\/\/XLatv+0SyJquzMNmU9msjKy1VoIi+Pf6Y8jI+bdorm+xh92j7XCoVPXKjbbVN5Vevq2v+vKs\/FgVCcCh5wXVmlkMGp08fLMvz0MwCa8l5RrDwY5W\/\/1cnCe1KxWca2vpTvUoY6c2YAYjJKq95VjPW\/HuAUla\/Cz\/\/6fyaqI6fOqsfnXaKmXDQhl88x8AkGBt\/42WjrC0C6E8urqz0AhDJV773l5pMeMLr6qJ+9+ob6Ly\/\/TnbueGn\/2MfWi9qtDgQh+igsR11Zyr76qKd37lf3fu9E22Cz0zgVPOC6bnHRI5+Bmy5VsFQT408ykou9gpkhUPmUUxeVzUVyYt9S0cunZd869YERt48x+mDPFn5U90X7wC0VmSHydzA+psIv+oCrKxrrxd\/d9rm9qsTjJDOnrqo4VC3PmYhUxzbXh8WrysSW5\/ZxmT5zG6nT78XG+oeiHAw2d75+2rmtQqFLr2K4Jn5UupIJuK41fipAqsjRwXbgpsvUsjunV3m1a8pydpCcujgcGDqZJK+OoS\/i5\/Zxmb5u+QgFB9dtHRrbryy8znn\/fFP7ODkVfcC1N8uLZro+N9XzHKqXPijvSAY9mHpgymBl+S7HSWaMrpSwA9\/ogaZ9DR5lPbCyfN6RXNZphkpwKZr92DPbbrojGdtfYDlVVs7cVqGQZ3IsFKeiD7h2Q9R7RyOPfaTpwKb2+zrYhrShtvHML2KCIJVJnLqobMbIKQq6mHepyugjb6H3Kbl9bOuTZWQqRuHl+D7XzMmp5AIuuClUwgbo1sGWO1kET8+4SnKSmVNXCJRhdSfEhe+hl5FNrLl9bOuLPZktBC996\/S9vMzJqeAB15VyX3SMw\/7sXZVjF3PmzMk4UuWYEXQ6X3tuv9r2ymjbxwiqHjMJlXJfdoykU45xdMqxoCJ+wpLuku+fUBMmTFDwbVX4+eYTfIxg\/bOHsu9Vz\/7Xb6gi28ylxk7hk15Shj7qwosvy77bfe7cuXFH\/6r2Aabf8vqg0PJ8c8qUv2PHDhSn4ArIn755Vv30s9NLOR\/7UbPgAdeVHVa2dm8utX157vkvcmD3g7Dl9HIalLc38X3uLRTZF0Kn2aG6RrSh9kcwI9WUsMvDWS9pHnxrVD3xqeudCSWY+uaVAT36TnC9moORVZcn2LaI8bGLn1Weg76\/\/8kvs2ALv7xs5BhwweLHVc4HJo\/uHuNk2epiHb2cnIo+4Loah7mn4iOjTXc6YEfo\/SsXFjE+5yQzp67QWMNFK67Op66N5gAztgx8bh9Pmfs5dfb9d5V28nVxlveqI6C3Eym3Vjg5lXzA1S6j3lvNG+FXp4e8wUlmTl0xeFZ3PtP6esatvNS1T7cjkBnj2VJuH1\/yme+o\/\/Sxm+TYX11CeXjPHBBSTII4OdUxARf8CkvM92\/dn7m47sjcDLTQ6dSV44FnSYrkJDOnrpicoS8KgNvYYLkNeFv1Z67kxMx5bh9z66vqt24ubw4O6wZe6O8\/dO+j6sTTD7JAGX3ArbMmbzoCOg+441RfBF+0XzUWaI9ny3TQYX30mh61duGNTifUsa\/pPkoInWBzHb2cHRZGV506NPVXXeyq6NWDTeAx8Fd3QK76Qlv51otvZIlY0E7Krk91ydL2+iyH8bGz0VYogNHns75VOFBUNoR9XDqB75DUp\/ttSCY0k7KKMNn4zM\/Vjw6dy260+lf\/dFKdfOLjFVhRv2jwgOvKUoaqlV30XfZ8xdY96ieQ2fbm2QwhuAh+1uU9auLEidnfo6Oj6tjb59qez3v\/xFagLcti1p1oVftcGYCd9jy2jxdU9Vfs5fP4ovd3wXYIvtMm96hbrpmsfnnkrYz7b4yczToa\/VzPaM1OsiyrPY\/79qDC1\/ucfII6YTLfY+dIVftS7IP0KZJNL53nNfB+8u+ea\/X3rx4bewbxANqF\/tgEJ6eCB1xXlnL9sUT7mzCSh9+u1061PbiiryfajyJQ1T2kHMwMgco+Tl1UNvuSAx0Q\/G\/X66fVkZEzWZDVv5Q5z+1jbn2++NBtcqG\/z+M94AADTPPTl5w+jj7gci1N2ISMWW\/MttkzHU4yY3SlhF1ZJ0lZD0pZtv+p64DxMWVwwehLCT9qfxTJSwkTjI+pOMUScM+cOaOWL1+uBgYGVH9\/f8t2TEVTcpyQeQwB02cYH1cls\/Cp3n66dI7FTBNOCaeq9kN1ynsPuCMjI2rx4sVq3759amhoqHLArVMpeSceBKgDrvApHt+GsISaT1AH4VQIT8aj0wenimrnNeDCqHHjxo1q\/vz5aunSpWrZsmUScOPhGYsllGQWPrG4LGollHyCigqnonY3i3HUnCoz2mvA1Yr1CDIv4LIgKkqCIkCdGCd8CurO4Mqp+WTOcqWPCu7eIAb44FReRYIG3CDIitLkESgKuMlXTCoQDAHhVDDou0oxacA9cOCAuueee9SxY8fUvHnz1Lp161Rvb29rj8QePXYV0lLZyggInypDJi84EBBOCUVCIkAacIsqIqPHkC7uPN3Cp87zaegaCadCe6A79EvA7Q4\/d1QtpXPsKHdGURnhVBRu6HgjWAKujeLu3bvVokWLsn9eu3atWrBgQRRA67N427Zty+xZsmRJllkd4++pp55SL7zwQmvZPgYbi5brfNsmfGqOsPCpHUPhlHCqOQLjJbAHXBhJDg4OqpUrV2bWrFmzRm3YsEH19fX5qF8lmdDpHD58OAuyesS7cOHCaAYEujI6sM2aNSuagGv6debMmWr9+vXqtttuazsGVskZyMLCJyRQJcWET+3gCKeEU80RyJfAHnBh5Aid8ebNm7OEqrwbqHxVtqpcsHP69OlRBVyYhT\/88MNqxowZau\/evdEEXPDr9u3b2VcEhE9VWd1eXvg0Hj\/hlHCqGQLFbwcJuMPDw1mggB8E3NmzZ0cV1MAue8bmywFV5cIsHH5XXnml0jjCwCX0DzopsOfIkSPZrWJmlrpP27Re4VM9lIVP+QFX+qh6fIK3hFMScCuxp+he1UpCPBSGpb8tW7aoBx54IAtqMQVcaGRbt25tW7ngGEilEHCFT9UbQyg+gaXCqer+0m9IH1WOXZAZbsxLyrHObPXIccWKFW0e5ZpJupqg2UnBjBt8DD\/fSWexL\/8Jn1zMyX8eik864EofVc9vMFCSPiqiGW7MCQlg2+rVq9WqVauiSOIqo7zdIdVrHnRvmX6dOnUq21aB8InGh8Kn8zgKp4RTNAiMl8I+w9UjSH0syP6CkK+KYuTCqHbTpk1tRWM6tmQaFlsHafuV80iVeYRD+IRh+vgywqd2TIRT9XgkfVRkS8rN3SgSBAFBQBAQBASB9BAIMsNNDyaxWBAQBAQBQUAQaIaABNxm+MnbgoAgIAgIAoIACgEJuCiYpJAgIAgIAoKAINAMAQm4zfCTtwUBQUAQEAQEARQCEnBRMEkhQUAQEAQEAUGgGQIScJvhJ28LAoKAICAICAIoBCTgomCSQoKAICAICAKCQDMEJOA2w0\/eFgQEAUFAEBAEUAh4DbgpfdAdhZYUCoqA8Cko\/B2pXDjVkW6NtlJeA24qH3SP1jtiWBsCwichBDUCwilqREVeGQJeA66tOMYPugs90kVA+JSu72K1XDgVq2c6wy62gBvzZ8o6w5XdVQvhU3f5m6O2wikOlLtbB0vALfoA99VXX93d6HdB7Q8ePEheS+ETOaTJCPTBJ6i8cCoZCpAb6otTeYZ6D7hlo0YIuK7KHjp0SF111VVOkKnLYWwDo6j1YuTFbBtgYtqHtdXp4H8pIHwaAwLDE2w5rI8odWJt880nkC+cEk654hC2f3KV8xpwXR90xzZ0VyXkebwIUPpY+BSvn7kso+STDrarV69Wq1atUn19feOqQa2PCyfRg0eA08deA67rg+6cFcXDLyUpEaD0sfCJ0jNpyqLkEyAgnEqTB5RWU3OqzDavAdcFCqaioZaxYtYbs232UiHGxy6eYJ9jdKWEXVm9KetBKavKUnEdvRgfY\/mCKYfRV6ceXL717Y+ieqSECcbHGK5gygQPuM8\/\/3xrjxacBD+9Zwt\/Hz16VM2ZMyf797znurx+Zr8fmzxzP9okZZH9JhZmXeHf89435e\/YsaMQO7Mh5um2G6rti7znefbdfvvtzn16DFExZaDhCJ\/G8h2ETxjGuMsIp87n0Ain3HxxlQgecLk2q11AyHM\/CHCOHjl1+UFLpLoQ4PYxtz5X\/eU5PQKcPpaAS+8\/kWggwElmTl3i5DAIcPuYW18YVLtbK6ePow+4Ke0FlNGWsh6UsuylYeo6cJIZoysl7Kh9USQvJUwwPqYMHxh9KeEnnBpDwPQZxsdUnJKAW4BkzI0oZttCkhnTcFLCTjrHsJ0jaBdOjWchpg1hyvge7GPbD8bHHRNwXUkuUNG8xCHtrNSeu5KmOu15bElTqfHFZW+n8aWorev2zsknHXClj2pPxOs0znFyKvoZLtXIQuSEQYBz9MipKwyaopXbx9z6xMP8CHD6OPqAm9LSBHYJo2m5lDDhJDNGV0rYNeUJdskuJUwwPqbssjH6UsJPOBV2m0ICruzhkt7LK3u4+YSKuVOO2baQfJI93PpcFk7lYxd9wKUcrYosfgQwMwQqqzh1UdkscqohwO1jbn3V0JDSFAhw+lgCLoXHREYhApxk5tQlLg+DALePufWFQbW7tXL6OHjAdWUAytWO4zME9XKNXO3Y3lFAwxE+ydWOlOFDOCVXO1LyKXjAdV3tmNJegCQkhE1IwIxUhU\/jWZoSJhgfU3aQGH0p4Sd9VNg+KvqAS9l4RBY\/ApgOi8oqTl1UNoucaghw+5hbXzU0pDQFApw+loBL4TGRIXu4wgEWBDg7R6gQtz4WEEVJGwKcPo4+4MpyjSwBYvsHTMMRPgmfsHzCBlzhlHAKy6ngAVeSXCTJBUtWVzlJcJEEFxdHqj4XTgmnqnKmrHzwgOtKmqKsrMjiRwAz66SyilMXlc0ipxoC3D7m1lcNDSlNgQCnj1kC7pkzZ9Ty5cvVwMCA6u\/vb2HEWVEKx4iM6gj48LHwqbofOuUNH3wCbIRTncKQ6vXwxak8S7wH3JGREbV48WK1b98+NTQ0VDngyv6I7I+YCAifxtCgbBeUsqhts+X56ByFU8IprpVWrwEXRo0bN25U8+fPV0uXLlXLli2TgOsYgGE6P0wZ3x1fWTV8fdxZ+HQedUoOUMryzTvqgCucEk5Rc6qsb\/QacLViPYLMC7iupCmQId\/DHZvRmFjYHVuMz4HI8KMePQqfzifaddq3Scu+h+uLT8BR4ZRwqvpidPU3ggdc6s64OgTyhk8EfIweyzpH4ZNPb4aX7YNProArnArvd58W+OJUns3RB9yUlruwy6xNy6WEiQ8yNwm4KWHXlCfY5d2UMPHBp6YBNyX8hFNjCPja9nINDCTgFiAUcyOK2TYOMkvAlaQpV8dW9blwSjhVlTN1ykcfcOtUSt6JBwEfM5ImnWM8yIgldRDwwaemM9w69ZB34kHAF6eCLSkXQctZ0Xjc212WcPqYU1d3eTGe2nL7mFtfPEh3jyWcPmaZ4ZYFXFeWsnwPV76Hi236cg2fXMOH5Qq2nHBKOIXlCqZc8IArGYAYN6VbhnP0yKkrXY+kbTm3j7n1pe2dNK3n9LEE3DQ5kozVnGTm1JWMAzrMUG4fc+vrMHclUR1OH0cfcFPKyC1jF2U9KGWBzT7lcZIZo8tnXbn8T+2zlDDB+Jiyl8foSwk\/Lo6mhAnGx1SckoBbgGTMhInZNjsYcJIZoysl7KRzHEMg1JlJ0C2cGs9CTBvClKEeONaVh\/FxRwTcKXM\/p+7+9KfVpMmTs\/qcPnUq+28n\/63rputbVNdOeb7lm99Ux575Syq+lsoRPp0qbDvCp3oUFE4Jp+oxJ\/+toDPcSz7zHTV16lTK+oisyBCALPOTT3ycxSrhEwvMQZVw8gkqKpwK6m4W5ZycChpwOafyLJ4TJeMQ4PQxpy5xdRgEuH3MrS8Mqt2tldPH0QfclPYCymhLWQ9KWXX3PbB15SQzRldK2GExblouJUwwPqYMHxh9KeHXlCvY\/iIlTDA+puJU9AGXqqIiJwwCnGTm1BUGTdHK7WNufeJhfgQ4fRw84LpumgL45Xu4aX4PF0a5t99+O\/n3cIuaJDQc4dP4m8k0XuasI8bvJ2Ps4+QT4CacGn\/TVKd9g5mTU8EDruumqZSWJmS5ZgyBUMc4MCNV4dN4lqaECcbHlHMkjL6U8JM+KmwfJQG3gIExN6KYbZOAm0+omH0Ws20h+aRnuDIpaOc0hi+YMtj9YN\/lMIMqqkFc9AGXqqIiJwwCnGTm1BUGTdHK7WNufeJhfgQ4fSwBl9+\/XaWRk8ycurrKiRFVltvH3PoigrprTOH0sQTcrqFVmIpykplTVxg0RSu3j7n1iYf5EeD0cfCA68oqle\/hyvdwsU0Qk1EqfBI+Yfmk93Clj2rnDOBinhzJO0Wi93Bdmec7duxQc+bMyVxSljmvn5m67b1dKGO277zn+n3TPslSrtIipGzUCHCOHjl1RQ16BxvH7WNufR3sumirxulj7zPc3bt3q0WLFmVgr127Vi1YsKAFPGdFo\/V2hxtG7WPhU4cTxlE9aj6BOuGUcMqViU6FkNeAOzIyogYHB9XKlSsze9esWaM2bNig+vr6sr8xjSel9PIyp1DWg1KWvexCXQeMj7FkFj6dR4qSA5SyUuIT2CqcEk5R9lGuvsxrwIWR4\/r169XmzZtVb2+vWr58uRoYGFD9\/f0ScAs8g+n8MGV8d3zYwExJZuGTdI6UfNKzW+mjxnhF2a9QyqK2zZZHzamyvtF7wB0eHlbr1q3LbICAO3v27NayMibJBd6Tqx3TvNoR\/As\/quUaCLjCp\/yEFbsTSfXqxqK2rlfEKPmkA65wSjhF1UcFn+G6yMxVURcQ8twPApSjR0zAFT758WMsUin5hA24wqlYvO\/HDmpOBZ3hNl2u8QOxSOVCgJLMFEvKXPUWPX4QoOQT1ZKyn7b+le8AAAb3SURBVJqKVC4EqDkVLOBKQsJ56Cn3NChl2UuRZWSpo5eSzMIn4RMlnwBN4ZRwippTwQKuHkHqY0FDQ0OthCm9JyPLNVzjuDB6qMlsHuEQPoXxaUit1HySPiqkN+PQ7YNTRTXzmjTlgpOzoi5b5LkfBDh9zKnLD1oi1YUAt4+59bnqL8\/pEeD0cfCAK9emybVpVE0Ik\/UuVzvK1Y5V+CacGv8BesBPrnaswqLzZYMHXNeScp19Q+p9yNjkpYQJ5+gRoysl7Lh4lxImGB\/X6wrz38LoSwk\/4dQYAqbPMD6m4lT0AZeqoiInDAKcZObUFQZN0crtY2594mF+BDh9LAGX379dpZGTzJy6usqJEVWW28fc+iKCumtM4fSxBNyuoVWYinKSmVNXGDRFK7ePufWJh\/kR4PRx8IDLD69o5EbAtU9PZQ80HPl1PgJcfAIkhVOdzyeoIRenggbc7nCl1FIQEAQEAUFAEFBKAq6wQBAQBAQBQUAQYEBAAi4DyKJCEBAEBAFBQBCQgCscEAQEAUFAEBAEGBAIEnDN+3DXrl3b+j4uQ31LVZw5cyb7Zu+2bduyckuWLFHLli0LbVau\/qeeekq98MIL2beGe3t7o7DxwIED6p577lHHjh1T8+bNY7NN+NTc\/cKndgyFU8Kp5giMl8AecF1f5\/BRSaxM6HQOHz6cBVmwc\/HixWrhwoXRDAh0PXRgmzVrFltQc2Fo+nXmzJkKPst42223tX2swiWjznPhUx3U2t8RPrXjIZwSTjVHIF8Ce8B1fdPUV0XryIWgMX369KgCLszCH374YTVjxgy1d+\/eaAIu+HX79u3sKwLCpzrMPv+O8Gk8fsIp4VQzBIrfDhJwh4eHs0ABP1jCnT17dlRBDeyyZ2y+HFBVLszC4XfllVcqjWMMS8rQSYE9R44cUfv27WNbUtZ6hU9VmTRWXviUH3Clj6rHJ+FUOW4ScHPw0Xu5AwMD3pdEq9Aalv62bNmiHnjggSyoxRRwoePeunWr2rx5c7anzDWQSiHgCp+qsPz8QCAEn0C7cKq6v\/Qb0kdFGHBhqdbsmGMKbLHObPXIccWKFW0e5UxOKqOS2UlBwAUfw8930lnsy3\/Cp3qddyg+6YArfVQ9v8HAW\/qoiJaUY05IANtWr16tVq1apfr6+uoxjuktu0NiUluoxvTr1KlT2Wa4wicazwufzuMonBJO0SAwXgr7krIeQS5atCizZmhoKJplWxjVbtq0qQ2lmI4tmYbF1kHafuU8UmUe4RA+1esqhE\/tuAmn6vFI+qjIlpSbu1EkCAKCgCAgCAgC6SEQZIabHkxisSAgCAgCgoAg0AwBCbjN8JO3BQFBQBAQBAQBFAIScA2YzIsufF96QX0TU6iLJ1As69JCwqcudbzHagunPILLIFoCbkHA9Ym9r+BIHcR9YtANsn0P2jSGwqduYNNYHYVTaftaAu6\/+M88PwaZyXCnMlzr+MEPflCtWbNGwb3Fjz32mJoyZYp69NFHM+LD5RNmNq55ef8NN9yQnTW2jxfpq\/TuvvtuBXcOm++A7CeffDL7d924dNa0mS2t73kG\/aYe89B5DLdPpd00mlkvfGqGn7w9HgHhVPqskIDrWFKGgAtfwPnsZz+bXT8JgRa+JgSBEX5f+MIX1COPPKLe8573ZB87gIse+vv7s3Lw1Rz7az5mUIT34UYmffGHOVMxv95y9OjRzAYI9BBgzXfMDy7YwTx9eqZdg7zlP+FT2j4Nbb1wKrQHmumXgIsIuDqowszT\/qLQ4OCgWrlypfr1r3+dBVk9q4XACjPjDRs2tM1yzaBadOWf\/nfzjmnd0PSM25arqyHLys0aBOXbRZ2j8IkS5e6SJZxK298ScBEB1wycZQFXX+ahRdpLxPDv5vvwd97ysL6LWH+XV8uD5Wv45J0Z2G36ScCNp0EWdY7Cp3h8lJolwqnUPNZurwRcwoCL+ZhAWYKLvu3ni1\/8ovqLv\/iL1lKz6bKimbPMcONriE07R+FTfD4NbZFwKrQHmumXgEsUcO09XPPrOWbilLmHC0vHekkalquL9nChHOwPL1y4UN1xxx1te8XmXi9UBb6VqxOymlFD3m6KQJPOUfjUFP3OfF84lbZfJeAa\/tNZgHlZynrPtGhJ2ZVxrNXYiU3mna12ZrN5tzMmG1qylONqjMKnuPzRCdYIp9L2ogTcAP6Tc5MBQO9glcKnDnZuoKoJp\/wALwHXD65OqdTJTb4aiLMiUiAKBIRPUbiho4wQTtG7UwIuPaYiURAQBAQBQUAQGIeABFwhhSAgCAgCgoAgwICABFwGkEWFICAICAKCgCAgAVc4IAgIAoKAICAIMCAgAZcBZFEhCAgCgoAgIAhIwBUOCAKCgCAgCAgCDAj8f7cGvSWoTYPxAAAAAElFTkSuQmCC","height":0,"width":0}}
%---
%[output:39d73856]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAdwAAAEfCAYAAADr33fvAAAAAXNSR0IArs4c6QAAIABJREFUeF7tfX2QFsWd\/9dgKauEww3lKSBZXyB6P6+ICBERC6+MMUpBJYQLLNS5x20QkxAuxxLepHgreXe9khgNwT0KPVkwOeKPLfRnTC5ngcipgVBYWgHFBZHzEkCOogCvuPNX31n7sXd2Zrp7uqefeeb5zD+7z0y\/fr6f6c9097e7L\/jkk08+IVxAAAgAASAABIBApghcAMHNFF8kDgSAABAAAkAgQACCCyIAASAABIAAEPCAAATXA8jIAggAASAABIAABBccAAJAAAgAASDgAQEIrgeQkQUQAAJAAAgAAQguOAAEgAAQAAJAwAMCEFwPICMLIFBkBHbu3EmXXnopDRo0iE6fPk1PP\/00HT58mO6\/\/366+uqri1x11C1DBIrIKwhuhoRB0kCgyAicP3+etm7dStu2baMpU6bQsGHD6JlnnqEhQ4bQ5ZdfTtu3b6dRo0ZRt27digwD6uYYgSLzCoLrmCxIDghUKgK7du2i1tZWWrFiBdXU1ATV2Lx5M82dOzf4f+PGjYGoiosbxjNnztBbb70V3OIe7oYNG+iiiy6iQ4cOBT3cvn37ViocKLcjBMCrz4CE4DoiFZIBApWMgBDW0aNHlwT3wIEDtGzZMmpubqb9+\/cHYrxkyZLg\/48\/\/piuv\/566t27N3GDKgR35cqVdN9991HPnj2D8Cy6F198cSVDg7JbIABedQYPgmtBJkQFApWCwIkTJ6ipqYnmzZtH\/fr1ozlz5tDw4cNp\/PjxgWByj5QvnjcTPVxuLMXvs2fPBvE5Hg8X846wPXr0CHqzQnAHDx5Mjz\/+OE2aNIkuueSSoHdcX19f6i1XClYopz4C4JU+VhwSgmuGF0IDgYpFgHuskydPpqNHj3YZHuZKyQLLQ8r8u729nWbPnk3csDY2Ngb\/y8PKHE8ILt\/fs2cPPffcc4ETFfeAx4wZU7F4oeB6CIBXejhBcPVxQkggUAgEeMiXBVeepxUVSyu4hQAGlbBCALzSgw89XD2cEAoIVDwCLKibNm0KPIr37dsX9FblK0pww0PKPCQ9YMCAiscCFXCHAHiljyUEVx8rhAQCFYsADwm3tbVRQ0NDUAduJPniOdy4Hm6U01RUz7hiQUHBrREAr8wghOCa4YXQQKCwCIR7uEKYeVlQnz59aP369ejdFtb62VUMvPoMWwhudjxDykAACACBikWAPdPZK51HRviaOnVql2mIiq1cmQoOwS0T8MgWCAABIJBnBKK81CdMmNBpGiLP5c9j2VIJrvjy4TV24SUCXEleJjBx4sSgvsuXLy8ZKG7XGhf35a8x3vGmpaWFamtryfS+WP6wd+9ekjcBcHWfMVHhl0eioExAAAhUNwLsiVxXVwfBtaCBseDKwhPe6o3LIS+E5t9ip5rjx4932bWGHTCOHDni5D7v6SrWDMrEkL\/SdO6LMLx+kIdTxEeFq\/sq\/CxsiahVioDp0N\/hE+eo9fUPafbddZkixvlMa32btn7\/pkzzQeLZIyC36\/BST4+3keDyi71mzRoaN24czZw5M3YRPIsT9zB58bwQLd7JJmqJwe7du63vz5gxg5588slOO+fwtnILFiwItqKTd9RJuj9r1iyaNm1aqV5CrHkZhbzoP+396dOnJ+J3zTXXpLdkjmIePHgwR6UpflFUQ39hXp38Rgt97swx+t9LelOPHavowmN\/yAQkOZ+ev5od5GlzgVc26KWPGzciV5T26je\/+Y23U62MBFeYTLXrjNgAncOLLeT4\/6hda1iIbe\/\/8Ic\/pC1btpR6ozykzaL\/2GOP0apVq7Tvr169mpYuXRpsf8dfccK7joV4\/vz51vfFkoo4\/JjASY3Ke++9l0iMrJ+zDXXyuPPOOxPrkf7VR0wdBMJDfzKvuGe78sX3aOv3bgp6n3w9em8v57wS+fx+\/q1046LtNPJLvekn9TfEFh+80rGs\/zBJPVu0V+b2gOAOGxbMOXMjlRfBDTc+\/FtcfL6ozXM5HSGg\/Fekm\/Scn\/GRa7wXrzjnVJRN\/ObnvNYTvRHzl9FFjKgGUm4Yx\/xkD11V270kfrUzfkttDf3otkHxm1noiGH43FvO57brLguGrTn+zT9upxOP\/BUE14WRPaXBXFq8eDEtXLgw8IcJXxBcc0NkIrgYUo4eghY7+6Tt4ZqbtzwxVC9ieUpV\/FyThv74A0gIH\/duR1zXK\/g9esMRenDUAKofekXsh5z8cSULq\/yxJf\/Pc7dffuhV+tnYK2jciI5eLQv7oq\/2pun3\/mXwO\/yhpvMbIyd+Oczt+Nq1aztlKjvBFuU991kP54ILp6kOz+04JytmLwTXb8NRDbnpDv2xEPIwr7i+\/+mwctJwryl+K19sp9bX\/8N5Pj4bRtM6V2N4\/ojqX9u9k50rDQfmavMz\/4\/++NQDXoruRHDDQw\/ysiDZk1ks\/wnvWuPivuypKS\/nMb0vexHLC71d3bcR3DRDezKLbOOLnkl4+DCcB3oiXt7dUia6Q388zMuX7DXMc60surbDvTInwsPWzLtdx2oS89HhJnjll1eq3Opu\/Ap9t\/nZwCcgij86NlW1JTbPVe0Vi+3Tr75PJ17bQkef\/0dVdZ08TyW4TnJGIl0QKMoXfFHqUSkU1R364x7J7Luv7rIciO+LYWYXdY5KT0fYVXmDVyqE\/D4X9ljz\/D5a9OtjJQ5FTQ\/ETUVEiWLcVIUIy3\/jfEh0piZEfB7t+eWkK8jnhxwE1y9HE3MTBI5zitKZSysngUXePgmcI\/PltiiCV3HCyg1P\/dArna3LjcuH7yf1pFUAQnBVCPl9LtuDe4vC+539A\/J87XjnJI15fE\/gOMi+Cz55BcHNETN8Gj7LahelHlli5DNttsdTv9pND\/zzPnpz0e2lrEVv4OFd54J7M4d177Q8yKSnIsI+u\/+CYP6Wew7hngg7aInlQSY9EXzI+WSLfl7h9zxqykI\/NT8h2aGPxfa2a3uVPPV9tlcQXD921spFZfi8z4mgYdQys\/dAzKupT\/wrvfLOR5G7PokhwbjepwnvohpdEZ\/nit8\/cS6yDDp5YOTEO3W0RuTkQPLoho5NbeZoVelHjfaFHfo4DZ+8guDmiMMQ3BwZo0BFYV7d2PQvpXWx4aq5Flyx\/lbkIxrGpHlcVePpu2EskPkzq0pUe8UfVa+8ezLwXNaxqW\/B5Q8CMZRcjg4CBDczOponrBJc8xTLE6Mo9SgPeu5zZXvwNotyQxPOxXZ+VaSX5IAl1uemnccFr9xzwybFOHuERc0mD5dx5Y8BOV2fvILgurSoZVpwmrIEENEjEeDlG6e+tpJ+94O62Dladpy657rutHzCZwcNmM7h\/u6Dc3T\/lg8Dx6i4OVrecYqFf1jvs0FZTbxNfQ79gUpqBJKEKrzeW51a9iHiPgQguNljn8scVIbP4xBNGEgM\/eWPWv2\/cg+dHjEr1kOYbSYcp6I2wNDlXdSGF2LYTggrN8Syw4pASycPCG6+uJXUQeAPK7HTWNyqC\/HBlfRc1Djqw0x81CV9uImtaHkdOHtRh535fG9Fix5ujjisEtwcFTWxKEWpR6XgrSpnn3v\/gS4f\/u3EHYHihttUacvPk5yiRDibna3AKxNrZB82yR46XMi+hJ\/lkDTM7ZNXEFyfVlfk5dPwWVa7KPXIEiOfaV9+309p2C23JJ5L62JjivAOU1F1tMkHvPLJGnVeSfYQa11dbqiiLlF0CPnkqqgQPnkFwU1rxQzi+TR8BsUvJVmUemSJkc+0WXCbJn1dubEF9wLYu5T3x01zxe1kFU4rbT7gVRqrZBdHZQ+dD7DsSvdZyqpyqOrhsowQXJdoWqalcpri5HE8nyXIVRZd9gxWzaVFLZkQcKl4d\/TU+eDkIeGBHLU9n+Av58Pzezf3\/WyjDRz7WHnEVAmVy+VmUeio5v05zi92vB048sX1tH37nEBwc8RzFYFVBMv6OUOlkwecW\/yTSj4wRD5CTWcpjrBp3BaPOjb\/4H8uC3bw0dnEPmo+TScP8Mo\/r5Jy1GmvvvnMh8H+3byFYvjSsbnNOl3Or2Hd6\/TRf18YO50Cwc0Xp7yWRkVgr4WxyKwo9bCAwGvUuCMx+dBwMZems\/bVxqEpzkM5Cog4T2UVaOCVCiG\/z3XsUW7nKZ2DOXTq4QpZ9HBdIekgHZ+Gd1Dc2CSKUo8sMXKZNvdu+cSglpYWqqmpoTlz5lB9\/afnMr\/YTqt\/+Tod++lfK7O0aRxZcOO2jgxnnFbYwSulCb0GUE2Bid5pWPRM13dzpUzWa8vhb1y0ndoa+iXG9zlyAsH1StHkzFQNStZDMKr0MaScI7JIRWHBbW1tpRUrVgR3WXCHDx9O48ePJ3YY2fXv\/554wLZs96gdp1S84Od\/\/\/zJIG\/5rF1RxHB8cbKM3OvWycNnw5hPS+erVLrtVZzTko7NbYaUOd\/LLjpPG6YMjQUOQ8r54pTX0uh+McY5v+gc3xe3A1DaL8io9NAweqUNJQkul8SEV7xhgfBUNumJcE9i9PU9SjtVqXjGwi7vfKUK77th9GvBysxNJbhyrVxtHWqClM5wsvx+mKSdNix6uGmRs4h39uzZTsN+IikTAltkn3nUotQjc6AcZZA0pGzaoKTZB1c4ZpmsuUyzNAi8ckQYg2TinPF88MqgmF2CmvgU+OQVBNfGqinisoNLY2Mj7d27lzZu3EjDhg0rpeLT8CmKrh2lKPXQrnCZAyY5TZk2jOzQFOdVGlfNNJscpBF28Mov0VzyKu28fdoa83By+NSquLR88gqCm9aiKeJxz3bNmjU0btw4mjlzJs2ePdtIcLOe81Clz1VWhcHQXwpiOIgi90RMP+Rkm0Z5EKtsLtY6mpynG16CpMoDvHJAEsMkbEdOZJuKHcbkjVV0bJ5mDpc\/AKdtelv7iECfU2AQXEMSugguerlhwb1i3EP0rW99i06fPk09evQoZcW\/eWMBvmq6d6ez584Ff8XFv8Ulng\/s06tT\/HC5RfqcNl\/8O5wP3xdpi\/xO\/td\/0cUXX1zKn50SRHz++5\/\/+Z+0Y8cO+vAX811AhTQcIKD6gpcbvqieiKphnLtpD73wzrnYvZqj4ofzUeUBwXVABMMkdH0D4pIN2zQ8eqJj8zSCy71bvtiBTycPCK4hMSoteJzgij1vo+pzlWK7vf61NXT4RMeRZ3y9f+IzETbBR5WPnFZcHiqvWJPyIKw9AiZOU8\/uvyA4VSXuiD25ARSOTuyhzB9mvHuUrvMdb0jAH5EvNd0aVFDHaeqWhS\/SX98zMjjeD1f2COgKrq4Tp1h29ui9HZ2BLJw8X9l7INjxTPgT6PAKgps9l8qaQ5zgmvREoiqg8zWX5otRzksnD58ELqshKyRzE16lGfq7q\/lV+uqNV8bu1RzFmfAhBjq8GrryDfrRN4cq94SuELPkvpguh5S5suFNWHRsbtpehQ8q0MnDZ3uFIeUy0B6CWwbQqzhLE8EVHscmc20qB6gkwZWXIKkaV16yZOIJXcUmd1J1W6epKLvLPgI6YqjiRPh5eM2vTh4QXCd0yW8iaQU3vzXqXDJVA18p9ShKOU3toRLQMC666x2j4vHwcNQ+u1HYp82nKHYsRz1snPGiymtzPKNO\/dNwxPT90ClHXBj0cG3QcxzXZK4tai4tizkRnTmQ8HyMzy9GxyYoZHKmvBIexN8e+EmAh4pX3PPk7fP69LxQew6XecXxhODq8IyHlHW2qCykEXNYqbRClUYUdavP3OVRE5MrbT1M8hBhIbhpUMsojsrwOsMjpkMwclVU6XNYVRh4k2ZEDotkTXll6kHMu0y9uej22BLGccZkeJE9odf92zsQXAseuI5qyiuRvxj2nTnss+MZo8qm09bI7V3U2ludNHx2ECC4rllokV5aAossdchlI8gQXAvjljGqKa9YcF9592Spp5DEKw67\/+jJkrexScMpC7uKu+yYxZvF\/PGpB8qIJLKWEVCNnIjRkSgvZjEqctugAV0+4sVoh5xXlPe7HI6fV8I5yxDcHL1DqoYxR0VNLEpR6lEpeKvKaWoPk3m2uI3pVWXi51GHGMTF43ze2PYUHX3+H3WSRhgPCJjySi5Smh3Nkqokr701rbpNPUzzguCaIpZheJ+Gz7Aapc3ys8wDaesjYMorkzN04w6t1ymdibBz76XHjlV0+LUXdJJGGA8ImPIqXCQWyajTpdIU3WZe2LYeJuWF4JqglXFYmyEaUTQeWolbiG46RMPhw84s27dvp3794s+X5OcNDQ108ODBjNFC8roIpOGVOM1HxSseGuQNL27u27Hzmc7Qn+AVb3zBmxTwJhsqXkFwda3tL5xKqFTTBCqRVMUXz8Nrb2UEdNLAHK4\/zuQqJ1sC65ALc7i5MrmTwojTp9ra2oL0pk6dGuzTLa40vJJ7rkm84kaTPZR5Li7uUsVnT+Vhvc+WxDoqHc6n13ON+JBzwhg3iaThlZwzz8vzFrRxO4fptmdJ0xo6aUBw3fCh4lJREbhSKlSUelQK3ps3b6b29vZAZMUa7wkTJgQH0POVxh46p7uYHIEWh2XUYQnhsCKfU\/80CYKbI1Km4ZVcfJMphbhqywcVpIXGth4m+WJI2QStjMP6NHyWVSlKPbLEKMu0V65cSXV1ddaCy3tlJ82xsRC+8s5HVvNwOsIuHGLebP4WBDdL4him7eI9t\/EB4OLaOEvpjgAZwpIYHILrEk3LtFQE1hkesRkyVqXP1VOFwTpcSxJYRpe34xswoGOYNw2v5N5HnM1drKcUnsq\/+0Fd7JAy5\/N\/vkC0eR4E15IeTqOn8Q0QbYgoyK5jNcQfXfLRjibLgoQPwbgRHQda5N3nBILrlIJ2ibkgMJym7GxQybHFXG59fX3kOctJznRh3ghPZRZCvqJ4xfOqU7\/Si+6\/5bOjIE2cpjhd0eD+329eGOuMx\/ksvb0brW76W\/Rwc0TQNB9ycvEFH9m+s+++usuhFKqPe54D\/tNZStxZSpWG7w4CBLeCCJyjoiYWRfUiVko98lrOAwcO0OTJk+no0aM0evRoWrFiBbHYNjU10bx580j0bG2HzJL2VBaHHMQdOm+CXVI+8hIl8MoE1ezDurQHc8CUSyovZ10EXNZDlScEV4WQx+c+DZ9ltYpSjywxcpk2DyMvXryYFi5cSLW1tV2STmuPJIcmk7W6qrombYIgO2alrYcqfzxPh4BLe8T1cuNKxtMMhz86Z7xvclR6LuuhQhKCq0LI43OV4XWGRzCH69FgOcmKnaTWrl3bqTTLly\/XdpqK45VwaIra81YWQltecj6nT5+mDVOGdkFUPrTc5\/KNnJg218Vw2V6ZnsMstnEUc7dxQOlw0yevILg5orRLAkdVS4d8SYLNaeqk4ZPAOTJfbouSlleiEYxyaJK9Q3U4kcSrpC0exYb0fHIReJUviqXllahFmDdhj2NbZ708tlcQ3BxxWOU0pTomTed52Isvi99oGHNEKslLOc5pKok37AUadUi87DAlHKq41lFOU6LhS3oe3kCDyyrvROXbuSVfFsxnaVSCa1pqMU0R5UAl0hIHa2z93k3Uv7ZjdzPby3U9kspjLLjygcTysJWciTzEJcKEd8PZuHFj4Enp6r7sSDJo0CBqaWkJ5rNM78v1Ew4pNTU15Op+FDYCO5+GtyVpUvyi1CNLjHymbWOPKIcm4TDFjd6I6z7zULapU1Q+4Xlim3rYlA1xoxHIqoMgthUNf8jxSMjTr75Pi77am+KWAaXpQPjsIBgJrrzGj02wbNkyam5u7uSowcLEosKCx0I1Z84c4mUKhw4dop07d3bxqNy9e7f1\/RkzZtCTTz5Jw4cPD+atuAytra20YMECWrJkifb9WbNm0bRp04Ide\/hjQOzgM2XKFGpsbLS+P3LkyEhsOC++VA2K7dCdbfw8DtGgMVQjYMOrqPnV8N61LnjFZ+qO\/FLvTtv8heeJfTaMalQRwoZXSW2JmGLgJWfLJ9xE\/IE3rfXtLk5SKt7lsb0yEtw4MRWCwRWUt5nj32LXG956Tux+I68XfPnll63vjxo1KhBYsSSCe7X8MfDggw\/S0qVLte+z2D722GOljwgh3N\/5znfokUcesb4\/ZMiQYCmH2OfWdEcgFcGyfp5HAqPZUyNg2zCGl2zozrWJkql4yeEa1r1OPXr06CS48h65GFJW29l3CFteJfGCRzce+Od9wbQCXzx8zFMb8qXDK1UY37wyFlwWNl73xxf3XkWvUgARJcp9+vQJRFXu4XKPkfd75cv2\/tixY+mNN97o1JOdOXMmPfHEE116vkn3H330UXrooYc69WQ3bdpEq1evph\/96EfW97\/73e8GZZJ7\/4yNEGAVgX2\/UGnzK0o90tY\/b\/Fs7RFe75i0bjZt3cNeqlHD1rb1SFu2ao1neyhGpeDmk1fOBVf0anmZAovJHXfcQZ\/\/\/Odp+vTpgUDziSY8x9q\/f\/9AIMeMGePk\/uDBgzttBnDq1KmgZ8uXvEmA6v7x48dp4sSJQTz++8EHHwQ92\/379zu5v27dumAJh4xNWHDTOLdweeE0VSmvuN9y2s61ffOZD+m2a3uRWB7E62a5tyHPl6WZOwtzlo\/qe3DUAKofegWteX4f\/ey1k\/TmotsDsHz3RPxaKJ+5ZXEoRh5rmmvBjZqflYeUw4CGh035edwWdK7uiyHl8Pyy6X0xpMw9ep6Plnvxoqdvcx9Dynl8\/YpXJlWDohp2Y\/Fb9OtjgcjyXBpf8qEGqviq50JQH951jrinyzsOhY9cg+CWn5dor+xtYNTD1XWaEmJ05MgRYkekVatWETtHiSPE5GHnl156yfo+z7tyHmIPWSYGX6JXrXufnaPE9nh8yLoYMr\/rrruc3P\/iF78YzDWzgMvYyJvM25s0HyngAPp82IFLwYJre538Rgt97swx+t9LelOPHavowmN\/sE0yMj7nw+fenvraSuqxY3WQp3yBV5nArkw07lAMZcQKCeCLV0aCy9jJy2PE0p7w1nLy0peo5T88nLp+\/fpgz1d5nsDmftT+stz7NL0v108+yNvV\/ShsKoSTKGYVIyA8RZOO63MBD+fjan2li\/IgjfgRSWBjjoCx4JpngRhAAAgAASCQdwRMD8XIe33yWD4Ibh6tgjIBASAABMqMgOpQjDIXryKzh+BWpNlQaCAABIBAtgioDsXINvdipg7BzYFddbbLLHcx5ZdPnmsvd7mQfzwC4BXYkQUC4FV6VCG46bFzElPH89tJRhaJxC3XskgSUTNGALzKGOAqTR68sjM8BNcOP+vYOttlWmdimUDUkgDLJBE9YwTAq4wBrtLkwSs7w0Nw7fCzji1vrsGJRW2XaZ2JZQLyEBInJS+Xskwa0TNCALzKCNgqTxa8siMABNcOP+vYlUBguZJieDm8h7Y1EEjAKQLglVM4kdinCIBXdlSA4NrhZx27EoZowpUUO3mJPaCtQUACzhEAr5xDigQ\/3fjIdHvfcgOXp\/YKgltmNlSCEwI33nyMIgssl1c+G7jM8CH7GATAK1AjCwTAKztUIbh2+DmJHbVdppOEHSYiLwvCHK5DYDNMCrzKENwqThq8Sm98CG567BATCAABIAAEgIA2AhBcbagQEAgAASAABIBAegQguOmxQ0wgAASAABAAAtoIQHC1oUJAIAAEgAAQAALpEYDgpscOMYEAEAACQAAIaCMAwdWGCgGBABAAAkAACKRHAIKbHjvEBAJAAAgAASCgjQAEVxsqBAQCQAAIAAEgkB4BCG567BATCAABIAAEgIA2AhBcbagqP+DOnTvp0ksvpUGDBtHp06fp6aefpsOHD9P9999PV199deVXEDUAAkCgMAgUsb2C4BaGnvEVOX\/+PG3dupW2bdtGU6ZMoWHDhtEzzzxDQ4YMocsvv5y2b99Oo0aNom7dulUBGqgiEAACeUagyO0VBDfPzIspm3xEVk1NTRBq8+bNNHfu3OD\/jRs3BqIqLibwmTNn6K233gpucQ93w4YNdNFFF9GhQ4eCHm7fvn0rEAkUGQgAgbwjgPbqMwtBcPPO1lD5hLCOHj2aVqxYQSy4Bw4coGXLllFzczPt37+fWltbacmSJcH\/H3\/8MV1\/\/fXUu3dvYuILweXDCO677z7q2bNnEJ5F9+KLL64wNFBcIAAE8owA2qvO1oHg5oyt8vFX\/fr1ozlz5pA47J0Fk3ukfPH8hhBcJrX4zQfENzU1BfF4uPiTTz6hHj16BL1ZIbiDBw+mxx9\/nCZNmkSXXHJJ0Duur68PxBsXEAACQEAXAbRXukh1hIPgmuHlJTT3WCdPnkxHjx7tMjzMBZAFlkWSf7e3tyvPqxWCy8PNe\/bsoeeeey5wouIe8JgxY7zUDZkAASBQLATQXunbE4Krj5XXkDzky4IrerFy5mkF12sFkBkQAAJVgwDaKz1TQ3D1cPIaigV106ZNgUfxvn37gp6rSnDDQ8rz5s2jAQMGeC03MgMCQKD6EEB7pW9zCK4+Vl5C8pxIW1sbNTQ0BPkxmfkaP358Kf9wDzfKaSqqZ+ylAsgECACBqkEA7ZWZqSG4ZnjlInRYcIUw87KgPn360Pr169G7zYWlUAggAATQXn3GAQgu3gcgAASAABDoggCveODVDjzixtfUqVO7TG8BNjMEILhmeCE0EAACQKAqEIha\/TBhwoRO01tVAYTDSkJwHYKJpOwQEF\/UvCZY3ilLpMrLmiZOnBj8XL58eenFj9tly8V9+Sufd+hqaWmh2tpaMr3Pc12NjY20d+9ekjctcXWfMVHhZ2cdxK52BNgTua6uDoJrQQQIrgV4plFNh2gOnzhHra9\/SLPvrjPNyig85zOt9W3a+v2bjOK5DCwLT3hrSs5HXmDPv8XOWsePH++yyxY7jB05csTJfd6DWqxxlhsc+etf574Iw+udeZhOfFS4uq\/Cz6WtkFb1ISC\/f1j9kN7+ENz02BnHVA3RXHPNNZ3SPPW1lcHv\/72kN\/XYsYouPPYH4zx1Ipz8Rgt97syxIJ+ev5od\/G9zHTx40Cg6f4hjnVqJAAAgAElEQVSsWbOGxo0bRzNnzgzmicI9XO7dsjhxD5M3+xCixTtvRS2J2r17t\/X9GTNm0JNPPtlppy\/eBnPBggXB1pnyDmBJ92fNmkXTpk0r1UvwgJd9ca9X1Dft\/enTpyfiF+aVkXFyFNiUVzkqekUXJW7kBLwyNysE1xwzZzHCQzRMYNGocM925Yvv0WMTbqBVL74X5JlFD1Tk8\/v5t9KXH3qVbru2F\/2k\/garOsr1MElI9NLiBJdFjXuvfIktL\/n\/qF22WIht7\/\/whz+kLVu2lHqjQvQfe+wxWrVqlfb91atX09KlS0msjRZemyzE8+fPt74vloDF4ZfWHia28xG2KPXwgZWrPJJ6tkWxh896QHBdMdMwnSgiy4Yf85M9dFVt95L41c74LbU19KPbBsVvZvHee+8lnmsb9Zzzue26y4Jha35+84\/b6cQjfxVbG5087rzzztKHgwksENyOPbFNhVhXcMO249\/i4vOQbZ7L6fD\/Im2RbtJzfsZHRPLe4eJcZjm+eM5r09HLNXmj7MLy+7h48WJauHBh4LcQvlRCpdNWJJ3DbRtf8FCVR9r2Kg26ENw0qFnGSRqi4QZFCN\/W791EI67rVfq96Ku9afq9f+msYeS5W+7V\/mzsFTRuxA1d8gk3nFk3jCrBxZBy9BC02IksbQ\/XtmGzjZ\/HhtHyFS9EdH7f1q5d26kusrNiHgS3259dSWMe30PclolL7jDocBOCWwi6RldCd4iGhZCHecX1\/da3g39th3vlUq18sZ1aX\/+PTvmEe9ZpTKF6EePSTBJcOE11eG7HOVkxpmkFN42NyxEnLa\/KUdZqyDMP9uCRvxHX9ipNt3GbxlNx8j2VLXzWAz1clTUcPtcdomHR40ues+W5VhbdpOFe06JG5SMIa5NPWgKHBSOMl7wsSPZkFst\/wrtsubgve5bLy3lM78texPIGAq7uQ3BN2Y\/wtgikfc9t8xXxWWxn33115CqOpGfh\/H3WA4Lryvoa6egO0XDvtn7olV2IxCTiXm\/\/2u6RuekMn8jzGZyeGLbmBDn+rmM1icKuk4fPIRoN2Ks+iKpB0bGpah7M5rngnioN8CpfVBa8spn7V\/kOiBqH5\/b\/\/vmTdPBPpwO\/Fr6i5v65fVt6ezf6+pevyo1vAAQ3RxwWBA4LoSgiCzF\/0dUPvcJJqePy4fvl6OE6qRQS6YKAqmGUGytZ9MKOS3ENq+wYFef0FJVWXEMpBDj8HIKbL3KrPuSyLC23hdxZiOt8cN5RI3hRZfJZDwhulqwwTJsN\/9SvdgdOAOGJf07q4V0djgEzh3Xv5I2cpmF8dv8FwfztLyd1iLfcUI7ecIRGfql3MF+cpqFEw2ho+IyD+2xQsqxKUeqRJUY+0y6XPUyGizkst2NJnRSf9YDg+mSoIi82\/NQn\/rWLI5OI5mJ+VaSV9PXHc8XvnziXet2vTwLnyHy5LUpR7FGUeuSWKIYFK4c9uG165d2TnRw9k4rN4dn\/JWnEzmc9ILiGJMsyOBv+xqZ\/CbKI2uRizfP7aNGvj8WSx2QuTl5\/K+ok4ic5aOnkgR5uliwxT1vVoOjYVDW\/avOca6RTBvDK3PZZxvDNK7GMUfRYVZwRvOK9BeJ6uZyGT15BcLNkpGHaTGDeZjHO8061MYWKgPLzqPlb8XzHOye7DGuHRTmuar4JbAhxVQb33TCGQVbxEoJbmbT0zSt5VzwdzogwYvpMXmYpt2cQ3Mrkn3WpheDKnsPhRF05TsU5TIn8dOY+4iqsehGtgUICRgjAacoILgTWRMD3ex41KqdZVEpq73zWAz1cXYt5CNf\/K\/fQ6RGz6Hc\/qIt1imLBHfTnF9KGKUNLJTJ1mvrdB+fooX\/rmAeJc4r65jMfBkuTvj3wkyAfE+9Tn1+MHsxS8Vn4bFCyBKso9cgSI59p+\/6Q47YvyslTp32SHU7D4X22VxBcnwxV5NXn3n+gy4d\/O9YhgMWRiRPn0KQauhPPo3aYCg\/RxO1spZOHTwLnyHy5LYpKqHRsajNHq0pfZ3gQUxX5o5dPXkXtgGfKq6h9DHzzCoKbIx5fMe4h+spfjUr0Dnax45TO+rQ4UdaBS\/Ui6qSBMO4QKIo9ilIPd5Ytb0o+7aGaAtNBIu40NJ\/1gODqWMpTmMvv+yk1Tfp64oHzLgQ3bicruZo2+fgksCfTVHQ2RbFHUepR0WSSCu\/LHjYf\/2Gsozb18VUPLgsEN0fsF2RI2tGHixt2aDKZw31l7wHijS3EurSkjS3EkYB9el6IOdwc8cS0KL7n2sQQsc7cGnwDTK2Zn\/C+hMrGWSpKcMNLhHzVA4KbH+4Gx0txz1PnLNo4T2XVnAY\/\/+B\/LtNe8hPlqayTB+Zwc0QsIlI1KDo2xRxuvmyah9L44BW3V9M2vR3p16Lirfjwk7kb9k3BHG4emJRhGeQTb+SzJZPWvoriCILFzUWoCCgEV5fAUUPPOnlAcDMkUIqkfTSMNoIc1TCGq+m7YUwBc9VF8cErPqSAr6iNgFRtURSvxFSZOATGN68wpOzxNYk707W2tpZ4nmL1L1+nYz\/9a2WJbLZeNJkPSXsGr+pFVFYQAZwiUBR7FKUeTo1bxsSynqrgw+X5o\/9nY6+gcSNuCGqaNAWm+5yn1B4cNSDYXxmCW0YCZZ019275iL6WlhaqqamhOXPmUH19x8HiLG4\/f+Fl+uNTDyiLYePQpOOhLAqQVtjRMCpN6DVAUexRlHp4NX6GmWVtDx71ixuNs6lWuCORdT3ksqKHa2M5w7gsuK2trbRixYogJgvu8OHDafz48cFv1RcjhxFHofH+oPKwiCiK6nzJ+7d8SEP6dqflE27q8sUovh5FPvyX53HljTi2b99O\/fr1y835koYmqMrgJryKOp5P5l3cc5l\/4Z5GFK\/CPRXwqvKoqRIq1ZCv6vldza9STffuscskVfEFD8PTHXKHBT3cyuOddol1BVcnQdVh9HFpmG7ZmCYf1YuoUz+EMUfg7NmznUZNRApFsUdR6mFu2fLFiPM5kTsIWZXOxdrbuLLJjqc+eYUeblZsiUg3aUjZlMCmwsnpC0\/opL2aw8VOk49PAns0X66zYv+AxsZG2rt3L23cuDGYpoDg5tpkuS9cks+JaXuVprIsilEHDqRJKxxHni7z2V5BcF1YTzMNlwTW2bwiXCzhCS2GonWKneawBJ8E1qlD0cNwz3bNmjU0btw4mjlzJs2ePRuCW3Sje6ifyw6CaXGjtnI0TUMVXoze3THkL+jgwYOq4E6eQ3CdwKifiDxEY9oTkecsojyIVXMaczftoRfeOZe4V3N4viO8BEmVh+85EX3kix9S9HJNBVfHpjbLflTps2VUYcAr\/\/y1nQLTsWkcr1gMZe\/kqNqr0lfxitu2e67rTpvnfQuC659e5c\/RxLmFDzF45d2TnU784RokOU3x4fVHT50PiBy1w0+Uc0vDutcDYMTpRDpu+bcsfJGG3XJL4p7Q5Ue7eCWIE1xuvIpyJW0MU5Q65qUeuoIbFr5wO2L6fM3z++hnr52ktoZ+JSh02yud9km0k9xpYQeqXs81QnDzQjqf5TAZik2zNMhkSZCoN6\/bXfnie4k7YIUx0tkT2ieu1ZJX2h5upeBj8n5USp3yXM5yDSn7GE4WuPPHaM9fzab2N1\/zYgoMKXuBWS8TVYMifylGzceqhlhuXLSd\/ubWq2IPR4iKHxZ2VR78nJcshfcr1UMAoWwQSCu4OjbFkLKNZSozrq3PSVpeCe\/kvt0+6nQueBhFVfocXhWG28Rjb75MH\/5ivhcjQXC9wKyXiYngcopRhxgkNYwqN\/socoaFXUVgIbgmntB66CCUCgEIrgohPDdFwJXPSVS+UW2JvBOeTluT1N7pCK5or3xNVUBwTRmYYXiV4IazNvUgVgluXNVMlwalWbubIaxVn7Qpr\/IKWFHqkVd8TcuVhT18DifLw8q+OggQXFOWZRhe5TQlOw7w\/8KDeOaw7kGpws9FUYUjAQ\/18q5R4bCq37z3qBiK1nFKGLryDa09oTOEEklLCJjyKswbFa+Eo56KR7bPcShGvmjtWnDT7BPgApHeD\/ycxo+8MZgGy\/qC4GaNsEH6KgKHh1iijpqKG2LhoZqnX32f3lx0e2yJ4oZw5HxUwzzs1fzC9t9p7QltAA2CWiBgyqtwViqb2z7n\/HTSgOBakCCDqK55FT5YRYcTLoaUffIKgpsBEdMmaUrgsAdxEkF5qObsuXP0UtOtqQT3\/RPngmU+qpeA9z\/d9+pvvDkhpMW6muKZ8gqCW03sSF9X1ciJGNHQXRYkDpr\/9sBPuhQq7bKgvO3RDcFNzzfnMVUNYzhDk6VBgsyz7+4YUja5TPLhYe4\/7nyWjj7\/jyZZIGyGCJjyKsOiWCVdlHpYgZCTyCanm+kWOa2PiW76ceF88gqCa2sth\/FNDS\/mPHQ87Ewdn+RqmQgu59Njxyo6\/NoLDpFBUjYImPLKJq8s4xalHlli5CttsXpBp+3RKVPao0B10laF8ckrCK7KGh6fpxmiEcfniWLG7TTFDlO8w9TNfbs6WHFc2fFF\/Oa\/ssMVv1yqIRoIrkfCaGaVhley7WV+4Hg+TdCrIJjK2Ug1\/SQ\/j+oQmMSPg1snDczhFpSs4vi0tra2oIZTp04NNpoXl+pLK4o88iEGceQSX6PyubZRECeRUyxBGtb7bOxidO4JL912gM481eBtq7SCUsVptdLwSi6ATqNlszGGEHdVGj4bRqcGKGhi\/e74Gzoz+O9id6Ez4U3UyUAm8SG4BSWZTbU2b95M7e3tgciKTQomTJjQ5QB6kzyiDjEIxw97\/5mkL8KGDzGISkPkc+qfJkFw04CcURyV4GaUrfNki1IP58CUKUG2x8lvtFjvKpdmy1mXVfbJKwwpu7ScYVorV66kuro6K8HV2evYBaF15lhEPm82+zt9wxDyqgzus0HJEuCi1CNLjHymzfa4Z+m2IEubNazlcpbSHVl0iSkE1yWaBmnJ+5QOGDAgiKlqUKKGWGSHprghGNndXjVsF\/dc5JM0LF2O464MIK\/aoGl4JYNlO7Snis95qcLwcwwp54vCglfhneWE3weXNun0Mn7+7P4LuuwPIMcXNcayoHzZvqJKI+Zy6+vrIw8K1123JirNhOejrPr0vDCS4Px80Vd70+gbepRwSkNgTmfd3d1o8MCrIo\/34+dLb+9Gq5v+FkPKOWIknKZyZIwCFUXwKm47RtVH1Ct7DxDvYhfn6ayKr3qexw859HAzfAEOHDhAkydPpqNHj9Lo0aNpxYoVxGLb1NRE8+bNI9GztR3aSFryY7J0SAVFUj7yMgFVj0qVD567RaAo9ihKPdxat3ypCXuI0a\/fz7+V+td2rILQuTgeH\/3J8cp5+eQVBNejpXkYefHixbRw4UKqra3tknNawyc5NLEQTtv0thNSyx7R4cLLjllp6+HRFFWVVVHsUZR6FIV8sj10nCrlenNHYMzje+i2a3tZzf+6wNInryC4LiymmQY7Sa1du7ZT6OXLl2s7TSXtdcxbLz56b68uS3ZcHnfFjlOnT5+mDVOGdqmxcJjiMmCuTZMQnoKpGhTV0FzWz\/M49OfJNBWdjcyrqF5uEm\/EEsLnpw+N7RXb8i6PvILg5ojyaRvGJIcm2UPZlsBJHtGyYxYEN0ekSumMJ9fAljeq+HlsGPNlwXyWJtxehXu5qnX9g\/78wsiPd1FbFW9Uz\/PIKwhujriscm5JOiaNd5IScyiCiGL+dupXetH9t3T0foUHYJTTlCAo\/417Lhy0bhvU4VnN6R09db7k\/ABv0hwR6tOi2PAqzIWonaZc8EqHl\/iQyxe3woJrsgVs1EYX5aqdqqPjslwQXJdoWqZlY\/iow+iFI5PLw5WjHKfC+6qmqceuXbto4sSJAYLyMLsMqTwkL8KEd+\/auHFj4Pnt6r7s+DZo0CBqaWkJ5t9N78v1Ew50NTU15Op+FDYCuzT2sKRyJtGLUo9MwClDolEfcjzSdfBPp4NVE3EdBG5DeJvZcSM6zp+NW5VRxA85CG4ZiBqXpapBSRpCiZpfzeJ8yRsXbaeRX+rdydEhPE9s2hOR1yQzNsuWLaPm5uZOjmUsTCwqLHgsVHPmzCFeVnXo0CHauXNnFw\/w3bt3W9+fMWMGPfnkkzR8+PBgnp3L0NraSgsWLKAlS5Zo3581axZNmzYt2GGMPwbEjmNTpkyhxsZG6\/sjR46MxIbz4suGV7rDcmnXd5sMH5ryKkevdiGLEscr\/vjvf1n3SJ8SedhZNSRs+1yXuz55BcHN0atg2zDyl6O8pi28w5QLAvMB8x\/994XB2bjiktfhpRlSjhNTIRicj7wtJv8Wu3TxVplity55ffPLL79sfX\/UqFGBwIolXNyr5Y+BBx98kJYuXap9n8X2scceK31ECOH+zne+Q4888oj1\/SFDhgRLz8S+3KY7mNnywjZ+HhvGHDULuS1KUnvFbREflCLO344abrbljSp+HnkFwc0RnW3n2ngel4eP+3b7KBjOEVum8W++XAzRrHl+Hy369bHSfLFYvC6GrdMKLgsbr1Pmi3uvolcpzBMlyn369AlEVe7hco+R96fmy\/b+2LFj6Y033ujUk505cyY98cQTXXq+SfcfffRReuihhzr1ZDdt2kSrV6+mH\/3oR9b3v\/vd7wZlknv\/jI0QYNWHXI5egcSiFKUelYK3zWEr7D8yrfVt2vHuyaC6vD63fuiVlOY87qzx8skrCG7W1jRI39bwYS\/BrBwT5PniqMXrpvUQPb4kwRW9Wl5WxWJyxx130Oc\/\/3maPn16INB8AhPPsfbv3z8QyDFjxji5P3jw4E6bl5w6dSro2fIlb2qiun\/8+PHSHDXPVX\/wwQdBz3b\/\/v1O7q9bty5YciZjExbcpLmyqN5AlCMTnKYMXugKD5rFYSt5hMS0vbKpAwTXBj3HcVWGVw2hyL1P\/rrkSx76VcVXPReN8sO7zhELLQ9fRy0FMJ0T0RlSDkMdHjbl53FbZrq6L4aUw\/PLpvflDwyej5Z78aKnb3MfQ8qOX0wkFyAAXtkTAYJrj6GzFFhwbS8+LutzZ44FyVyy+5\/owmN\/sE0yMj7n0+u5Rjr1tZXUY8fqUp4i8MGDB7Xz1XWaEmJ05MgRYkekVatWETtHiSMPZeF+6aWXrO\/zvCvnIfa85gaHL9Gr1r3PzlFiO89+\/fqVhszvuusuJ\/e\/+MUvBnPNPEIgYyMfiqFtjJwHNOFVzqtSUcWLO2yloiqRUFhfvILgFoUxn9ZDzJ3IPdssqsj5mOybqiqDvDxGLO0Jb4UpL32JWv7Dw6nr168P9qiW559s7kfth829T9P7cv2mTp1aml91dT8KGxXmeA4EdBCIGyHSiYswnRGA4IIRQAAIAAEgEPkRmXTYCiAzRwCCa44ZYgABIAAECo+A6rCVwgOQQQUhuBmAiiSBABAAApWOgOqwlUqvXznKD8EtB+qhPHW2NSx3MeWXT54TLXe5kH88AuAV2JEFAuBVelQhuOmxcxJTx0PXSUYWicBpwgK8MkUFr8oEfMGzBa\/sDAzBtcPPOnaaNajWmRomELUkwDAJBPeMAHjlGfAqyQ68sjM0BNcOP+vYurssWWdkkYA8hMTJyMtaLJJF1AwRAK8yBLeKkwav7IwPwbXDzzp2JRBYrqQYXg7vdWwNBBJwigB45RROJPYpAuCVHRUguHb4WceuhCGacCXFjktir15rEJCAcwTAK+eQIkGi4IjKqGMy5ZO98gZUntorCG6Z2VEJTgj8kvFxdyywXF75DNcyw4fsYxAAr0CNLBAAr+xQheDa4eckdtS2hk4SdpiIvCwIc7gOgc0wKfAqQ3CrOGnwKr3xIbjpsUNMIAAEgAAQAALaCEBwtaFCQCAABIAAEAAC6RGA4KbHDjGBABAAAkAACGgjAMHVhgoBgQAQAAJAAAikRwCCmx47xAQCQAAIAAEgoI0ABFcbKgQEAkAACAABIJAeAQhueuwQEwgAASAABICANgIQXG2oEBAIAAEgAASAQHoEILjpsUNMIAAEgAAQAALaCEBwtaGq\/IA7d+6kSy+9lAYNGkSnT5+mp59+mg4fPkz3338\/XX311ZVfQdSgLAiAV2WBHZlWIAIQ3Ao0mmmRz58\/T1u3bqVt27bRlClTiDcaf+aZZ2jIkCF0+eWX0\/bt22nUqFHUrVs306QRvooRAK+q2PioeioEILipYCtvJPmIrJqamqAwmzdvprlz5wb\/b9y4MRBVcXHDeObMGXrrrbeCW9zD3bBhA1100UV06NChoIfbt2\/f8lYKuZcdAfCq7CZAAQqOAAS3wgwshHX06NG0YsUKYsE9cOAALVu2jJqbm2n\/\/v3U2tpKS5YsCf7\/+OOP6frrr6fevXsHR2sJweXDCO677z7q2bNnEJ5F9+KLL64wNFBcVwiAV66QRDpAIB4BCG7O2CEff9WvXz+aM2cOicPeWTC5R8oXz5sJweXGUvzmA+KbmpqCeDxc\/Mknn1CPHj2C3qwQ3MGDB9Pjjz9OkyZNoksuuSToHdfX1wfijauYCIBXxbQralVZCEBwc2gv7rFOnjyZjh492mV4mIsrCyyLJP9ub29XnlcrBJeHm\/fs2UPPPfdc4ETFPeAxY8bkEAkUySUC4JVLNJEWEDBHAIJrjpmXGDzky4IrerFypmkF10vBkUmuEQCvcm0eFK7gCEBwc2hgFtRNmzYFHsX79u0Leq4qwQ0PKc+bN48GDBiQw9qhSOVCALwqF\/LIFwh0IADBzRkTeK6tra2NGhoagpJxI8nX+PHjSyUN93CjnKaiesY5qyqK4xEB8Moj2MgKCMQgAMGtQGqEBVcIMy8L6tOnD61fvx692wq0a7mLDF6V2wLIv+gIQHA9Wpg9iNl7mHuwfE2dOrXLcLHH4iCrgiAAXhXEkJrVYOfHiRMnBqHjPrB5RKOxsZH27t1L8hJC0\/tyXsuXL+800qZZXASTEIDgeqRDlDfxhAkTQGKPNihiVuBVEa0aXSd5eRf7aIh5+ZaWFqqtrS1FYue4urq6YPUBf+Tzsj9enWByf+DAgcESQ\/YH4Uus9ZfzqR7k3dQUgusGx1SpCPLL87OpEkIkICAhAF5VDx1k\/w0hhKIXy86WLLLig4ydMLnXq3t\/5MiRgUCzmPPyQ1m4qwdhtzWF4LrFUzu18JeqdkQEBAIJCIBX1UUPeXRD1DyqF8yrGGbNmkXz588Peqyid5x0f+zYsbRly5ZgaSJf8iY81YWyu9pCcN1hqZ2SmHMTwzwi4jXXXKOdRp4DHjx4MM\/FK2zZwKvCmjayYlF7X3NAX4KL9sqcbxBcc8ysYiT1QJjARRCrotTDytCeI4NXngEvc3ZRPVu5h2sydBw31KwaUla95++9917isZ9ZP2c8dPK48847vbW7EFyPLw43iosXL6aFCxd2cnCQe7hJgqtDnqRzbW3j55HAHs2X26zAq9yaJpOCRa3ND2dk4hwV50ylcppSCW4mlc8gUZ\/1gOBmYMC4JPklWLt2bafHsqu9yvC2gmkbH4LrkSwGWVUCr7r92ZU05vE9dPjEuVLNTjzyV6X\/dbjpsydiAL\/XoPJ+2CJjPm6THZvWrVtH3CtlAZWX\/8jLD03vy8uCwsd+qtorr8BYZOazHhBcC0O5jlo747fUv7Y7\/X7+ra6T9pbejndO0tglm+iPTz3gLU9klIyAzwYlriTM7RHX9qKt378pCLLyxXZa+eJ7ne6p7JiHeqjKWE3PVfbQ+YiqthE5CG6O3hAm8D1Lt1Hr6x+S\/PWfoyImFkU0ohcdfoU+\/MX8Sil24cupahizBoDFdvbdV9Psu+u6ZJX0LBy43PXIGqdKS5\/t8Zvf\/KY0T8sCy5cQ0e3btxMfMSp+2zwPx+V8OP3bb789yDPquXwv6bnPkRMIbo5YLhqUhnWvU9vbp2nr926iEdf1KpFJkDn85SjIpPNcVDfqJQinE0XkpJfoyw+9Sj8eeT7YB7oIzl85ooZVUQSvbHjDfEmKH8erv3\/+JB3802lqa+jXqTGWG0AW3aW3d6Ovf\/mqxMYZvLKigfPIRfkA8lkPCK5zGqZPUDb8mJ\/soR3vnqyYnq7cU\/FJ4PRoV0\/MctqDP8L4w5GnSuIu5jpfYrg5Llw561E9bNGvaVHs4bMeEFx9fmUeMmx4boiuqu1OP6m\/IfO8bTIIN5g+CWxT7mqJWy57mAwXc1jmef3QK2LNUq56VAtPTOupsgfmcLsiCsE1ZVmG4aMIzA2RGFrOMOvUSfN8Mzu\/yL0Y1YuYOjNETIVAOezx\/da36ZV3T2o7AHJ4le9COeqRCvAqiaSyBwQXgpvrVyGKwHLvMY8E5g8C7pWIXjiX0acTQq4NmpPC+W4YeekPDyWLHquKtwwTh7n5x+2xvVzwKidkkoqh4lX+ShxdIp\/1QA83R6yIM7zo5fbt9lHudm7hhlVexoSGMUeE+rQoqgZFJYimz8Woh+CFKr4Q3Gf3X0Ctr\/9HZK8YvMonr5K8lMOewTa\/47yM4zygBad0nqveD5fIQ3BdommZVpzhdZ1KLLM3jh7u3YoEfBLYuNBVGMG3PZivt113WeQyIBX8SVMovuuhKmu1P1fZQ\/WhlfXzsOhG2cv3hxwEN0dvTdLyDR5uk5cJyQvGo74c457bfGUKAvPfXcdqgrnbX07qcHKRvyQxpJwjUhGRalmQbDsXvOJRjyhehHkS9fvhXR07Uc0c1uHVDF7li0tyaSC45raB4JpjllmMJAJzr+Gyi87ThilDY\/P3+cUY50Ht+4sxM2MUKGGfDWMUL1S8DPdEuJfLw9HyUiLwKn+EVPEqfyWOLpHPekBwc8SKJMPzvBh7coYbonIUn7dvnLapoyxRl08Cl6P+lZanT3u48KrnHvJt1\/bqshzOZz0qzcblKG9R7OGzHhDccjA1Jk+V4bkhqh96Zaq5MZfVVM0pq+rhsixIS42AL3vw1p5xTk\/qUnYOwcId3t7UVz1My1qt4TECrU4AABGqSURBVNke2NrRzPoQXDO8Mg2talBEL7ec+yyLJR9Ja4NV9cgURCTeBQFf9rBxlgoXOmojDF\/1AIX0EFDZQzWVkPXz8FRFVK18T1VAcPW45SWUDoFVaxWzPn1j7qY99MI752KHk30T2IthKjwTHV7Z8uaD\/7ksdppB1bBGNYw8fcIX1nfnl3wqXuW35J1L5rMeENwcsUJleG64hBdn1HaPqobN9jlDdVfzq\/TVG+OHtSG4OSLUp0XR4ZWt4PIhBXxF7Yes4l2U4IZ9FsCryuNV\/kocXSLV++GyHhBcl2hapqW7fCPsxelzWRDn\/bOxV9C4ER37O0fljWVBlkRwHF2XV2Fh1OUVHy7P\/gUqXnC14jYiiMpr9IYj9OCoAcFOZhBcx6RwkJxKqFQfWlk\/x5CyAyMXOQkVgUXduXHj80WTNnrPAieVs5TIU7ceWZQRaXZFIGt7qLzW09okPKycdT3SlrNa46nskbWgqtKH4FYrMzXrrUtgbojeP3Guy\/CdioC2z8O926hqoSeiaWyPwXR5FVckFW94mqGme\/fY4\/VU8eMaRtlJELzySBjNrFReyq422RFnMcsjJGHOhPMyea56PzTh0AqGIWUtmNwF2rVrF02cODFIcPny5TR+\/PhS4iaGj9ocwF0po1MK75scl59JPbIuc7Wk74pXafBysfY2Ll95NAe86owS27y1tZVWrFhBNTU1nR6eOHGCGhsbae\/evTR69OhSGNP75eRVGi6mieOTVxDcNBZKGYfJ3tTURPPmzQtSWLZsGTU3N1NtbW3w28TwcZsDpCyaMprJ2bwm9VBmjABKBFzySplZRADdD7E0acujOeDVZwhu3ryZ5s6d20lMZXxXrlxJdXV1NGbMGJozZw7V19fTsGHDyOT+wIEDrdor1chG1s8ZD508fPqcQHDTtAIp4\/DXIhO+paUl+CKVXwRTwfW9JtekF4OGMSVBUkZzySvTIph8iJmmLcKL0Zw7hvwFHTx4MG0yhYnH9j506FBQn507d3bp4Ype7OzZswORZXFub2+nKVOmBL1e3fsjR460aq90xM7WOz4pPgS3MJRPVxF5CIhTYMEdPnx4aVhZJVRhAodFMCuCy8et6eTh84sxnSWKFcs1r8LoJNlcd17fpmHkHvQ913WnzfO+BcGVjMNCGie4YiRtwIABgeByuFmzZtH8+fODETad+2PHjqUtW7YEgp6mvaqUt0zV7rqsB3q4LtFUpKXbMMYtz+DkhQMB\/xW9C3GySvi5KI5wKJCLF7U8Qw4nP79\/y4c0sE\/H3rY6jhAQXI+k4pObpLm8pIZRl1dRvJF5J56veX4f\/ey1k9TW0K9UYRNeCb6Knkjcbx5W5o++Xs81QnBzJLj8sVWUy9fufRBcj4xxPfTnY1hZZyvHMIQ+vxg9mi+3WbnmlW5FfQwni7Jw497zV7Op\/c3XdItX+HBJPVyToeO4oWYMKbunEATXPaaxKdo6t0QN7cmbvOsM95rOmYQ3pNfJAz1cj6Qioix4JdcgzuZiSqNvt49KG1pE1VzFGdHDTeLmjYu207E3X6YPfzHfL7g5zi1OcLnIJs5Rcc5Utk5TOYauU9F8dhAguJ5ZIbvZb9y4MXBqEJfK8FENl7w5gKphS\/M8fEKRThoQXM+k+nRYWSw3c8ErleDKH2I6nLCZwxWCzPuI+xr6829B8xzDgssiy71SblPk5T9Tp04NHKX4Mr1v016Z16g8MVTtrstSQXBdommZVhrDZz2szIL72IQbaMR1vbRrl6Ye2okjoDECWdjD53CyqLCJp7wxSIhgjEAWvDIuhIMIPusBwXVgMFdJpDV8Vls96m7lGK5\/2nq4whHpdEbAtT3SzOu7sEnvB35O40fe2OVgehdpIw1zBFS8sh35sI0vRkZU02g+R+QguOY8yyxGWgKLzQEevbeX1VyaatlRHgmcmTEKlHBaXgkIwrxIM6\/vYkjZZ8NYIPNnVhVs7WgOLQTXHLPMYoiGMc3yDZ7fEqe1JMUXhddZvsE9519OuiKIIsJv376d+vXrF3vqCz9vaGjA8o3MWGKesA2vZL4IXomD5r898JMuhdHhlfhwA6\/MbZmnGKoPuTyVNaksPusBwc0RK2wM73pOzSY9m3rkyByFKAqPfvz8hZfpj0894Kw+5ZpLBa+cmdBJQip72A4J28bP44gcBNcJ9dwkoiKwKhd5iZAqrOq5TaNqWw9V2fBcHwE+Om\/M43uceffGnVSlX6L0IcGr9NhlEVNlD1vBtI0Pwc3C6gVK05bALJK8G1TcObm6BOZG9ZV3T9Lv59\/aBV2dNDDXli9SqpyNdGwqhoqjOGYSPw4ZnTTAq3zxStVe5au08aXxWQ\/0cHPECpXhVY1Sw7rX6aP\/vjD1uaScfrc\/u5J47jZuvaOqDPwcDWOOSEVE\/e74Gzoz+O+sbCoEN+pkIB1OwGkqX5xwURo4TZmjCME1xyyzGCrBVWXsYvhQPqhAlV\/cc9t6pM0X8aIRYHuc\/EZL4uiHDnZpl4nppK0TBrzSQclfmKLYw2c9ILj++KnMyYXhbZyduIC28TkNF\/VQgoUA2giwPe5Zui0Iz1MOaS+bef20ecrxwCsXKLpLQ2UP25EP2\/hcU500fI7IQXDd8c86JRfLN3YdqyGeg5WHhJl04Stu+cboDUdo0Vd707gRHQ2ziItlQdbmLVsCglfiXNn+td072ZZ\/RJ0GJPPm2f0X0NOvvk9vLrq9VA8TXolIUbzjZ1huVjZ6pM5YJbipE\/Yc0Wc9ILiejZuUncrwOl9r3KDxPNtt13Ycpydfqvh3Nb9KNd27x84B5\/GLMUfmy21RBK\/iRi9UvHhl7wHiDzGbeX3M4eaWHqkLpmqvUifsOaLPekBwPRvXRnBNihruzejEdTVk6JPAOvWq9jDCHmLfbfY+F71cHWxczOvr5KMKA16pEPL7XGUP1Ydc1s\/z2EGA4PrlaGJuKgKbFDWulxuXhkuHGJf1MKkzwkYjINvDlBe8bzKv440aMfGNN3jlG\/Hk\/FReyjrTBKpd68Tz8NQWl4zTv\/32jimOqOfyvaTnmMPNF6+8lUbVoJh8EYrezNbv3VQ66ScpPjfE8+\/oVZq7jau0Thl8EtibcSo4I5lXUb3cJJty+KXbDtDz04fG9op1OIEh5QomUEzRVe1VpdTYZz3Qw80RK1SGN23Ywr3WuPii1zNzWPfEww\/yOESTI\/PltihhXoV7uaoPsUF\/fiFtmDI0tn6mvIxKSCcNfMjli2Kq9ipfpY0vjc96QHBzxArXhhfHqM2++2qafXddZE1ZlA9\/dC5yV6m00LiuR9pyIF4HAmF7mJyhHLXRRblwBa\/KhXx0vkWxh896QHBzxGHVsiB5SYU8RBeenwj3FtgZ6nc\/qCst\/eAqc3w+Zo2Xeqz95hV026ABkXMhUWmL+KLHG\/6NnkiOSCUJrswL\/tA6+KfT1NbQ+eQnmVfMG3ECVdToRtIcWRregFf54o2qNCqh0hm1UJ1Va\/M8jyNyEFwVqzw+z4rAokcz9Su9aPmEm0h4nXLV5P2SVS9IHgns0TwVm1Ucr7j32v+y7hR1jrI87Kzihe1z8KoyqaVymnL50RXn9BS3rjvMqaT4qnbXpXUguC7RtExLZXibho23fVz32wPU9vbpkvNL+HACVfpoGC0NXKboSbziXuzNfbvTS00dB1VEDTereGH7HLwqEzEss1W1V5bJe4vusx4QXG9mJTp79izNmTOH2traglynTp1Ks2fPLpXAp+GzrHZR6pElRi7TtuEVz\/NPa32bdrx7MigSr8+tH3pl7Jy\/y3KbpgVedSAm23vQoEHU0tJCtbW1neA8ceIENTY20t69e2n06NG0YsUKqqmpIdP7u3btookTJwZpL1++nMaPH6\/dXtl+iNnGz+OHHATX9K23CL9582Zqb28PRFYQf8KECSUSF6VBKUo9LEztNSp45RXusmcm23vlypVUV1fXSQi5gOL+mDFjgo\/8+vp6GjZsmNH9gQMHUlNTE82bNy+o87Jly6i5ubkk7qr33FYwbeNDcMtO1XwVIPyygMD5sk+llga8qlTLqcsterfDhw8PRJZ7oK2traUeLKcgPub5w55FVgj0lClTgl6v7v2RI0cGAs09aO4dy8LN+ajaK3Vt8hHCZz3Qwy2TzfmlEF+PAwZ0eAiz4YtyHTx4sChVqah6gFcVZS7jwgrBFT1WFlwhimJYOcwBFtydO3fSrFmzaP78+UGPldsc1f2xY8fSli1bAjHniwVXCD3aK2PTBREguOlws4oVfmmsEkNkIPApAuBV8amQJ8EtPtruawjBdY9pKcUDBw7Q5MmT6ejRoyXHBX5hwj3bDIuApAuIAHhVQKNqVilPQ8qaRUYwCQEIrkc68FDP4sWLaeHChV28Cj0WA1kVDAHwqmAGVVQnL05T1YW6m9pCcN3gqJUKz7WsXbu2U9iwq71WQggEBCQEwKvqooO8LEhe8sM8YEcndpSSl\/\/Iyw9N78vLgjZu3BikjSs9AhDc9Ng5i5m01s1ZJpYJyY16nz59aP369YHjBa78IgBe5dc2lVwy8Cq99SC46bFzElP2KOQEw2vdnGRimQiccSwBLEN08KoMoFdBluCVnZEhuHb4WceW3fqj1rpZZ+AggailJg6SRRIZIgBeZQhuFScNXtkZH4Jrh591bHnhOicWXutmnYGDBOQhJE4uvCWlgyyQhGMEwCvHgCK5AAHwyo4IEFw7\/KxjVwKB5UqGlyVYA4AEMkEAvMoE1qpPFLyyowAE1w4\/69iVMEQTriQ7UPElH7xgDQQScIoAeOUUTiT2KQLglR0VILh2+FnHrgQnBH7JXn755U6HLoj9WK0BQAKZIABeZQJr1ScKXtlRAIJrh5+T2JWw1k1eFoQ5XCdmzzwR8CpziKsyA\/AqvdkhuOmxQ0wgAASAABAAAtoIQHC1oUJAIAAEgAAQAALpEYDgpscOMYEAEAACQAAIaCMAwdWGyk9A+fDw8EHirksg773qIm3ZucpFekjDHQLglTsskdJnCIBXZmyA4JrhlXnorEVWVCArcXQt4pkDXiUZgFdVYmjP1QSvzACH4JrhlWloPnZr7ty5QR58ilB7ezvV1dXR4MGDgz2Wb775ZnrkkUeIDw94+OGHicm+d+\/eTjs\/yWelDho0iFpaWrocBcibVyxdupQaGhqCAwjkOOGDCWTvZPlkI\/nUETkfTmvDhg304IMPEm9Viav8CIBX5bdBEUsAXplbFYJrjlmmMaKGaFhw+SD7H\/zgBzR+\/PhAaNva2oITe\/iaNWsWrVq1ir7whS9QY2NjsF6Wj9HicEePHqUVK1Z0Ej9ZFDk+bydZX18fxJF7vvxC7dy5M4h\/5MiRoAws9Cywchz5fM6wmGcKFhLXRgC80oYKAQ0QAK8MwCIiCK4ZXpmHjiOwEFXukcoCJy9EP378eCCyolfLwhp1+pAsqnEnAUVt4SjKJnrczc3NXXrPDBCGlTOniXEG4JUxZIiggQB4pQGSFASCa4ZX5qHjCCwLZ5LgTpw4sVMZo86uleNz4KjhYXFyEfek5Ys3veBDrmVhD4MCwc2cJsYZgFfGkCGCBgLglQZIEFwzkHyGtiVwa2trlyHkcPmTHKbE5uQLFiygJUuWlIaa5TTies4iDATXJ2P08gKv9HBCKDMEwCszvNDDNcMr89A2BA7P4XJPdtOmTV0cp+Q5XB46bmpqonnz5gUOVHFzuByO54cnTJhAd911V6e5YnmulwGSHbIyBwwZaCEAXmnBhECGCIBXZoBBcM3wyjy08PyL8lIWc6ZxQ8oqj2NR+LBjk7w3atizWfZSlvdQjvOGhpdy5hRJlQF4lQo2RFIgAF6ZUQSCa4ZXYUJjHW5hTJmrioBXuTJHYQpTFF5BcAtDSfOKuJ5rzeqlMK8ZYpQTAfCqnOgXN+8i8AqCW1x+omZAAAgAASCQIwQguDkyBooCBIAAEAACxUUAgltc26JmQAAIAAEgkCMEILg5MgaKAgSAABAAAsVFAIJbXNuiZkAACAABIJAjBCC4OTIGigIEgAAQAALFReD\/A2zGBNtvUZoCAAAAAElFTkSuQmCC","height":0,"width":0}}
%---
%[output:70d065db]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"ans(:,:,1) =\n\n   1.0e-15 *\n\n   0.0000 + 0.1110i\n   0.1110 + 0.0000i\n   0.0000 + 0.0000i\n\n\nans(:,:,2) =\n\n   1.0e-15 *\n\n   0.0000 + 0.0000i\n   0.0000 + 0.0000i\n   0.0000 - 0.1110i\n\n\nans(:,:,3) =\n\n   0.0000 + 0.0000i\n   0.8622 + 0.0000i\n   0.0000 + 0.0000i\n\n\nans(:,:,4) =\n\n   1.0e-15 *\n\n   0.0000 + 0.0000i\n   0.0000 + 0.0000i\n   0.0000 + 0.1110i\n\n\nans(:,:,5) =\n\n   1.0e-15 *\n\n   0.0000 - 0.1110i\n   0.1110 + 0.0000i\n   0.0000 + 0.0000i\n"}}
%---
%[output:54145e5e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: An error occurred while drawing the scene: SceneModel error in command compositeCommand: TypeError: Cannot read properties of undefined (reading 'deferred')\n    at https:\/\/127.0.0.1:31516\/toolbox\/matlab\/uitools\/figurelibjs\/release\/bundle.mwBundle.gbtfigure-lib.js?mre=https"}}
%---
%[output:0c5099ec]
%   data: {"dataType":"textualVariable","outputData":{"header":"logical","name":"skip_non_poly_LMI","value":"   1\n"}}
%---
%[output:635343b0]
%   data: {"dataType":"matrix","outputData":{"columns":2,"name":"ans","rows":1,"type":"double","value":[["33","33"]]}}
%---
%[output:2850b3b5]
%   data: {"dataType":"matrix","outputData":{"columns":2,"name":"ans","rows":1,"type":"double","value":[["66","66"]]}}
%---
%[output:163145a2]
%   data: {"dataType":"matrix","outputData":{"columns":2,"name":"ans","rows":1,"type":"double","value":[["11","11"]]}}
%---
%[output:4f20a75e]
%   data: {"dataType":"text","outputData":{"text":"+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n|   ID|                            Constraint|   Coefficient range|\n+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n|   #1|   Matrix inequality (complex) 186x186|              1 to 1|\n|   #2|   Matrix inequality (complex) 186x186|   0.20412 to 800000|\n|   #3|   Matrix inequality (complex) 186x186|   0.20412 to 800000|\n|   #4|   Matrix inequality (complex) 279x279|             1 to 10|\n|   #5|           Element-wise inequality 1x1|       1 to 10000000|\n+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n","truncated":false}}
%---
%[output:2007d292]
%   data: {"dataType":"text","outputData":{"text":"\nMOSEK Version 11.0.24 (Build date: 2025-6-25 11:19:55)\nCopyright (c) MOSEK ApS, Denmark WWW: mosek.com\nPlatform: Windows\/64-X86\n\nProblem\n  Name                   :                 \n  Objective sense        : minimize        \n  Type                   : CONIC (conic optimization problem)\n  Constraints            : 496             \n  Affine conic cons.     : 0               \n  Disjunctive cons.      : 0               \n  Cones                  : 0               \n  Scalar variables       : 1               \n  Matrix variables       : 4 (scalarized: 364095)\n  Integer variables      : 0               \n\nOptimizer started.\nPresolve started.\nLinear dependency checker started.\nLinear dependency checker terminated.\nEliminator started.\nFreed constraints in eliminator : 0\nEliminator terminated.\nEliminator - tries                  : 1                 time                   : 0.00            \nLin. dep.  - tries                  : 1                 time                   : 0.01            \nLin. dep.  - primal attempts        : 1                 successes              : 1               \nLin. dep.  - dual attempts          : 0                 successes              : 0               \nLin. dep.  - primal deps.           : 0                 dual deps.             : 0               \nPresolve terminated. Time: 0.03    \nOptimizer  - threads                : 8               \nOptimizer  - solved problem         : the primal      \nOptimizer  - Constraints            : 495             \nOptimizer  - Cones                  : 0               \nOptimizer  - Scalar variables       : 0                 conic                  : 0               \nOptimizer  - Semi-definite variables: 4                 scalarized             : 364095          \nFactor     - setup time             : 0.14            \nFactor     - dense det. time        : 0.00              GP order time          : 0.00            \nFactor     - nonzeros before factor : 1.23e+05          after factor           : 1.23e+05        \nFactor     - dense dim.             : 0                 flops                  : 1.16e+10        \nITE PFEAS    DFEAS    GFEAS    PRSTATUS   POBJ              DOBJ              MU       TIME  \n0   1.0e+07  2.0e+00  7.4e+02  0.00e+00   -7.440000000e+02  0.000000000e+00   1.0e+00  0.54  \n1   1.1e+06  2.3e-01  2.5e+02  -1.00e+00  -2.716669792e+03  -1.964891699e+03  1.1e-01  1.21  \n2   2.2e+05  4.3e-02  1.1e+02  -1.00e+00  -1.599144393e+04  -1.520305153e+04  2.2e-02  1.87  \n3   3.5e+04  7.1e-03  4.4e+01  -9.97e-01  -1.019276766e+05  -1.009136508e+05  3.5e-03  2.53  \n4   3.5e+03  6.9e-04  1.3e+01  -9.68e-01  -8.869180261e+05  -8.839394454e+05  3.5e-04  3.26  \n5   8.9e+02  1.8e-04  4.7e+00  -7.00e-01  -2.032785577e+06  -2.027243247e+06  8.9e-05  3.98  \n6   1.3e+02  2.5e-05  1.5e+00  -4.77e-01  -3.301761001e+07  -3.299259770e+07  1.3e-05  4.63  \n7   2.1e+01  4.3e-06  5.8e-01  -9.50e-01  -2.006392744e+08  -2.005040730e+08  2.1e-06  5.35  \n8   1.2e+01  2.4e-06  3.2e-01  -6.65e-01  -1.358862633e+08  -1.357512269e+08  1.2e-06  6.17  \n9   2.5e+00  4.9e-07  5.4e-02  -1.15e-01  -4.164412815e+07  -4.155693709e+07  2.5e-07  6.95  \n10  3.5e-01  7.1e-08  2.8e-03  6.85e-01   -8.930599022e+06  -8.919212721e+06  3.5e-08  7.81  \n11  7.7e-02  1.5e-08  2.4e-04  1.15e+00   -2.281918625e+06  -2.280093496e+06  7.7e-09  8.52  \n12  1.6e-02  3.1e-09  2.2e-05  1.15e+00   -9.998561813e+05  -9.994943199e+05  1.6e-09  9.19  \n13  5.4e-03  1.1e-09  5.2e-06  8.74e-01   -7.562405202e+05  -7.560716505e+05  5.4e-10  9.89  \n14  1.5e-03  3.0e-10  8.4e-07  9.02e-01   -6.512060401e+05  -6.511474567e+05  1.5e-10  10.56 \n15  6.1e-04  1.2e-10  2.9e-07  6.33e-01   -6.184457240e+05  -6.184051307e+05  6.1e-11  11.20 \n16  1.1e-04  2.3e-11  2.7e-08  8.11e-01   -5.948578391e+05  -5.948476286e+05  1.1e-11  11.87 \n17  2.3e-05  4.6e-12  2.7e-09  8.49e-01   -5.882472211e+05  -5.882446761e+05  2.3e-12  12.52 \n18  1.6e-06  3.6e-13  4.9e-11  9.61e-01   -5.864753155e+05  -5.864751383e+05  1.6e-13  13.20 \n19  5.7e-09  1.1e-13  1.1e-14  1.00e+00   -5.863470493e+05  -5.863470487e+05  5.7e-16  13.91 \nOptimizer terminated. Time: 13.93   \n\n\nInterior-point solution summary\n  Problem status  : PRIMAL_AND_DUAL_FEASIBLE\n  Solution status : OPTIMAL\n  Primal.  obj: -5.8634704931e+05   nrm: 2e+05    Viol.  con: 6e-02    var: 0e+00    barvar: 0e+00  \n  Dual.    obj: -5.8634704867e+05   nrm: 5e+06    Viol.  con: 0e+00    var: 0e+00    barvar: 2e-08  \nOptimizer summary\n  Optimizer                 -                        time: 13.93   \n    Interior-point          - iterations : 19        time: 13.93   \n    Simplex                 - iterations : 0         time: 0.00    \n    Mixed integer           - relaxations: 0         time: 0.00    \n\n","truncated":false}}
%---
%[output:5b17b169]
%   data: {"dataType":"textualVariable","outputData":{"header":"struct with fields:","name":"ans","value":"    yalmipversion: '20250626'\n    matlabversion: '25.2.0.2998904 (R2025b)'\n       yalmiptime: 3.6547\n       solvertime: 14.4563\n             info: 'Successfully solved (MOSEK-SDP)'\n          problem: 0\n"}}
%---
%[output:7eb77596]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAjAAAAFRCAYAAABqsZcNAAAAAXNSR0IArs4c6QAAIABJREFUeF7svQt0VsW5N\/5wsSZcE0hdUBIlXMqHwP9UqDQs8fOIWg5UFKHxgJF7QbxwehoQFsKngCKKGHr5IyIWkeOSW9Fl\/QiVIu06DZBCEQopAQSiEMQe0iQQCsFy5FvPvGd25t3vvs47e969k2evxdK8s2fmmd\/85pnfPDN772Z\/\/\/vfrwNdhAAhQAgQAoQAIUAIRAiBZiRgItRbZCohQAgQAoQAIUAIMARIwBARCAFCgBAgBAgBQiByCJCAiVyXkcGEACFACBAChAAhQAKGOEAIEAKEACFACBACkUOABEzkuowMJgQIAUKAECAECAESMMQBQoAQIAQIAUKAEIgcAiRgItdlZDAhQAgQAoQAIUAIkIAhDhAChAAhQAgQAoRA5BAgARO5LiODvSJw5coVeOWVV6CgoAC6d+\/uNVtK7jt58iRMnz4dFi9eDAMHDkyJDVjp3r17YfLkyaz+YcOGwaJFiyA9Pd23PTU1NfDEE0\/AD3\/4Qxg9erTv\/FHK4AWzLVu2wHPPPcewLSwstG1eUVERbNu2DV5\/\/XVfnPXK9abUL1HiENkqhwAJGDncKFcEEMCJZd68eb4ng1Q0LQwCBifBZ599Fjp37uw4yXrBp6lMlF4xQwGDoiQrK8tW1HGRg\/j7FTBR4roX\/tA9hIAXBEjAeEGJ7okkAlFy6mEQMCpFh8qywkw+r+3kAmb8+PFw6NChhMgWj6BgW\/\/whz+QgAlzp5NtoUGABExouoIM8YuAGLrHvOLKla9meZli6J5POocPH2bJ5hUvFxOzZs2CtWvXAr8Py3j88cdZlALD\/Hh53WYx22NXJ05w69atg3PnzrHyzVsOfMXP67e6xw5HJ0zMWGIZa9asSdjO4tgNHz6cTcR2ONjd169fP3jttdcgMzPTMJPjzduMCQsXLozberK6x8o+pzZiuZheWloKbdu2hU2bNiX0vRk7My5if3vFjNeLURXcIkS7n3766bgtImzfjh07WIQG7xMjMGa+mvvcqs2cpz169IDf\/\/73jMNo+09+8hOYOXOmEQXCLSu0R8SSt8ttu8vveKX7CQHVCJCAUY0olacFAe5kRceLjvxXv\/qVMUFaRWD4RIjOnJ9F4KtjPmnwe7Ah\/DdxsuJ18onl9ttvd9xyMduF5eLEsW\/fPsNWcYK2K5+LF8zPz6Z4jQCYz1bwfNnZ2b7KEidTLjKc7MKJk7fH6j6ryJO5b\/k9eEaIn6ex6lsvbfR6FoWLDjy3Yra\/srLS6Dev+HOOoY0oUPPy8uIEGqbn5OTAmTNn4gSMVfle8BCFrjhGzOXx+3ibsN14dknkhZYBTZUQAhIIkICRAI2ypB4BK1FgtXo2n4ExCwfMY55Yv\/jiC3agVpww7cQKloeRA7vDrnYTnHnitpqU0DZxou7QoYPUwVi77SmzUPAyGdvhYK7D7j5zv1nhZz5X4mUr0GsbzWLVjslu\/ca54QUzMQKDgvjgwYMsCsQ5g+1duXIlTJo0CXbu3BknYOx4bsbNjJGVWEQ73ATR559\/LnWIOPUegSxoigiQgGmKvd4I2ixGLMzbDbx5ZqfuFDERJ4rq6uqEJ4JkBYwZah6y57\/z1bHdBGyecHh+q60Yu261EwDmsr1Mxnb3mPGxu89JQJi3ZPgWhhj1sdvW8NpGL8IXcbTrD7O48oKZWcDg3\/Pnz4cXXniBbSOh7SUlJSyK54SPeRvNvJ0linW7w8VO\/YLRJrzsxlMjcBvUhEaGAAmYRtahTak5VudBROdrJ2D4mRYzVlwUqBYw4hkFfvYF6xYfm8bJSZzUuG1WE455oncTM3aTYjICBs9SiI97e53Y7bbr+PkX7D88X2P1NJRZ\/IlixmsbvQoYqy1KMVrHn9SSETDf+ta3WPv4NhK2a\/DgwQxPczvM5194m71GYMxPlLkJUGyj+YxSU\/Ip1NZoIUACJlr9RdY6IMAnOC5i7ASM27tJrFbfshEYu60huy0k83tg3CZIXg4e\/rSbeLxGJ9zqstuCEH\/n54G8RGD4RI75xS04L48mm8+yeG2jVwETZAQGoy78MDEKweXLl8MzzzzDDjabBYzVlifipVrAiOKQDu+Sm40KAiRgotJTZKcrAmaRYXfQ0+rMijgh8DMwopiQFTB2K3k+AZu3kMRzN9hgL+c\/3M51OJ2vwcnKfGjYSeD5PQNjLku0lZ\/pMR+Ctjpobe58L+eWOH5iG70KGLd2JnMGBgUMj7hhtOn8+fOWB8rtBB63TTxoa3cGxksEhnMUhT8eJBbxch10dAMhkEIESMCkEHyqWh4Bq0nb7ukVUYhYTY5e8skKGKv6rJ44snryyZzXyga7w5pmZL08oeMnAiM+XWQ1ofqJwIhP9YjbgjwSYCUCrUSZlzZ6FTCIn1lkmp\/YwYiJF8x4WeLj0XZPCVlFYMxv5uXREvEMjDli5PUMjFXf2UV95Ecr5SQEgkGABEwwuFKpGhAQz5bw6sRHRsVJQnT25nMFdu9kURGBQbvMhy\/xzMqLL77Itg149IHfY34PjPlApdl2LF\/2XTR2ZXuJwPzzP\/+z8X4RtMG87eBFwGAkwuocE\/YhHmoVHzO3eueK1WFTMyfM9\/gRMGIEh\/PLazvN9LcS3FZCwe4+xIRf2Ca87AQR8mHu3LmwZMmShLcqi\/3CzxrhU1FW753BOug8jAZHRlVII0ACRho6ykgIND0EvEYcmh4y1GJCgBDQjQAJGN2IU32EQIQRIAET4c4j0wmBRoYACZhG1qHUHEIgSARIwASJLpVNCBACfhAgAeMHLbqXECAECAFCgBAgBEKBAAmYUHQDGUEIEAKEACFACBACfhAgAeMHLbqXECAECAFCgBAgBEKBAAmYUHQDGUEIEAKEACFACBACfhAgAeMHLbqXECAECAFCgBAgBEKBAAmYUHQDGUEIEAKEACFACBACfhAgAeMHLbq3ySJw9OhR9mr5Xr16wciRI6Fly5ZNFgtqeDQQIM5Go5\/ISnkESMDIY0c5mwgCdXV1sHHjRhg3bhzs3r0bWrduDQMHDmwiradmRhEB4mwUe41s9otASgQMrQz8dhPdn0oE8DtFpaWlUFBQwL5rdOLECRg6dGgqTaK6CQFHBIizRJCmgIB2AeNnZbDpk\/NQWlEX64dmAPm3ZcGg3HZNoV+ojZoQwG2hzz\/\/HAoLC1mN4scF8aOL+DG76upqOHbsGODH79wEDHFWU8c14WqIs02486npcQhoFzBeVwY4EVTWXoXCIdmGwUU7K6Gy9isoGtWNupEQSBoB\/BowfuVX\/MKwODlg+i233AJ33XUXbN26FfBL0YcOHYKqqiq45557EuonzibdJVSACwLEWaIIIdCAgFIBo2plsKfiIuA\/Ubxwk1HEYBSGIjFE42QQQK7m5ORASUkJKwYjMDz6kpeXB6NHj4a9e\/fC5s2bYeHChVBcXMyEy6VLl2DixImQlZUVVz1xNpneoLxeECDOekGJ7mlKCCgTMCpXBoVbTkHRaPsoS+F7pygK05RYGmBbkbdmAZOfn88O6aKAWb58OdtGyszMdLSCOBtgJ1HRcQgQZ4kQhEAMASUCRuXKoG\/fvtQ3hEDSCJw6dQouX77sWo6KyYA46woz3eARgbKyMtc7ibOuENENGhHwwtmgzFEiYLhxqgZW\/39bAz8e0IIVm52dDQN+8Rl88FDsvRvN23WCuVvPwuI7W7C0yspKAxs\/f\/u5l9vB6\/KT18+9ydSTqrxhbd+ECRPAy8Cy4qx5C2nRokWQnp5uOwZRwNw5dyMsuDe2rZSbmwsdCn8H+2d0ZX+3aN8Zprx1AN4Y1YmlVVRUGGX5+dvPvdwOXpefvH7uTaaeVOUNa\/vwXBVxtmF8+OknP\/eminfJ1BvW9nnlbJMSMEvXbYX9lfUw7XsZCZPBpuPNoHvrKzCgSxpNBg6ToS7C66rH7+D3OrBEAYN1WB3ixfMwThcKmGQ4+\/Piw7D\/bD2rok2bNnB3DjB++22zU1\/46Sc\/96bKxmTqDWv7iLMx8U+iO+ZtZLHQmdcrZyMrYGRWsxj+f3J9OdzcIR3mDO3KVrMH5w+C9fu+ZDjgb\/xCsmNH8wvvKTlRww5b4mQw9vbOMLhHRhx+5jxios40nXVhG+3qC4sdKm3s1q2b79Us1i8+Rj1s2DBwi75gHhQwspxFvp6urmec5v3w8kefwenqK7BibG9bnhNn48d9qvAgzsb6gTjbwEDys0HJlcRyA9tCwqpkV7M4GeBVcqIWdp2shZc\/qoCxt3eCwT0y2X\/tnBVNBvYiRaWjtROPQU4ifkWWVwGjYqhxAeOXs5zfXJCLbcQJ4Y7uGYb49tv+oPqI7EhkjKoJizgbw1Y1x1SXJ2tjWOxQORfo5KyVrw5UwMiuZj\/++GPDVvN5AruwGobgT\/49HR7+9vW4vEiaN\/5YCwOy0+CHg3tLhSexQNlwXljD1Y19u0FnaBMFjAxnF+yogren3m57JsacrotLuupJZlwlkzes7SPOkp+1Ox9HnLVeaioVMKpXs7w83EKqLrrbsniuanHLySnkLqaHRQmTHY1vNauSsw+sOAC\/fvK2QFalUV9FhmXsRHU1K0YNibMNRxDCFEUOKooaVc4GHoEJs4ChyYDOwKjgp1UZqiYDUajjuZin1peTgDGdcUuFU9cl9nSG44mzwWxXqRQHqeC634WBTs42KQFDk0Fw0Q0aWPHYyk4G5jMwImfpDIzzBOPX0RJnibNWE6BOHumsqzGK7sgIGJnzBOYzMPjuGNx2ojMwsW73s4ca1L1mO\/zU4zdvFM4TYJuWldZDu2axVwYgZ\/Fpu5Xby1mf8dcIyDxW6hcv2TNeydSTqrxB8i6Zs2XEWToDQ2dg\/MXJQ3kGRkbA0GQgP\/iTmUhoMog9Ri3LWcTvVyXl7L1Hq\/bWsqfserW\/BiN6t0kQnrqw1lVPMrxLJm9Y2xcVAUOclXtXDXHWnzjxcncoBQx\/jJo3wMshXn6v30evzSDpDPPprMsppBgWO1TaqHNvVnYLycw9PzwX8+rsP511qeQDx0un\/X7rIs7Geskvbm59q7o8WRvDYofKcaWTs1aCJpQCJpnVLA\/BiVtIvOGpCJOHdbWXTKg7FTj6Xb1EaTVLnI19WkHXWNFVD3GW+pT8rJc4ivw9JGD+Bzs\/Ti2oe\/06PFWDI5l6\/WCRTD1+85KAkd9SDGufquJ7WNtHnCXO0hkYf2ImlAImmS0kmW0nCsc33dey+xsu1nfTFlIwoX+VoW63bQbZbQGVNuoMx6eKs6n41EtYtm7CYkdUORuZLaRUCBgaWHpe5uR3opDtF5oMghEVYXHCYbEjqpNBkALGbsy6fepFdqy75VPJFae6nNpt930+mfKQc7L57IS8bHk6\/WxkBIyXMzAfll+CYxdaso824oVPbYifCrA6A1NalQ6\/OXiG3S9+9RfLutKibdxnCPCr10dOn4cF92aBua6pd\/eELi1qDDztQtJ2Nsp+fVi0HysX7XCy0a8d2B4nG+1wRJtkbTRvEfGyDp08xx4n5ulzNxyAL+qusVfwO9nYWMLxYeGsXWjb79aeOFYaI2edxo7buAojZ\/36WbMvxTbjmP3T2Xr4bpc0WDLmtoRPZ2A6v8R0Jx\/M\/T36gHMXr8WVi5+OuXg9DWblxb7obuVbnM7x2fkwbJu5Lm4jfj3eyo6\/\/A2gT8fYKxHMPuxfvpMDdv4NyxPnJI4j+j5zGpbrZIcbjv+4MdPWDicbcW7UydlGI2C8kNYsYMIysJwIjYSQmbCsBk\/UB9a32raEjMxM229bYVptTY2lY9A9sLw+Rh3lyQAf7cZ3LZnFv1m0igsDJ0FLk0HsvAcX5GHkrF8\/i\/zA1wGIIoQLhWnvfQmPfS8jbpHJJ3R8jcCqP9bCG6M6JXxzzk7coEgRfal5EYn+D\/0DTvh+OGu3mN1Rdo4JMKt3M4l9KNqBeKz7cz2M\/6c0JjpEIeUksjCNL9LM74CySjMLI+SSWaDZ4YjfW\/vzX6\/B+wUNH0nmed1sxO8Lzh7\/AygrK7PSFlp+i9wZGK9f8BUfSXXKg58YmD001\/Krv5hv6UcVlq9zxzel8gu\/JszDlPj7rhM1cEePTDB\/ZZiXx+sTQ5sy5WH95s\/Ym8s8XX3F+EYUT3OzwyqPbF1e8lnhiP2S0yEtwXa8161f8CvO47\/fX9vA8hKOdwufY7vCyFm0S+wLkefIE\/xCPH72wDwGnNLsxocTV6LAWfze2q6TtexlhHh5GYu8zWHjrIyfdfoeHabhxb9XJ2KDaWeq6y39rNNYd\/LdON7W7z3HyvTKWSdeOtVl1zaOB\/+vuc1OeKhOs8MRbbu5Q7rx5XuvNqJ9mHfbPBIwcUrNbTLw+tFGcTIIy8ByIrTsQHVzDFYDoTEMLDuHl4qB5cbZKE8G3HYuav1MynZCmCaDxHed6J4M3DgbhJ8Vx6zon5EPeFl9uNRJ3Nj5Ps5ZFCSigHFbTMnwkvsbK3\/EMeTf4fMq2tzabOf7ZPK5zQVOfhbbVfbqaG0LRauQTigjME5nYPDV63xv0+rsyX2v7mGhSHELCcNkVmE1JBSmfXHxGvx25iCm1DHf\/hldWQgOy8LLKrQ5YfU+ls+chvejjXg2x1wn1oVnN7iNPPRnZYc5Df\/GvOazCHZ2oP12aW521HzVEn42PLZn68dGEUcxbzI2Yl7EUexPDEenp6VZ2oj1YvrJFQXaBpbbFpIbZ63aGBXOmnnJxyTa78RZTLMak42Js+axYzeuoshZL35WHLO4NYP\/7Pwsbi\/hNgvyB8WNmw\/mfvb4Fw1bTzwfcu\/50f1gylsHmH9286VeOSv6e9FGM2ft7DDPSSIfzGnmOYl\/FsdpvrKzw2kuM88FXm1E+zCvTj8bGQFjZSj9Rgj4RUDX3iwKGLoIARUIEGdVoEhl6ERAF2cjIWDcgN9TcRHwX+GQ7IRbi3ZWwqDcduyfeLnl2VNRB4Ny2xplVtZehU2fnIfK2q8gO+MblnXl\/7IcCod0SagL68W8mw9UweYpvRPsKNp51jKfTHlYeOF7p2xttEtDPOzskClP1g7M59RutGXPqYuwZ9Z3GI5e+sWOA268CjId21E0qptlFdimn2w5ZckVJ57r4izng9W4wrQzNVcTbOd8sEuzGx9OPIoCZ9FGHPdW7XYaV2HkrJvPtOKDVf\/xMcvShmQz343\/0Lfilde1LTzc\/5sJfszrWLcaBzieRH8uDjwnzjrxEvvIbszFfFP8XMHtl8njVJ5smtNcNmjZQRjUrZ3ho8Q+s2sX788g\/aaXskO3heTFaLMzMA8SqzLc8ugYWDhAZCYe2cET9YGF7eaXV4cXloHlR0CHeTJAoeskUuwEvtuELePYZfLIOnynfI11MjD7TTefaeer7Xypm2\/364PtRBFy0s4PJMNLJwGG9dnZ74SHzjSn\/kQfZCUsndrl1p860iMpYGRBTfXA4qsNVQOVl+eGR5QHllWkjQ8M2f7UMbAay2SAHCv9rC4uyud1hea2epNx3jJ5ZMeHU77GOBlYjYuwjDG\/djj1TzK8TIXvUFmnXxxV1h1EWaEVMHv37oXJkyfHtXnNmjUwcOBAWxxOnjwJ06dPh3PnzrG8hYWF7N6amhqYOXMmYJmdO3eGWbNmwbJly9h9\/fr1g9deew0yMzONcq9cuQIvvPACVFdXw4svvsjStmzZAs899xy7Z9iwYfDggw\/CwoULYfr8ZbC9ZB\/sWv9Tltbplp6Q2aoljBkzBkaPHp1gK7cxKyuLlf3MM8\/A4cOHmR1FRUXwzPMvw77\/3BGXD+vBsvxiohKPV155BUaMGAHvvvsus9mMCfbL6dOnYdWqVXDw4EEDq969Y9tobniktc2E\/vmF8Nu1S+Hi2ePQpVtvWLvq5wyTbdu2KcEjiAEklunWP1bOY0D7i5Hg7KQfTYMWrTKg5wM\/hrL3fsr6yAtn3\/rg9\/DqvKc8j2PibNAsjS\/fjbNW1qjsI69+9vXXX\/ftV2Q56xcTlXgE7Wft5p0o+VmRk6EUMAgmihVRsHBSicJEbAiKjmeffRby8vLYZI+Co7S0FBYtWgS1tbUwbtw4NvHii7awDC4K1q1bxybdCxcuxI3V69evs79vvfVWePzxx2HJkiWAg6h79+5MeGzfvh1WrFgBaWlpMH\/+fCZ4MA3Lw38jR46EXbt2xYkjbmN+fj4rW7QD7V27di088sgjMHbsWCZWSkpKWN3YrsrKSiZ0vGKiGg+0FzHJyMiAd955B86fPw\/z5s1jmHTo0IEJFDztPmPGDIZF2PAIelogzhJnVY9h4iz5WfKzzqMgdAIGoyU8MiFGRbAZ+Ajxo48+miA2zJMrTqR4cRHw0EMPwdSpU2H16tXQqVMnWLlyJUyaNIlFEbC+uXPnwqeffgpPPPEEEz\/mCMzOnTtZeTyignbw8r7++mvYsWMHPPbYY+weMQ0neyyTC4+ePXsabcN7xXaK+dB+tAHVeEFBATRv3tyoj7cN86Pt06ZNg\/Ly2EuixEsUGyrwMK8MREzQjgULFrDqJ0yYAPv27UsJHhwTO\/4ENSEQZ1cz8UqcbfALXvEgzjZEusnPxny+GOkmPxsxAcMjBxhxELd2cJJAMXDTTTcxRzlnzhwW8eCX6DzF31HEbNiwAQ4dOsQiMEOGDIlDBNM3b97MIjXFxcUsooCracxz5MgRJmywrsuXLxuTMoYM8f76+np2P5IMoz33338\/sxE\/cIVRGFGA8RV6+\/bt4fnnn0+IBGF9b775JkyZMsWIwPAIx9WrV9kWWLt27Vh9vFwuYFClYzSoT58+geHx\/vvvs4G1dOlSJsjMmKAoXL9+PWzcuBF2796dMjwQ\/+zsbNY\/6enpQWmWuHKJs8RZXNzIjmHiLPlZPu+Qn\/XnskMXgeHmi2dO+G\/iWRAuOsRJiu9FLl68OO6sDIoUXJXfcMMN8NJLLxnCB3\/nIoGLHj4Z4TkOvj2CDqaqqor9jRffMsKtEx5hwd\/xfA1O8OI5EXN38K0wPD9gPgMj\/o35cLuoY8eORn3i2RIRk5ycHEOEBY0H4sTFpBmT2bNnGwInlXhYnT3yNyzk7ibOEmdlxzBxNuZLyc\/G5h3ys958cGgFjDfz6S5CgBAgBAgBQoAQaIoIkIBpir3eCNrsdPJf9xmYRgAnNUEDAsRZDSBTFUoRCDtnScAo7W69hfGtHDyTYnVZPSKu18JganN6wgrPJZGACQZ3FaUSZxOfkiTOqmBWcGUQZ8PLWRIwwfFeWclOAwifNmrbti078Ny\/f3+jTrc8rVq1YoehzReW5zeNP3LerFkzJeVhIW528Me5zU9Y4dNlJGCUUU+6IDf+EWcbnpIkzkrTTGlG4mw8nE5PsoaFsyRglA6B4Arjj3s\/\/fTTcU9fYY3ik1TiIV6nPKrTVJeH7bIr0+mJMzzcjRd\/2V5wPUIluyFAnPX2lCRx1o1J+tKJs9HiLAkYfWMj6ZrshIpTwU55VKepLs9JnDk9cbZ8+fKEtysnDT4VIIUAcbbhUX7irBSFtGcizkaHsyRgtA8PqpAQIAQIAUKAECAEkkWABEyyCFJ+QoAQIAQIAUKAENCOAAkY7ZBThYQAIUAIEAKEACGQLAIkYJJFkPITAoQAIUAIEAKEgHYESMBoh5wqJAQIAUKAECAECIFkESABkyyClJ8QIAQIAUKAECAEtCNAAkY75FQhIUAIEAKEACFACCSLQEoEzNGjRwG\/3NurVy8YOXIktGzZMtl2UH5CgBAgBAgBQoAQaEIIaBcwdXV1sHHjRhg3bhzs3r0bWrduDQMHDmxCkFNTCQFCgBAIHgFaKAaPMdWQWgS0Cxh8G2VpaSkUFBQA\/v+JEydg6NChvlHY9Ml5KK2oi+VrBpB\/WxYMym3nuxzKQAh4QUDFZECc9YI03aMCAVULReKsit6gMoJCQKmAwW2hzz\/\/HAoLC5m9\/KvB27ZtA\/5l5Orqajh27BgMHz5cWsDgoKqsvQqFQ7INXIp2VkJl7VdQNKpbUFhRuU0UARWTAXG2iZInRc1WsVAkzqao86hazwgoEzBFRUWwZs0amDx5siFgREGD6bfccgvcddddsHXrVhg\/fjwcOnQIqqqq4J577vFs8J6Ki4D\/RPHCM6OIwSgMRWI8w0k3ekAg2cmAOOsBZLrFMwI6ForEWc\/dQTemEAElAgYHVE5ODpSUlLCmYASGR1\/y8vJg9OjRxheTFy5cCMXFxUy4XLp0CSZOnAhZWVmeISjccgqKRttHWQrfO0VRGM9o0o06JgPiLPFMFQK6ForEWVU9RuUEiYASAWNEQIqKEgRMfn4+O6SLX\/j08pXgvn37BtleKruJIHDq1Cm4fPmyY2tVTQbE2SZCqoCb6cZZlQtF4mzAndlEinfjbNAwhFLAvP3220a7s7OzobKy0vj7zaNp8KP\/Vc\/+xrQBv\/gMPngo9hh283adYO7Ws7D4zhZGuphXLMtcrp+\/g7qXt4nb7KeeVOUNq40TJkyAsrIy2\/GjejL4+OOPjbpyc3OhoqLC+HtZaT3Myktjf2Nah8Lfwf4ZXdnfLdp3hilvHYA3RnUy0sW8Ylnmcv38HdS9vE3cZj\/1pCpvWG3ErXQnzqpeKJKfjc0rYfVhqZiv\/M4jbn428gLGvIW0aNEiSE9Pt20XrgycJoOz\/50JW\/9UAdO+l5EwGWw63gy6t74CA7o0TBY0GcSg9uO0g7pX54QVpsngr82y4HdH\/gvG9G6eILqLz7aBm2+sg75ZzVydqS5Hq6sev87SyS4\/Nvu5V6eNXicDjBziJW7Vy0S6yc\/GFhl+\/J1OH5aKxYvf9nn1s0EJmcAiMGiw1SFePA\/jdKGAwbCUeKEIwc7k15Pry+HmDukwZ2hXtpo9OH8QrN\/3JUvG3\/Ay53EqL1VpYbExLHY49ZtfG7t166Z1NUuc1Tvm\/PKBj3Gd+fzWlQxnZRaKxFnirHkuDoqzkRQw4mPUw4YNA7foCzbSi4DB+0pO1MKuk7Xw8kcVMPb2TjC4Ryb7bypcWwg9AAAgAElEQVQclaxg8kuWoNoWFjvCImCCmAyIs\/GLkGQXDcTZWAQmyIUicZY46zZOvYruSAgYFUa6bSGZw2rieQI\/oUA\/92K7ZMN5uupJxsZk8oa1fV5Dm2I4PpnJwCkcT5xtOA\/kxBddXNJVj99xJctZ2YUicZa2kMwcDYqzKrSBVRlKt5BUGOk1AsPrQgFTXXR3QtWNcYVGEZiGbURVKwOzgJGdDNzC8aK9xNn44SozVmXyyEZKZfP5tVHnapb8bIyDfvsoGR+ssy7Ztvm1USdnScD8DwJ+OykZ0joRiexIpKQdJn6x0jmwaDKgyUA8oyfrL4izwfDIr++Q7T+3fGGxQ+WcpJOzJGBIwMQdhnYbcLIqPoh8fge\/zoFFAiaYiUelo00F14mzas+QyPoVv\/0QFFfCYofKcaXTz0ZGwNDeLO3NIlmTOS\/h9TyB1aDw+xud22p4701Yz5fInmEz8zDI9hFn6axhql\/74ZfvOjlLAkZ4wViQjigKzjIZcRCF9ukcWCRgSMBw55rMuCLOkoAhAeNv+UiHeE146Qzz6axLZdgwqBCrSht1hjZpC4m2kOgMjPfzbLJbQbL5yM967xu\/WOn0sxSBoQiMwYFkVooUgYkfShSBoQgMRWCsH5WnSHf8Z0VkfacuHGkLyV8EKOHuICeD0qp0+M3BM6zONm3awN054OmzA347VZakydSTqrw0sGIvXwzq3BZxVv07ZIizxNmo+WjirLWwaDJbSPipgdPV9exTAzxM9vJHn8Hp6iuwYmxvAx2\/IbRktlN01qVyeyaZNusKA+sMbQa1hUSctX9PR1jGjspxRZwNZisyLFwJix1R5WxktpBUvxSMf3bA6jtJKGLu6J4Bg3tkMHx0kkxnXSpJq1PA4CRecqIGLl26xKJmY2\/vbPSVkx1RnwyIs86TWVjGjspx1dQ56zTWZf2Abp+ukg86\/axsXTo5GxkBw8PxH5ZfgmMXWrLJC68RvdvADwf3ZiKDb5l4+ZTAgh1VsODeLNvPAWD621NvN8o1b8eY7Zh6d0\/o0qKG2eBko13az4sPw\/6z9Sy\/eSvLactATMO8TnZ4sdHODrdwpS4b\/3FjJhw6eY59eZz3ydwNB+CLumusv5xwTNUTHcTZ+K3ZpsZZ9E9O4yqMnDX7O\/Pf6A\/MftaqHch9Nz\/7L9\/JMbbxRR+Gec9dvAZLxtxm+OFNx5vBkdPn2Ta\/Oe2NP9bCxetpMCsvzfANdk\/weGmfirxuY98K1yDmAqt6ePtU26jTz0ZKwFgR2oq0KgbWtPe+hMn\/O1fpwOrV\/prloPvL3wD6dATLSRkH6pUWbeHhb19PmLDNaXiD0wB3G\/xOdjiJA+wXHTaiqPzzX6\/B+wUNH+jkwgp5kJGZCbU1NZY4ohPVObD4GRjibOwxWC4ymxpnkZdO4wonbztBnirOui3CrBaKKEqt2rFqby3sn9HVdqGIfvbevp0N\/8Z92I6yc\/DdLmlsLJvPpnAumYURz9u99RX44uI1y4VuEOLAbjHrNPZxLrBasNrhiAs02bkAcdJlo27ORkbArNv+CfvStJctH\/G7MnZnBvCL1eL3ksQ8eC7mgdcOsK0J8\/mYXSdq4I4embZ2iOdnxJC2+WwNT8NtgaUfVcDsoblsG8RLHuw0u\/Jk09zs4EQJGg8n+59cXw43d0g3tvdErLj9v37yNmaqGUfcEhz\/\/f5QVlZmxXnlv6GAIc4m9oPf8eHEhyhw1s1GJzxSwVncqvd6zor7TKetTRyz4hav6GedxuwDKw4YPtE8nrFMvPg5RfPWIeblPlo82yj6bvOZx8E9Mm3PQ8qkOc0TaF9OhzRmv519Vj5MduzI+G5ZG3VzNjICZtjirY4Ha5HUnNBeBhYKALy4IBIHljlNHCCyA8tu0HG7+X\/FupwGquo0JzuCcDQy9rthdaa6HqwEDPYz5t027wdaBQxxNl7AyPQ57zurCasxcNZpIk4FZ2VEt+h7zRMv+o71+84l+Ga8LxlfajfWsb6nNpTDwfmD4hYyyQhJGeFg1zYu9niZQc4tCICM75a1MRV+NjIC5s65G9leKg9fmvcn73t1D7wxqhMLN\/ItJLdzLAN+8RnMGZrLwpf4\/0j6ldvL4cOjl+DDCdmWr62fsHofs8FqXxfTMHTJ7RBtXFZaD8e\/qE1I4zZy+8X2YRqW99uZg+LO+GC5Vmk8r50diI1dmpsdNV+1hJ8Nj507MduIf4vnhXjIN2gbsc8wioZ4YCg6PS3N0ka0D9NPrijQKmCIs\/GfvwiaD2ZehpGzXm0MK2fRf3Df59XPjni7EsYNyvHtZ0Vfyuuy8mHikYEfF9fClfp6Wz\/L7ffjw8z+TfSzTnOBVZqbn8Xzj\/wMj465QKWNqeBsZASMlaH0GyHgFwGdW0h+baP7CQErBIizxIuoIaCLs5EQMGjknoqL7F\/hkOwEm4t2VsKg3Hbsn3gVvncKikZ1s+z7ytqr8JMtp2DzlIb3vfAbnfJt+uQ8bD5QZZsvO+MbljZimVZp2KainWehcEgXS\/v9lodtsKvLKc3Jjvxfllvah+WpxsPN\/kHLDsKgbu2MfsV+RBsqa7+yxBfLs+NH0E6BOBuPsAwvo85Z2fGdKs7K+EwZnnvx6Xsq6mBQblvDn\/KxjnnN452nbf6kCvbM+k7C0MZ2oY+18vmYdqbmqq1Pl0mz84t8DFjNV0H4UhnfLWtjqjhr7uzQvchOFBbihC4S2krYhGlgOQ06mYHqVJ5smp0dbuJAt\/04+LFv0S688rq2hYf7fzNBuLnxI2gBYzX5utlEnG1YAHiZsKLAWdnxYeXTguasLP\/M4tSN5159OrfHPNZFASSm5WTeaLnQlRWSsqIbJ3OrfsfFqp0gssuDwkbWp8v47mRsTAVnIyNg7EiLk5fdFZaBxW20G5B+B6pbeW5Y+bXDDUfd9tv1t5MdQTt\/VTa5YS3LddV9pLq8xsxZWaxSwVlZ\/smOPdl8fseBrJC0i\/a4iQqczK3aVvpZXVy02Cz2ZLlil8+tP1XbmArOinWGNgKzd+9emDx5chw+a9asgYEDB9pidvLkSZj0o2lQff6v0HXwKPj29yeyFft9uS1h5syZgGV27twZZs2aBcuWLYNz585Bv3794LXXXoOjtS2Mlf61r+qh6ve\/hPbNLsOLL74ImZmZsGXLFnjuuedY3cOGDYMHH3wQFi5cCK+\/\/jocPHjQSOvdO7ZNNWbMGBg9enSCrWjj9OnTISsri5X9zDPPwOHDh5kdRUVF7N+2bdvi8mE9WJZfTHhd2E7EsrCwkJVbU1PjG48Lu\/8D\/n3Sw\/Duu+8ym82YYL+cPn0aVq1aFVo8gh5sfvsH7SHOxvcKcTZolsaXb8XZ8XOXw+yC+xz9LPowGb+CfoNfV65cgRdeeAGqq6uT8rPf6n9fQpR2QPuLbC5o0SoDej7wYyh776dw8exxw88+8\/zLsO8\/d1j62bc++D28Ou8py7nHSgBgXW54tM68CbrfNxGOf\/QW1F84b8w7ZjxeeeUVGDFiRFJ+1g4Pp3nHCQ8Zv6aLxaEUMDiJo1gRBQsHUZyIRZBwMDz77LOQl5fHJnsUHKWlpbBo0SKora2FcePGsYkXT35jGVwUrFu3jk26Fy5ciMP8+vXYy+RuvfVWePzxx2HJkiVMrHTv3p0Ntu3bt8OKFSsgLS0N5s+fzwYipmF5+G\/kyJGwa9cuJo44SbmN+fn5rGzRDrR37dq18Mgjj8DYsWOZWCkpKWF1Y7sqKyuZ0PGKiWo80F7EJCMjA9555x04f\/48zJs3j2HSoUMHJtjwqYEZM2YwLMKGR9ADijhLnFU9homz5GfJzzqPgtAJGIwO8MiEqE6xGfgI7aOPPpogNsyTK06keHER8NBDD8HUqVNh9erV0KlTJ1i5ciVMmjSJCQusb+7cufDpp5\/CE088wcSPeWWwc+dOVh6PqKAdvLyvv\/4aduzYAY899hi7R0zDyR7L5MKjZ8+eRtvwXrGdYj60H21ANV5QUADNmzc36uNt45GUadOmQXl57GVP4iWKDRV4mFcGIiaI4YIFC1j1EyZMgH379qUED46JHX+CmhCIs6uZeCXONvgFr3gQZxsiMORnYz5fjMCQn42YgOGRA4w4iNELnCRQDNx0003MUc6ZM4dFPPglOk\/xdxQxGzZsgEOHDrEIzJAhQ+IQwfTNmzezSE1xcTGLKOBqGvMcOXKECRus6\/Lly8akjGFuvL++vp7djyTDaM\/999\/PbMQPXGEURhRgfIXevn17eP755xMiQVjfm2++CVOmTDEiMDzCcfXqVbbl065dO1YfLxcxQQGDKh2jQX369AkMj\/fff58NrKVLlzJBZsYEReH69eth48aNsHv37pThgfhnZ2ez\/klPTw9Ks8SVS5wlzuLiRnYME2fJz\/J5h\/ysP5cduggMN188c8J\/E8+CcNEhTlJ8\/3zx4sVxZ2VQpOCq\/IYbboCXXnrJED74OxcJXPTwyQjPtfDtEXQwVVVV7G+8+JYRbp3wCAv+judrcIIXz4mYu4NvheGZF\/MZGPFvzIfbRR07djTqE8\/aiJjk5OQYIixoPBAnLibNmMyePdsQOKnEw+rskb9hIXc3cZY4KzuGibMxX0p+NjbvkJ\/15oNDK2C8mU93EQKEACFACBAChEBTRIAETFPs9UbQZqenVXSfgWkEcFITNCBAnNUAMlWhFIGwc5YEjNLu1lsY38rBMylWF39E3HwYWq+V6mtzesIKzyWRgFGPuaoSibOJT0kSZ1WxK5hyiLPh5SwJmGA4r7RUpwGETxu1bduWHXju37+\/Ua9bnlatWrHD0OYLy\/Obxh85b9asmZLysBA3O\/jj3OYnrPDpMhIwSuknVZgb\/4izDU9JEmelKKY8E3E2HlKnJ1nDwlkSMMqHQTAF8se9n3766binr7A28Ukq8RCvUx7VaarLw3bZlen0xBke7saLv2wvmN6gUr0gQJz19pQkcdYLm\/TcQ5yNFmdJwOgZF0pqsRMqToU75VGdpro8J3Hm9MTZ8uXL4x7BVwI+FSKFAHG24VF+4qwUhbRnIs5Gh7MkYLQPD6qQECAECAFCgBAgBJJFgARMsghSfkKAECAECAFCgBDQjgAJGO2QU4WEACFACBAChAAhkCwCJGCSRZDyEwKEACFACBAChIB2BEjAaIecKiQECAFCgBAgBAiBZBEgAZMsgpSfECAECAFCgBAgBLQjkBIBc\/ToUcAP3\/Xq1QtGjhwJLVu21N5wqpAQIAQIgcaMAPnZxty71DZEQLuAqaurg40bN8K4ceNg9+7d0Lp167gvR1O3EAKEACFACCSHAPnZ5PCj3NFAQLuAwZc5lZaWQkFBAeD\/nzhxAoYOHRoNtMjKJosArWabbNdHsuHkZyPZbWS0TwSUChjcFvr888+hsLCQmcE\/urdt2zbgHxasrq6GY8eOwfDhw10FzKZPzkNpRV2sSc0A8m\/LgkG57Xw2kW4nBJJDwM9qljibHNaU2x0B8rPuGNEdTQMBZQKmqKgI1qxZA5MnTzYEjDjQMP2WW26Bu+66C7Zu3Qrjx4+HQ4cOQVVVFdxzzz0JaONEUFl7FQqHZBtpRTsrobL2Kyga1a1p9A61MhQIeF3NEmdD0V2N2gjys426e6lxPhFQImBQqOTk5EBJSQmrHiMwPPqSlxf7FDf\/vsTChQuhuLiYCZdLly7BxIkTISsrK87sPRUXAf+J4oXfgCIGozAUifHZ03S7JQKqVrPEWSJY0AiQnw0aYSo\/aggoETCGuCgqShAw+fn57JAuChivH9kr3HIKikbbR1kK3ztFUZioMS2E9qpczRJnQ9jBjdQk5K15oUh+tpF2NjXLEYHQCZi+fftSlxECShAoKyuzLUflapY4q6S7mnwhp06dgsuXL7vioELAEGddYaYbPCLg5Gc9FiF9W+ACxryFtGjRIkhPb\/hct9lyHFh3zt0IC+6NbSvl5uZCh8Lfwf4ZXdnfLdp3hilvHYA3RnViaRUVFUYRfv72cy+3g9flJ6+fe5OpJ1V5w9o+PFflZWCpmgz6\/9sa+PGAFoyL2dnZMOAXn8EHD8Xeb9S8XSeYu\/UsLL6zBUurrKw0OOvnbz\/3cjt4XX7y+rk3mXpSlTes7ZswYYI0Z8nPhn8uSMZHR93PSisUl4yBCRis1+oQL56HcbpQwCxdtxX2V9bDtO9lJAiYTcebQffWV2BAlzRLAfPz4sOw\/2w9q6JNmzZwdw6we5MhTzJ5w0o8J7v82Ozn3mRw9JtXt4BJhrNeRbgurHXV47dPibMxzymKbvKzsQUscdZewAWJjVc\/G0kBIz5GPWzYMHCLvmAjUcBgKPXJ9eVwc4d0mDO0K4vAHJw\/CNbv+5LhgL\/xC50\/dhBemH66up6l899f\/ugzOF19BVaM7W2ZxwysWF7QaTrrwrbY1RcWO1Ta2K1bN62rWVnOct6WnKhhh9pRdI+9vTMM7pERR7+w9BHZkeiKVY0rGc6iNeRnY32iqh+s5hYVc0FYxo5KrLxyNhICRoWRXMBgWSUnamHXyVp4+aMKGHt7JxjcI5P9V7w4Kfi9XNyIZEERc0f3DGNSCAuRyI7wTQayUUMUMH45S6LbeeJR6WiDmpRU2qhzMiA\/K889pz5XyQfirLuiULqF5F6d+x04sD7++GPjRvMZGDEcJv7\/gh1V8PbU223PxJjTgwyr2dmIjfJTr597zWXryqurHr\/t8xraNIfjZVezMpzFrc6Tf0+Hh799PY7vKGzf+GMtDMhOgx8O7m1wWhfWuurx26e0heTuP73eQX62YcslCnwPq41e\/axXXvq9L5QChq9meWNwC6m66G7LtvEoBm45OW0TPbDiAPz6ydtYGRT5iIcyLHioXL2kajWrkrMip8PSR2RH6qOGfp281f1iBEYlZ8nPhn8rK6p+1orHjVbAiKIHz8U8tb6cBIxwXkgkQ1gmpagOLJoMKBzPz+ElM66iKLrJzwYnaGkLyV2qNxoBYz4DIw4sOgPjPMGQgHEfKE53qBIwNBnQZODl0f\/k2BrLLctZ8rPyYp38rArmxpcRSgEjc54Am7WstB7aNYs9fo3v4cAnl1ZuL2ct5o9ky7zLBfPLnmsJ694lnSdQN5hkzxOYz8AgZ3GrlM7AxPrGz9gJ6t5kxr7fvDrPE8hylvys\/Fzglw+yc04y9fjNq5OzVh67UQkY7PBflZSzd8is2lvLnljq1f4ajOjdJsEh+nF4fjs1CsQjAZN6AUOTAU0G4juAdE4GyQgY8rP6DwDrmq\/8znU6ORsZASNziNfcOC8Hf60A0Rnm01kXtjXs70lQaWMUzhNw\/vl5XQBxVj6E78SvINL8ju8ocZbzkPxs4ogkP6tucehWUqOLwPAVjRiO5yCkIjISVuVMERi3oeE9PdnVLHFW79tUwzomda5mibP6oyh+oxupmK\/82qiTs5GJwMiegRE7nASMnm+D0GQQOxBJnNUrQvw6WpoM4t0\/cZYEjNXC3u+4IgFjklWyp+NpC0k+tO431M2x1pnPb10Ujpfng8qtvKC44pcPQdmhEiviLHE2alvEOjkbmQgMnYHx\/7I9\/A6U7Pd0aDLwvl1kdSeJbrmJR4azMnnEPtPJdb916ZwMiLNynE1GCPvlQzJ1OQnrqIruyAgYL+H4D8svwbELLdkH8PDCJ43E16773UIqrUqH3xw8Y2A09e6e0KVFDfvbXJdTmmiHnY1OX8wW7TB\/TdvORqzn3MVrsGTMbcZr5\/Gr3UdOn4cF92Yl2M9ttLPDbVvIyUanvDI4ytqoM7TpNRxPnAXgY0eGs\/+4MRMOnTzHXonAQ91zNxyAL+quOX5GBO8lzsa7\/1Rx1s23yGz1uY0r87aIGx9UzwWyPiwIzlphgZjL2qjTzzYaAYNgmyds\/HbMxetpMCsvjbXTSsDYEQIHwJUWbeO+ScMFwIAuabbiwCqN24GPb1vZ+Je\/AfTpGHsvjdkJY3miHUgs7qDNaZgXbdxRdg6+2yXN8j03mJdforhBG53swG9K2RHajJVoI4olGZFlh6OsjdhmnQPLy2RAnI0xUZaz+C2zP\/\/1Grxf0PAxVz7ZIZ9\/8N1c+OT4GcB36+Alin\/ibKLrTwVnzcIB\/052oeg0rtAH++WDnZ\/FxaDMXCDrw4LgrN1i1m1h4LTg1ulnIyNgnLaQvH512vx4H4ae8ZMC+LVqHsrDN\/TuOlEDd\/TIZL\/jZf6K9enqK8Y3lpJNQ9uXflQBs4fmsi9jJ1se2ovfHuHlme3n9Zm\/AeVmByeKX6xk81lhLGsjLytM4XjibPy4kuEsfhfq5g7pxlflzeF4LJOPYz\/jmzjbMC2IPlM1Z3FcDu6RaemDZdOcfDfyIadDGvPdfvig09\/Lck9lPhxXu07Wspe+Ws1\/dnWlws82CgHj9tFGnu51MDo5UywLL\/6RSNFpyqRx2\/h\/ky0PbXOz40x1fcI3oJzssBM9smLJLZ+d\/bI28s9GjP9+fwjLa9mJs\/ECJgjOPrWh3NIJy4gl4mzMrzh9HNevn8VxqVIcOPURF1+8PtHPhsXf6\/azdvU5LQycbEyFn42MgHE6A4OfC+DbRFb7qfe9ugfeGNUpbgsJw88YPrPaX8U0vHDbRHwrJt47YfU++OLiNVaeOa9dGpaFNh7\/ojYhH9aF9XAb8V5eLqZhXb+dOciwwymN5zXbgcJt\/4yurFysBy+z\/W521HzVEn42PLbFZbbRjJVXGzGfVR\/Y4ShrI6tnRxX8Ycm\/ahUwxNn4x6it+KySs+IW8bT3voT0tDTirJWHt\/nNbQtJt5918g9WadzP4hlIs19x8x2yPszvXOBmR5j8LPphxNHs7+1sTIWfjYyA8TEO6VZCwBYBnREY6gZCQAUCxFkVKFIZOhHQxdlICBg34PdUXAT8VzgkO+HWop2VMCi3HfsnXoXvnYKiUd0si970yXnYfKAKNk\/pnZCO+bIzvmFZl0wa2l208ywUDuliaaNMXdjmPRV1MCi3rWFnZe1VwHZV1n5lab+THfm\/LLe0D8Fxwko2nx2OsjZi23+y5ZRlf7pxK6h04mw8sjKcxRIGLTsIg7q1M8Yy5\/nmT6pgz6zv+B7fxFl7xqvmLI7zMzVXbf2sTJqdP+I+xWouCIu\/l+We6nzYzzj\/WeHvVFdY\/GzoPiXgZRIxT3rckWFeK2HjNhjtBACWZxYBYl0yaUHUhW3mbUSb8Mrr2hYe7v9NsMPKr+jBMmUnHqd8ThjL2mjluLzwKsh7iLM4lmLCmo9Tv5xlIia3HeO6yPOczBsdFzWyPJIZq7J1NXbOyiz4GE8cFpF2fgUXiXaCSNYXyc4FsnzQmc9uYWC3AOZzQRg4G0kBgwDaOT+7SchtAnEqT3Wa6vLcJl67+mQmEBmx5DZhOfWnrI1umKQinTjbIKzd8NeFVRB+gTjbELVOdsHnJhzs\/FHpZ3Vx0Wev4tltbpHx3bJ80JnPamHgtAC2Cxa4jWvV6aEVMHv37oXJkyfHtXfNmjUwcOBAWwxOnjwJ06dPh3PnzrG8hYWF7N6amhqYOXMmYJmtM2+C7vdNhOMfvQX1F85Dv3794LXXXoPMzEyj3CtXrsALL7wA1dXV8OKLL7K0LVu2wHPPPcfuGTZsGDz44IOwcOFCeP311+HgwYNGWu\/esa2oMWPGwOjRoxNs5TZmZWWxsp955hk4fPgws6OoqIj927ZtW1w+rAfL8ouJFzw6d+4Ms2bNgmXLljHcOB5Ha1sYK91rX9XDhd3\/Af8+6WF49913mc1mTLBfTp8+DY\/\/n1dhe8k+2LX+p6wNnW7pCZmtWgaCx1sf\/B5enfeUL46oHkBieX77B\/N66SPibPwYJs6qY7Fuzpr9StXvfwntm11mvhDT1r67yfAdffKGwIyJDzM\/O33+Mt9+ZdKPpkGLVhnQ84EfQ9l7P4WLZ483eT+b1jYT+ucXwm\/XLmV4dOnWG9au+rnjvBM2PyuyP5QCBidxFCuiYOEDTRQmYkNQdDz77LOQl5fHJnsUHKWlpbBo0SKora2FcePGsYkXX3CFZXBRsG7dOli1ahVcuHAhzitcv36d\/X3rrbfC448\/DkuWLGFipXv37mywbd++HVasWAFpaWkwf\/58JngwDcvDfyNHjoRdu3bFiSNuY35+PitbtAPtXbt2LTzyyCMwduxYJlZKSkpY3diuyspKJnS8YqIaD7QXMcnIyIB33nkHzp8\/D\/PmzWOYdOjQgQkUfCJpxowZDIuw4aHO5VuXRJwlzqoew8RZ8rPkZ51HQegEDEZLeGRCjIpgM\/Ax50cffTRBbJgnV5xI8eIi4KGHHoKpU6fC6tWroVOnTrBy5UqYNGkSiyJgfXPnzoVPP\/0UnnjiCSZ+zBGYnTt3svJ4RAXt4OV9\/fXXsGPHDnjsscfYPWIaTvZYJhcePXv2NNqG94rtFPOh\/WjDK6+8AgUFBdC8eXOjPt42zI+2T5s2DcrLY++rES9RbKjAA20ZMWKEEYERMUE7FixYwKqfMGEC7Nu3LyV4cEzs+BPUhECcXc3EK3G2wS94xYM42xDpJj8b8\/nkZ7176tAJGB45wIiDuLWDkwSKgZtuuok5yjlz5rCIB79E5yn+jiJmw4YNcOjQIRaBGTJkSBw6mL5582YWqSkuLmYRBVxNY54jR44wYYN1Xb582ZiUMeyP99fX17P7ceBhtOf+++9nNuJbYDEKIwowvkJv3749PP\/88wmRIKzvzTffhClTphgRGB7huHr1KtsCa9euHauPl8sFDKp0jAb16dMnMDzef\/99NrCWLl3KBJkZExSF69evh40bN8Lu3btThgfin52dzfonPT3d+0hI4k7iLHEWFzeyY5g4S36WzzvkZ\/054tAJGG6+eOaE\/yaeBeGiQ5yk+HmCxYsXx52VQZGCq\/IbbrgBXnrpJUP44O9cJHDRwycjPNfCt0fQwVRVVbG\/8eJbRrh1wiMs+DvuzeMEL54TMXcH3wrDsybmMzDi35gPt4s6duxo1AWmlu0AACAASURBVCeetRExycnJMURY0HggTlxMmjGZPXu2IXBSiYfV2SN\/w0LubuIscVZ2DBNnY76U\/Gxs3iE\/680Hh1bAeDOf7iIECAFCgBAgBAiBpogACZim2OuNoM1OT+\/oPgPTCOCkJmhAgDirAWSqQikCYecsCRil3a23ML6Vg2dSrC6rR8T1WhhMbU5PWOG5JBIwweCuolTibOJTksRZFcwKrgzibHg5SwImON4rK9lpAOHTRm3btmUHnvv372\/U6ZanVatW7DC0+cLy\/KbxR86bNWumpDwsxM0O\/ji3+QkrfLqMBIwy6kkX5MY\/4mzDU5LEWWmaKc1InI2H0+lJ1rBwlgSM0iEQXGH8ce+nn3467ukrrFF8kko8xOuUR3Wa6vKwXXZlOj1xhoe78eIv2wuuR6hkNwSIs96ekiTOujFJXzpxNlqcJQGjb2wkXZOdUHEq2CmP6jTV5TmJM6cnzpYvX57wduWkwacCpBAgzjY8yk+claKQ9kzE2ehwlgSM9uFBFRIChAAhQAgQAoRAsgiQgEkWQcpPCBAChAAhQAgQAtoRIAGjHXKqkBAgBAgBQoAQIASSRYAETLIIUn5CgBAgBAgBQoAQ0I4ACRjtkFOFhAAhQAgQAoQAIZAsAiRgkkWQ8hMChAAhQAgQAoSAdgRIwGiHnCokBAgBQoAQIAQIgWQRSImAOXr0KOCXe3v16gUjR46Eli1bJtsOyk8IEAKEACFACBACTQgB7QKmrq4ONm7cCOPGjYPdu3dD69atYeDAgU0IcmoqIUAIEALBI0ALxeAxphpSi4B2AYNvoywtLYWCggLA\/z9x4gQMHTrUNwqbPjkPpRV1sXzNAPJvy4JBue18l0MZCAEvCKiYDIizXpCme1QgoGqhSJxV0RtURlAIKBUwuC30+eefQ2FhIbOXfzV427ZtwL+MXF1dDceOHYPhw4dLCxgcVJW1V6FwSLaBS9HOSqis\/QqKRnULCisqt4kioGIyIM42UfKkqNkqForE2RR1HlXrGQFlAqaoqAjWrFkDkydPNgSMKGgw\/ZZbboG77roLtm7dCuPHj4dDhw5BVVUV3HPPPZ4N3lNxEfCfKF54ZhQxGIWhSIxnOOlGDwgkOxkQZz2ATLd4RkDHQpE467k76MYUIqBEwOCAysnJgZKSEtYUjMDw6EteXh6MHj3a+GLywoULobi4mAmXS5cuwcSJEyErK8szBIVbTkHRaPsoS+F7pygK4xlNulHHZECcJZ6pQkDXQpE4q6rHqJwgEVAiYIwISFFRgoDJz89nh3TxC59evhLct2\/fINtLZTcRBE6dOgWXL192bK2qyYA420RIFXAz3TircqFInA24M5tI8W6cDRqGUAqYt99+22h3dnY2VFZWGn+\/eTQNfvS\/6tnfmDbgF5\/BBw\/FHsNu3q4TzN16Fhbf2cJIF\/OKZZnL9fN3UPfyNnGb\/dSTqrxhtXHChAlQVlZmO35UTwYff\/yxUVdubi5UVFQYfy8rrYdZeWnsb0zrUPg72D+jK\/u7RfvOMOWtA\/DGqE5GuphXLMtcrp+\/g7qXt4nb7KeeVOUNq424le7EWdULRfKzsXklrD4sFfOV33nEzc9GXsCYt5AWLVoE6enptu3ClYHTZHD2vzNh658qYNr3MhImg03Hm0H31ldgQJeGyYImgxjUfpx2UPfqnLDCNBn8tVkW\/O7If8GY3s0TRHfx2TZw84110Dermasz1eVoddXj11k62eXHZj\/36rTR62SAkUO8xK16mUg3+dnYIsOPv9Ppw1KxePHbPq9+NighE1gEBg22OsSL52GcLhQwGJYSLxQh2Jn8enJ9OdzcIR3mDO3KVrMH5w+C9fu+ZMn4G17mPE7lpSotLDaGxQ6nfvNrY7du3bSuZomzesecXz7wMa4zn9+6kuGszEKROEucNc\/FQXE2kgJGfIx62LBh4BZ9wUZ6ETB4X8mJWth1shZe\/qgCxt7eCQb3yGT\/TYWjkhVMfskSVNvCYkdYBEwQkwFxNn4RkuyigTgbi8AEuVAkzhJn3capV9EdCQGjwki3LSRzWE08T+AnFOjnXmyXbDhPVz3J2JhM3rC2z2toUwzHJzMZOIXjibMN54Gc+KKLS7rq8TuuZDkru1AkztIWkpmjQXFWhTawKkPpFpIKI71GYHhdKGCqi+5OqLoxrtAoAtOwjahqZWAWMLKTgVs4XrSXOBs\/XGXGqkwe2UipbD6\/NupczZKfjXHQbx8l44N11iXbNr826uQsCZj\/QcBvJyVDWicikR2JlLTDxC9WOgcWTQY0GYhn9GT9BXE2GB759R2y\/eeWLyx2qJyTdHKWBAwJmLjD0G4DTlbFB5HP7+DXObBIwAQz8ah0tKngOnFW7RkSWb\/itx+C4kpY7FA5rnT62cgIGNqbpb1ZJGsy5yW8niewGhR+f6NzWw3vvQnr+RLZM2xmHgbZPuIsnTVM9Ws\/\/PJdJ2dJwAgvGAvSEUXBWSYjDqLQPp0DiwQMCRjuXJMZV8RZEjAkYPwtH+kQrwkvnWE+nXWpDBsGFWJVaaPO0CZtIdEWEp2B8X6eTXYrSDYf+VnvfeMXK51+liIwFIExOJDMSpEiMPFDiSIwFIGhCIz1o\/IU6Y7\/rIis79SFI20h+YsAJdwd5GRQWpUOvzl4htXZpk0buDsHPH12wG+nypI0mXpSlZcGVuzli0Gd2yLOqn+HDHGWOBs1H02ctRYWTWYLCT81cLq6nn1qgIfJXv7oMzhdfQVWjO1toOM3hJbMdorOulRuzyTTZl1hYJ2hzaC2kIiz9u\/pCMvYUTmuiLPBbEWGhSthsSOqnI3MFpLql4Lxzw5YfScJRcwd3TNgcI8Mho8VyXAiKTlRA5cuXWKRm7G3dzbuT2YyDwuhw2JHVAdWEAKGOOs8mRFnkwt1E2eDEUsqfVgyc0tjXChGRsDwcPyH5Zfg2IWWTDjgNaJ3G\/jh4N5MZPAtE\/OnBH5efBj2n62P2ybCchbcm2X7OYAFO6rgX76TY2wvYeapd\/eELi1qAPOeu3gNloy5zagXv3p95PR5VqYXG+1Olpu3ffBv2S0Dsx3cfizTzUa\/209ONoppIo5Wdog2WmGBYVOr\/sSvjbuFVKP0RAdxNja+ibP3ePoAaXLSJZZb3PZ08w841sLsZ+3st\/Mdbv7OyYe55XXzS+atq6jbqNPPRkrAWAmHN\/5YCxevp8GsvDRLAYPEO3TyHEz7XiyagmSZu+EArNpbC\/tndLUVMNPe+xLu7dsZHv72dQMjFCk7ys7Bd7uksfLMxMNy+SWKG9FGvwML77\/Soq1hB7f\/i7prTCzZDSwnkYWTvVmAcRt7tb+WIPa4OLAbWE42Yl7RfsSHiz0rO9yE4D9uzLTsT8Tj7am324obrFfnwEpmMiDOJi4MiLMqJIpzGZyzUfez6MOs\/Ntf\/gbQpyMkzAVWvlQ8D2n2b6IPs1qwehHddr7UbuxHyUadfjYyAmbd9k\/Yl6a9bPnw78o4hdyfXF8et+0jfosG8y39qAJ+\/eRtDB8xNP3AigMwe2iu5faSUz7cluKX+czNrhM1cEePzISzOOLvZju8lmeVTzzjY25bToc0dv7HfCYIv+xtdV7IzUa7usxnjUQ77NqGfYYcODh\/UEK\/OOHBbUjFeQKvZ1aIsw1jgDjbcA5PN2cbg5+18jncN3Pf7cXfoE9x8292Ph1t8OIzRT8rW5dsviBsxKMX47\/fX1vUMDICZtjirY4Ha3Fy4wdv+WQg\/mZ2ikjo9fvOJeTB+5xECpaJF69LHAiYdqa63lL42Ikbp4ElK5Zk7Odijw9+8wB3Ej12gs4NKzsc7bDC8m7ukG6cTxJtdBOPugcWrmZlJgPibOK5MzseEWfVRmWQs1H3s3Zc4eOK\/9er75D1wXYLNN3+3s0v+hV7TnggGxHfbfN+QAJGHJo4sO6cu5FtmfBtIPMZkvte3QNvjOrEtnX43iyeY8FtBbvzJiPeroRxg3LY9syAX3zGVvYrt5fDh0cvwYcTsi3PVExYvQ++uHgtoS6sF23Ai9sh1rustB6Of1GbkMZtxHLN7cM0vMQ28G0rTEM7fjtzUEL7sCy8rM74mO3nNnI7OI4izpiGIVW+TWdOk7VRxFE8w+Rmox1WNV+1hJ8Nj02AZp5gG\/6w5F+1DSwvnOXtIM7mAnE2GpyNup+1829mXyr6WbN\/477Fzk9hXjs+u\/l7Jz8bdhsRFzx6cXJFgTY\/G5kIjNq1BpXWVBEoKyvT0nQUMHQRAioQIM6qQJHK0ImALs5GQsCgkXsqLrJ\/hUOyE2wu2lkJg3LbsX\/iJZPHS117KupgUG5bw5bK2quw6ZPzUFn7FWRnfMPSxvxflkPhkC4JNha+d4r9\/pMtp2DzlIZ3z6AdWObmA1UJv2OaXXlu+bA+Kxv571Y4YtqZmquWdjjZaFcX2uiUZtc27E\/Ew8oWJzywf6zwDXpAYxuLRnWzrMbOJuJs\/Bhw4gpxVj2DZfgnkycoP2vnV9DGop1nLX2wrC91ymfnM3X7e5022s3D6lnqXGLoXmTHzTWTkwsH5uQshI2V8\/OSx0s+PmhRtOCV17UtPNz\/mwkTs5u4cRpYSAgZseSUD201Cy20EQe3nUhxEhsydaHosbPDTQgOWnYQBnVrZwgDN3yxnlQNLFnHLsNz4my8UyPOyk8bMvyTyRMEZ538impf6rRglRFSQfh7nTY6zcPybPSfM7QCRlTtZuHg1Ew7seEGjep8doPcbmBxQvgVS2757HAs\/awuLjojij0r0cOFCIpHJ6xk0pwcIkaJsEyv4jGVA0vWsavmnizXibPxCxQnH9RYOBt1P+tkv2pf6tTndj5Tt7\/XaaObn9GRHloBs3fvXpg8eXIcBmvWrIGBAwfa4nLy5EmYPn06nDt3juUtLCxk99bU1MDMmTMBy+zcuTPMmjULli1bxu7r168fvPbaa5CZmWmUe+XKFXjhhReguroaXnzxRZa2ZcsWeO6559g9w4YNgwcffBAWLlwIr7\/+Ohw8eNBI6907FhYfM2YMfKv\/fQmT74D2F2HSj6ZBi1YZ0POBH0PZez+Fi2ePMzuKiorYv23btsW1EesZPXo0vPXB7+HVeU95xsQLHq0zb4Lu902E4x+9BfUXzht4HK1tYdh+7at6uLD7P+DfJz0M7777LjzzzDMJmGC\/nD59GlatWmWLB7bBfHEb09pmQv\/8Qvjt2qUMjy7desPaVT9XikfQA8qKs+PnLofZBfcRZwUEnMYxcTZolsaX3xj8rJNfycrKYj4cfdbhw4cNP\/vM8y\/Dvv\/cYelnnTCxEkXo0\/m80+\/eh6Hz\/36Ulfv\/dfwa\/vPtF9m8Y+dnzfPOK6+8AiNGjGB+9r6CGVBW3Rz27Pi\/UP7r\/5+V6cXP2s07aKOdn5XFQy9bE2sLpYDBSRydnOjoOKlEYSI2B0XHs88+C3l5eWyyR8FRWloKixYtgtraWhg3bhwjMZ78xjK4KFi3bh2bdC9cuBCHzvXrsZfa3XrrrfD444\/DkiVLmFjp3r07GxDbt2+HFStWQFpaGsyfP58JHkzD8vDfyJEjYdeuXXHiiNuYn5\/PyhbtQHvXrl0LjzzyCIwdO5aRvqSkhNWN7aqsrGQD0CsmqvFAexGTjIwMeOedd+D8+fMwb948hkmHDh2YYMMT+TNmzGBYhA2PoAcacZY4q3oME2fJz5KfdR4FoRMwGC3hillUp9gMfPz20UcfTRAb5skVJ1K8uAh46KGHYOrUqbB69Wro1KkTrFy5EiZNmsSiCFjf3Llz4dNPP4UnnniCiR9zBGbnzp2sPK700Q5e3tdffw07duyAxx57jN0jpuFkj2Vy4dGzZ0+jbXiv2E4xH9qPNqAaLygogObNmxv18bZhfrR92rRpUF4ee1+NeIliQwUe4soAhaCICdqxYMECVv2ECRNg3759KcGDY2LHn6AmBOLsaiZeibMNfsErHsTZhkg3+dmYz+cRGPKz7h47dAKGRw4w4iBu7eAkgWLgpptuYo5yzpw5LOLBL9F5ir+jiNmwYQMcOnSIRWCGDBkShwqmb968mUVqiouLWUQBV9OY58iRI0zYYF2XL182JmUMc+P99fX17H4ceBjtuf\/++5mN+EZNjMKIAoyv0Nu3bw\/PP\/98QiQI63vzzTdhypQpRgSGRziuXr3KtsDatWvH6uPlcgGDKh2jQX369AkMj\/fff58NrKVLlzJBZsYEReH69eth48aNsHv37pThgfhnZ2ez\/klPT3cfAQruIM4SZ3FxIzuGibPkZ\/m8Q37Wn0MOnYDh5otnTvhvfNtHFB3iJMX3zxcvXhx3Vgbvx1X5DTfcAC+99JIhfPB3LhK46OGTEZ5r4dsj6GCqqqrY33jxLSPcOuERFvwdz9fgBC+eEzF3B98KwzMv5r1Z8W\/Mh9tFHTt2NOoTz9qImOTk5BgiLGg8ECcuJs2YzJ492xA4qcTDak\/c37CQu5s4S5yVHcPE2ZgvJT8bm3fIz3rzwaEVMN7Mp7sIAUKAECAECAFCoCkiQAKmKfZ6I2iz09Mqus\/ANAI4qQkaECDOagCZqlCKQNg5SwJGaXfrLYxv5eCZFKvL6hFxvRYGU5vTE1Z4LokETDC4qyiVOJv4lCRxVgWzgiuDOBtezpKACY73ykp2GkD4tFHbtm3Zgef+\/fsbdbrladWqFTsMbb6wPL9p\/JHzZs2aKSkPC3Gzgz\/ObX7CCp8uIwGjjHrSBbnxjzjb8JQkcVaaZkozEmfj4XR6kjUsnCUBo3QIBFcYf9z76aefjnv6Cmu0O9TslEd1murysF12ZTo9cYZPlOHFX7YXXI9QyW4IEGe9PSVJnHVjkr504my0OEsCRt\/YSLomO6HiVLBTHtVpqstzEmdOT5wtX7484e3KSYNPBUghQJxteJSfOCtFIe2ZiLPR4SwJGO3DgyokBAgBQoAQIAQIgWQRIAGTLIKUnxAgBAgBQoAQIAS0I0ACRjvkVCEhQAgQAoQAIUAIJIsACZhkEaT8hAAhQAgQAoQAIaAdARIw2iGnCgkBQoAQIAQIAUIgWQRIwCSLIOUnBAgBQoAQIAQIAe0IpETAHD16FPDDd7169YKRI0dCy5YttTecKiQE\/CBAnPWDFt1LCBAChEDwCGgXMHV1dbBx40YYN24c7N69G1q3bh335ejgm0w1EAL+ECDO+sOL7iYECAFCQAcC2gUMvsyptLQUCgoKAP\/\/xIkTMHToUB1tpToIASkEiLNSsFGmFCNAUcMUdwBVHzgCSgUMbgt9\/vnnUFhYyAznH93btm0b8A8LVldXw7Fjx2D48OGuAmbTJ+ehtKIuBkIzgPzbsmBQbrvAQaEKmg4CxNmm09dNqaV+oobkZ5sSMxpXW5UJmKKiIlizZg1MnjzZEDDi5IDpt9xyC9x1112wdetWGD9+PBw6dAiqqqrgnnvuSUAVB1Vl7VUoHJJtpBXtrITK2q+gaFS3xtUL1JqUIECcTQnsVKkGBLxGDcnPaugMqiIwBJQIGBQqOTk5UFJSwgzFCAyPvuTlxT7Fzb8vsXDhQiguLmbC5dKlSzBx4kTIysqKa+CeiouA\/0Txwm9AEYNRGIrEBMaJJlEwcbZJdHOjbKSqqCH52UZJjybVKCUCxhAXRUUJAiY\/P58d0kUB4\/Uje4VbTkHRaPsoS+F7pygK06RoGlxjMQpjFt3E2eDwppKTQ0Bl1JD8bHJ9QblTj0DoBEzfvn1TjwpZ0CgQKCsrc22HCgFDnHWFmW7wiIATZ1VGDYmzHjuEbnNE4NSpU3D58uWUoRS4gDFvIS1atAjS0xs+121uOQ6sO+duhAX3xraVcnNzoUPh72D\/jK7s7xbtO8OUtw7AG6M6sbSKigqjCD9\/+7mX28Hr8pPXz73J1JOqvGFtH56rkhUwxNnwj6tk+N4YOes3aoh+tv+\/rYEfD2jB\/Gd2djYM+MVn8MFDsXdyNW\/XCeZuPQuL72zB0iorKw0\/6+dvP\/dyO3hdfvL6uTeZelKVN6ztmzBhgic\/G5TCCUzAoMFWh3jxPIzThQNr6bqtsL+yHqZ9LyNBwGw63gy6t74CA7qkWQqYnxcfhv1n61kVbdq0gbtzgN2bjMNLJm9YnaWTXX5s9nNvMjj6zSsjYIizscVAWPuUOBvznKqihsn4Wa8LR11c0lWPXz9EnA1KusTKDVTAiI9RDxs2DNyiL2gQChgMSz25vhxu7pAOc4Z2ZRGYg\/MHwfp9XzKj8Td+4UBCkuCF6aer61k6\/\/3ljz6D09VXYMXY3pZ5zPCK5QWdprMubItdfWGxQ6WN3bp187QyECcDrJ84G2N92LlCnHV+WMLN1ybjZ7mvLTlRwx7EwIXi2Ns7w+AeGXEuMyx9RHYkighV49urnw1KxigVMCqM5AMLyyo5UQu7TtbCyx9VwNjbO8HgHpnsv+LFO4Lfy8WN2EEoYu7onmEMMCJ0fE+FBQ+VE6fOgUWcdRY9Tv0qm9aUOWsW3bKRblwo+vWztFCU53pT5qwKbWBVRigFzMcff2zYaj4DI4bkxP9fsKMK3p56u+2ZGHO6rpCjrnoQMDtszGluf\/ux2c+9Om30uoWkYmChgCHORmf7KeqcVRU1lOEsbs+f\/Hs6PPzt63E+GifnN\/5YCwOy0+CHg3sbflgX1rrq0enDUuHP\/bZPp5+NjIDhKwNuMG4hVRfdbTnXcFWLW05O20QPrDgAv37yNlZGWJQw2dE4QptiBIY4G9vONV+qua66vCD8gl8bUxU1VMlZ0Q\/7bT+3Q3U+1eXJciUsdjjZ79dGnZxtUgJGFD14Luap9eUkYITzQiIZ\/JI2KEcT1YGlSsAQZ70LH+JscrFDVZw19wMtFMN\/fiyqfrZRCxjzGRhxMqAzMLGuV3VwiwRM\/FCSnQyIs868JM4mJ1Kccsty1hzpJtFNotvL6yqCYnKjOQODAC0rrYd2zWKPX+M7DfDJpZXbyxl2\/JFsnMTx0rVnqqsec5t01aurHr\/t07k3K3sGhjhL57a4P0IuRIGz5jMw6Gdxe5\/OwMSmaD\/+MKh7\/frKZB711snZyERgZA6X8U77VUk5e4fMqr217ImlXu2vwYjebRLI5Yc8yRBCVz3J2JhM3rC2T+fASkbAIH7EWVpUREXAkOgm0Z0q0R0ZASNziNfcOC8Hf60A0bm3rrMuCscHFcRseHeRWIMM\/2TyOPVrEGnE2eC2DHQeiJTdQuKt9\/OKC\/KzMQTCMnZUzgU6ORsZAZNMBIarQzG0yRueisfSwhqhSCZsmAoc\/UaJohSBIc5SBCZKERhx\/JOf1fPZjbDOIzr9LAkY4dtJugihqx6\/EzwJGHURmWS3kEjAkIAhARMbj7KLI\/Kz9kIqSGxIwJjmkWRDm7w4Cse7h7rx0ws6XwcuU59MHmy5ztAmcVZ\/iJzC8ckJcOIscRaFjfnyO650+tnIRGC8nIFxm9isBIxTHp1pOuvCTreq70xNveN3o1Tb6PadKhkbOaGtBp3OgRXkZKC6H0QnYIWb27jy6+Cc+sitLqu9epk8bm0Oi42NhbNOeKYiTSVnZezXzVkZG2Xz6ORsZASM2xkY\/OL0uYvXYMmY24xXVuNrrC9eT4NZebEvT5v3Zj8svwRXWrQ1XoGN6nPuhgPwRd019rVqMQ3z41evj5w+z9LMdTmlcTvw6ScrG\/\/yN4A+HWOPdfOQqZ0dQdmIn1X481+vwfsFDd+V4mFGtD8jMxNqa2qU2bij7Bx8t0ua5aPs2HZ+if3pZuMPvpsLXVrUGHnNYVKdoU2vW0jIwWMXWrIP4OGFT8eJr13XwdkF92aBnR1O4wr5bPeV99KqdPjNwTOsTeYvwItpmD717p6s39AGu3FlZ+M\/bsyEQyfPWfISPyPi9CX6sNgoPsHRGDmLfSyDNeYz85JzhQ9yu60QJz7LcNZpy8XJRis7dHNWFke7seO2\/aTTzzYKAXP2vzNh658qEiZDLjq6t77CRIc4GSCJ8dFqcYIUJ2wUMVbfUeLCAh2qeW\/WLk0UP+Z8aMe6P9fD+H9KYzZyAcPfoWC2IygbURzcevM3gWMl2sFt\/Nnwhi\/LOtnhxcb7Xt0Dj+G7ebqkJeCIjxCv+mMtvDGqU1yak41oL77zx2lS1TmwvAgYv6I7KM5aiUUUraKwNjutHxfXQuY3rgHyWYXoDkLQRkV0f3L8jK0QjDpn+eSpeqEoI7plF4pYl4zotlroui3CVHPWbcEdxMJAJ2cjI2CctpDcvnnE08UtJKc8mIYX\/46SGG5Uncbt4P8Nsi5sk539bnacqa63\/OyCajy4jVb1OdmI+fCV5Xf0yAT8+jjHEd+4fLr6CutLnaFNty0kr19KD5qzaMfSjyos+xbxnD00N+GL7dx2jqvIWRFv7BOvaXZ1YRl2NiIfbu6QbnxVXqzLqV1hsbEpcNaN55xDfriC\/ccv81jfdaLG8AFimZwPnM9eeSlTF9ZrxzHdnHXjuhWOaOOuk7Xspa9W\/WKHfSr8bJMQMPxbHH4mA6cJW2Uan5S5jWYBo7IuL+LACiv8DS+rD1+i\/SpsFPvGrj4zVmIetwnrju4ZMP77\/UHXK67dBEwQoltFP3gR0G6cDUL8+xW0brwMg43I2ac2lNtOFE2BszL9ICO63RZodnbI1BWVhWIQCwPdnI2MgHE6A4NbB\/ycC4a6cWLbP6Mra1uL9p1hylsH2HaEuIWEoTxxO8ec9sXFa\/DbmYPYChLTsDwse8LqfYBpfHuD1+WUhnagjce\/aNgW4fnQjudH97O10WwHf0U35lNpI25P\/O4MGDaKeEx770t47sHexvkSN6y82IjbJ7g\/bIWjeXuJY+XHRsRc3PZAvP6w5F+1Chg\/nDWfg0AMdHAW68HLvF3nh7Mi1mZeuqXxdKdxZWcj1oXbvFZYIWfT09LAatszCjYiLo2Rs7J8ED\/3YuaKOHbQz+J5MvNWvZkrbnZw32HFFZGzrJ8sjhPY2ZhqznrBkduIbRC3iDEvptV81dJyXKWCdMDDPAAAIABJREFUs5ERMFaG0m+EgF8EdEZg\/NpG9xMCVggQZ4kXUUNAF2cjIWC8dF7he6cgO+MbUDgkm91eWXsVNn1ynv0\/\/81cjlOeytqvbMtTnbanog4G5ba1tF11XTFs7Ns2KLcd7Km4yO7BK69rW3i4\/zchCKywX7Auv\/VZ2Vj6WR0UjepmSRXkwk+2nILNU3p7oZKWe3i7rbhZtLMSsI34L2jOmrkg1od2WHGzaOdZOFNz1RJPM0\/E8pzS7OrCMexk46BlB2FQt3ZG3\/Nx75QnLDZu\/qQK9sz6TqPmrBPPZfsh\/5flUDiki+X4QL5sPlCVwE20A3lrlc\/JDpm62Jxjmo\/ETtbJWZm2IVaIodUYd8IjLH42dF+j9jqj2E2GTvmd8uhM01kX4hEWrGT7xpxPVhB45VYQ98mIbre+k+GRmx1WZaJgtFswyIpulYI2CqI7J\/NGNg79itgguOi1TDeuWJWjevEjI7pR2MgsFGXrCtNC0W48ql4YOC28vPJLxX2RFTAqGk9lRBcBGeea6tbKCMkgbJaxQ0YsuQkwVYKWlxN2G5sKZ1X3gxtudvXJ2CFblyzXZWx0q0tVpNttYWC32xGEz7IrM7QCZu\/evTB58uQ4u9esWQMDBw60xefkyZMwffp0OHfuHMtbWFjI7q2pqYGZM2cCltm5c2eYNWsWLFu2jN3Xr18\/eO211yAzM9Mo98qVK\/DCCy9AdXU1vPjiiyxty5Yt8Nxzz7F7hg0bBg8++CAsXLgQXn\/9dTh48KCR1rt3bOtizJgxMHr06ARbuY1ZWVms7GeeeQYOHz7M7CgqKmL\/tm3bFpcP68Gy\/GKiEo9XXnkFRowYAe+++y6z2YwJ9svp06dh1apV2vB464Pfw6vznvLFkSAHl9\/+QVtU9hFxNn4ME2fd2R4Vzk6fvwy2l+yDXet\/yhrV6ZaekNmqZSB+1q9fUTmGg+ZsWttM6J9fCL9duxQunj0OXbr1hrWrfu447\/jFw5116u4IpYDBSRzFiihY+EAThYkIA4qOZ599FvLy8thkj4KjtLQUFi1aBLW1tTBu3Dg28eKbQrEMLgrWrVvHJt0LFy7EoXr9+nX296233gqPP\/44LFmyhImV7t27M+Gxfft2WLFiBaSlpcH8+fOZ4ME0LA\/\/jRw5Enbt2hUnjriN+fn5rGzRDrR37dq18Mgjj8DYsWOZWCkpKWF1Y7sqKyuZ0PGKiWo80F7EJCMjA9555x04f\/48zJs3j2HSoUMH5kjwNP+MGTMYFmHDQ92QsS6JOEucVT2GibPkZ8nPOo+C0AkYjJbwyIQYFcFm4KNdjz76aILYME+u\/CNVXAQ89NBDMHXqVFi9ejV06tQJVq5cCZMmTWJRBKxv7ty58Omnn8ITTzzBxI85ArNz506GIo+ooB28vK+\/\/hp27NgBjz32GLtHTMPJHsvkwqNnz55G2\/BesZ1iPrQfbUA1XlBQAM2bNzfqEz\/AhbZPmzYNystjL+MTL1FsqMDDvDIQMUE7FixYwKqfMGEC7Nu3LyV4YP1O\/AlqQiDOrmbilTjb4Be84kGcbYh0k5+N+Xwx0k1+NmIChkcOMOIgbu3gJIFi4KabbmKOcs6cOSziwS\/ReYq\/o4jZsGEDHDp0iEVghgwZEocIpm\/evJlFaoqLi1lEAVfTmOfIkSNM2GBdly9fNiZlDBni\/fX19ex+JBlGe+6\/\/35mI74FFqMwogDjK\/T27dvD888\/nxAJwvrefPNNmDJlihGB4RGOq1evsi2wdu3asfp4uVzAoErHaFCfPn0Cw+P9999nA2vp0qVMkJkxQVG4fv162LhxI+zevTtleCD+2dnZrH\/S09OD0ixx5RJnibO4uJEdw8RZ8rN83iE\/689lhy4Cw80Xz5zw38SzIFx0iJMU34tcvHhx3FkZFCkY7bjhhhvgpZdeMoQP\/s5FAhc9fDLCcy18ewQdTFVVFfsbL75lhFsnPMKCv+P5GpzgxXMi5u7gW2F45sV8Bkb8G\/PhdlHHjh2N+sSzNiImOTk5hggLGg\/EiYtJMyazZ882BE4q8bA6e+RvWMjdTZwlzsqOYeJszJeSn43NO+Rnvfng0AoYb+bTXYQAIUAIEAKEACHQFBEgAdMUe70RtNnp5L\/dGapG0GxqQoQRIM5GuPOaqOlh5ywJmAgTk2\/l4JkUq8vqEfEIN9cw3ekJKzyXRAImvL1MnE18SpI4G16+omXE2fBylgRMuMcOs85pAOHTRm3btmUHnvv372+0xi1Pq1at2GFo84Xl+U3jj5w3a9ZMSXlYiJsd\/HFu8xNW+HQZCZjUk9qNf8RZMF6VQJxNPV\/Jz1rPBWH3syRgwjF2XK3gj3s\/\/fTTcU9fYUbxSSrxEK9THtVpqsvjDgUfcTe32emJMzzcjRd\/2Z4rsHRDYAgQZ709JUmcDYyCvgsmzkaLsyRgfFM8dRnshIqTRU55VKepLs9JnDk9cbZ8+fKEtyunrteads3E2YZH+Ymz0RgLxNnocJYETDTGFFlJCBAChAAhQAgQAgICJGCIDoQAIUAIEAKEACEQOQRIwESuy8hgQoAQIAQIAUKAECABQxwgBAgBQoAQIAQIgcghQAImcl1GBhMChAAhQAgQAoQACRjiACFACBAChAAhQAhEDgESMJHrMjKYECAECAFCgBAgBFIiYI4ePQr45d5evXrByJEjoWXLltQThAAhQAgQAoQAIUAIeEZAu4Cpq6uDjRs3wrhx42D37t3QunVrGDhwoGeD6UZCgBAgBAgBdwRooeiOEd0RbQS0Cxh8G2VpaSkUFBQA\/v+JEydg6NChvlHc9Ml5KK2oi+VrBpB\/WxYMym3nuxzKQAh4QUDFZECc9YI03aMCAVULReKsit6gMoJCQKmAwW2hzz\/\/HAoLC5m9\/KvB27ZtA\/5l5Orqajh27BgMHz5cWsDgoKqsvQqFQ7INXIp2VkJl7VdQNKpbUFhRuU0UARWTAXG2iZInRc1WsVAkzqao86hazwgoEzBFRUWwZs0amDx5siFgREGD6bfccgvcddddsHXrVhg\/fjwcOnQIqqqq4J577vFs8J6Ki4D\/RPHCM6OIwSgMRWI8w0k3ekAg2cmAOOsBZLrFMwI6ForEWc\/dQTemEAElAgYHVE5ODpSUlLCmYASGR1\/y8vJg9OjRxheTFy5cCMXFxUy4XLp0CSZOnAhZWVmeISjccgqKRttHWQrfO0VRGM9o0o06JgPiLPFMFQK6ForEWVU9RuUEiYASAWNEQIqKEgRMfn4+O6SLX\/j08pXgvn37BtleKruJIHDq1Cm4fPmyY2tVTQbE2SZCqoCb6cZZlQtF4mzAndlEinfjbNAwhFLAfPzxx3HtrqyshOzs2HmXBTuqYMG9DRGbAb\/4DPbP6MrSvrh4jaW\/MaoTiHnMIIYljexIpLcdJn6xwm3JsrIy2\/GjejLAgSxeFRUVkJuby356cn05rBjb20juUPg7qC66m\/19uroenlpfDr9+8jYQ85gND0sa2ZFIKTtM\/GLVrVs3R86qXiiSnwWt84RfH8b7W2c+v3W5+dnICxjzFtKiRYsgPT3dtl24MhAHFk4C6Aj4dfa\/M2Hrnypg2vcy2ASBkwEXMJuON4Pura\/AgC5p7HZzXvFvpzSVeXXVY7ZZV7266vHbPq8DC6MweInbnjJRw7ffftvgKIptdAT8+muzLPjdkf+CMb2bMyGOovuDh2LvPio+2wZuvrEO+mY1Y3+b84p\/O6WpzKurHrPNuurVVY\/f9k2YMCFlAob8bMMc48en+bnXrw9LxXzl10avfjYoIRNYBAYNtjrEi+dhnC43AYOdOmH1PvhW25awZMxtTMB8OCEb\/m\/5JcjIzISHv33dKJ4ETPCDUtcADmpgqRIwTqKbOOuNh7q4pKsenZxVvVAkzhJn3RbymN6oBYz4GPWwYcPALfqCgKCAcQrHc3VScqIWdp2shZc\/qoCxt3eCwT0y2X\/55Td8m4p8YbExLHZgH4QhHC8zGRBn7fvOqV9l04izsahhMgtF4ixx1hxM8DuuvG57RiICo8JIrwKG1yWeJxDr99sRJGBi5zWSJXRQOPrtT68DS4zABD0ZEGfVccwvH4LiZapFN9Yf5EKROEucdRo7Xv2sCm1gVYbSLSQVRnrZQuJnYsxnYPyEhv3ci+2S3Y\/UVU8yNiaTN6zt8xraNAsY2cnAbQuJOBvzDrSta7814ZWz5GcTuRRWP6SK72Ftn07ORkbAeAlt0sqAVgZhWRlQ1DDWEzqjIjrrkm2bXxt1rmaJs8RZ\/pSkKAzCzFkSMP+DgN9OCir8THYkUlL3GRhVq1kS3SRgzFzyO75JwAQjKvz2Q2P396nY9lThZ0nAkIAx3k2SjOoOaoBHdWDRajaYiUclH4iz8e6fOEucpQhMALKKzsA0vPcmrPuesueBkC668urcmyXOEme5K0zmzANxVt4\/RMFXJuP\/wto+nZyNTAQGD0TyDnP6LzaIv8hOVJPigUmrVRwPKVr91+yIrOq3Ws2Z6xTrdarPrn2c7E72u9khU29YcbTCw6l9OgcWX82KoWonPvAn5\/z0j04+cM6JdVr9Zg7N6xxXfvlgFs5hHFe6t5DIz8YOVHOfZ8dn8rOxt4pbjW+dfjYyAobOE9B5AjNZ\/e5j654MiLPEWeJsPAJ+x2xQ23xkR+LUH8WzhpERMPRIasPKgK\/kzeFHt791hRx11eM3\/KpzZRDkFlJpVTr85uAZNnbbtGkDd+eAp09l+MVL19ZeMlssqvISZ2MvDCU\/S3422XlEp59t8gKGJgP1r8emySC4yeDD8ktwpUVb4\/MYiPXcDQfgi7pr7IOmfrAP6t5khFKq8vrBQqeNOieDIAUM+Vnys1ZiI4jfQvkiuyDC8ev3fcm+\/DtnaFdjL+\/ljz6D09VX4r4UrDPcqLMuqz3\/oEK2TnXJpvnFKupbSPxTGchXM2bI2zu6Z8DgHhkJabxPke8lJ2rg0qVLLHIz9vbOxv3J9LvffkimrjBwNhkc\/WIVdc5if5GfVfc5lKDGjspxpZOzkYnAJCNgrBwONhy\/m+Q0GZypqbd1+Mk4Mb9kka1LJp9MHpFEfh20TH1OeZzSdA6sZB9JtWrH+n3nHIX1k+vL2fe\/rERKMpOIHaay\/SCbz855y5Ynyz2nRY9qrKLEWav+IdEdQ0XV+RISMO4xm1BGYGT3ZjF0eejkOZj2vdjKlIfcV+2thf0zuto+wjvtvS\/h3r6d475kvel4Mzhy+jw7b3Du4jX25Wt+HoWnYRgfw\/zHLrRkK128RvRuAz8c3Jvda5f28+LDsP9sPbtfPNfgtmUghmYx79S7e0KXFjWsHicbORbieRo7rHB74u2pt4OdjViWU4hY1kYrrP5xY6Zlf6KN\/\/KdHNs07JeohONVc3ZH2Tn4bpc0NgbM51pw+2lAdhp8cfGaJWexz808euOPtfCXvwH06QgJ4wr7AceH3TaXOQ2543Vc+eWsEx+cbESuWHH2yy\/Pwf7KelscuWsV\/YIbVmHlrJMP4\/3An\/bkf1v5ByzHvLUpcnDBjio2bvmZLj8+zIuNqT4v6GajmdN+tjL93OtUj4yNTnOBTj8bmQgMFzBuYGOn8oGFguDk39PjRAjvdBQoj30vwxAWYgdjPhxYZQvuNASKKH6szhtwYWTnxC5eT4Ne7a\/5mgz+dLaeTTyiQ+T2o2PEKyMzM0FkyUxY32rbMqEssS6sp7amRtmE5WajFY7YJ3\/+6zV4v6DhC+OijYjxrLw0g9NiGk7Ss8f\/AMrKytwlvII7xPMEqjiL7RfF8IBffAbVRXczjiJnV\/2xFn47c1ACZ+97dQ\/jOk7aVgdzcSyYeWaeeMV8WNe6P9fD+H9KSzg8jPm44BUfG8f\/t0rzOq6sRJYbZ5344GSj3bjCr9zjosdqMvhVSTnD\/41RneIwdsMqjJy1Wvxg34m2in6WL2BULhRl\/INoo9+FouwiDPOZ63JaRHIbcS6wWrAirqoXiqptdFro6l4oRkbA4BaS1zA4f6cGhtRXjO1ttFEM42EaXjxd\/IL1AysOsLRfP3kb+69TPjENw6VLP6qwzGc+W8Pz8Tyzh+ay8whieWgH\/91sh1NdTvmwHEy\/o0dm3Nkf\/C2nQ5qBR7Ltwnrs2sxt8Ns27LObO6QbZz289gvWh3m3zdMrYFRzFvtc3EYyc9YOTzPXzX371IZyODh\/UALX7XjExxX\/7\/9r71pDrEiu8NFdk1HxPRE1I8aV+aFiWCZEFEwQySRM0Kjxgc81anwSSRxHxUd86xifG4gvNOoa8RF1B7NZJYOILL7Q4IqKikaN7mz8oY6PDTpK0PDVcO7U7dtdfbtu950ePQ3LOre6qk59\/dWpr09VdaX7HEx2mNJsOGtbl22\/Qn1fVVal9H1brGqLs7vKL\/hOrfM6K+afaZoI7dfXWumctcXaz8+yw3eubTz1r8cpfo99FKZevaYGbdL0upy+W\/e13HeisgPrOb3st7HRlIfX4X3004KsvSjWGQFj07FMAgadQHfcQQYDN0fFDseUpgsmJq6tg\/NymGyHW134zdlu3Q6TOAjaLj87\/AYYm8HAy0YeAK+sHZS1joUITNicRTs+XHZGDQhwzuAshAeEPYQNixA30a1jk6lYZ87CEUPkOwWMbf\/wyhcFZ21tDPrSY4tVbXG2aPnnvuusnC99fn42bNFt8n1e4sb0omh60bJN8xL\/LPZ4o4jed2zrss1nY6Pfy3G2XxTrjID50Zz9ai4Vl9vc35it5xNzrTyFhJA71m5wKFsPuWPOv98nFTT7Zx3VFAzSMABsKr9Gn13\/L302Ji9RD9J4vQzqQV4OFetf\/UWoHpczjIzf1pytohv\/SQ0xs43Ii3x6+5CGunhaQLcfYX9cblMGNjayHYyjs66GOTn0x59XryMy2eiXxumZ2MhY6TbC\/sev3k\/Y6HzWSL+1YWRWBUwUnEU7FvbvTJ\/\/8w5hHdfwH7ZRU5Nq2lPbRq1zFiFphLm9OOs2veTGWb1fLR3Ujcbv+FKV6XwOXpx143Mm\/crEWRMfbPsVpomQl32Ds+\/rOKaLVV3jLPc9fQrJ6WedPhp+dnTP9r5+VvctNv4hHT\/L\/s3pp\/C3PlbwGIO2OdN0H4Z\/u63xQT1uaX7+Husf3abB3eyIm41oL8albPrZOiNg3AyV3wSBoAhkcw1MUNvkfkHADQHhrPCiriGQLc7WCQEDI4s\/vU3rfvmB63OsePKSph+6TQfG16x34RuRL6\/5t6i4T576Cff+9cID9W\/+za1QU76KJ69cy3T+rpfrLI\/Tztx5RuuOf03Ffb5LPTs2TTJl3fEKOnPnG+rZsUmK\/aa6TPkOXHhIZ0o+dMURNp65\/SyRzljZtIufmY69Xqlt23quuUg9P2ia4IL+PL2ei9+zjso5xImz4Dq4hv+AE64e32tC7Vt8W\/3m1hdMz8iLl9V9zL1\/2KbZcNa2Lr9+xX3UiePQgu8oH+Xma2yxMvmnqDjLHPHiA9rv9FM2eWC\/Xz4b3zfkz9dcfSmeDXys2ziBtK8ev\/QcP2zSMMYc+PJhSpnMETcco7DDVKaNjV558DzhL9zaFRVXvcqN3TbqdMhuAs7NccPh+F2mfF5pNsLH5OC8Bh6Tw+QB22bAYkz0Qc6vLtsBy7ZteNZuA4jOE6f9fs86inQ\/B51tznq10U\/ke3Hdpn\/4PaMwOWtTlx\/X\/URF2FhFwUu\/Mv34EPSFz\/ZF0cY\/eAlQ04ui18ul30uYKZ+X+MfLqkkQeb3wxcVG00tNbb0oOvkYSwHjRqZ0oyl+HTbsdBvHHoXIsh2wvPLZtMtvEDFhb4tJ2M8zk\/JsBoNM6rPNG3ess41j3PGwfc7p5rNpv02eKPxD2JGwsF\/Qzv77G89ZAdu6bPN5iUSTjSZhmS6\/orwvtgLm3LlzNG7cuKS2fzRnPc0aWeiJx61bt2jy5Ml0\/\/59lbe4uFjd+\/jxY5oxYwahzLZt21JJSQmtWbNG3detWzfauHEjtWjRIlHuixcvaNmyZVRZWUkrVqxQaYcOHaKFCxeqe4qKiqh\/\/\/60ePFi2rx5M128eDGR1rlz9dTWsGHDaNCgQSm2so25ubmq7Llz59Lly5eVHevWrVP\/HT16NCkf6kFZbphs376dunfv7oqJjke3nwyltj8epe77fqvX9MUnKwLhsXr1aurXrx\/t2bNH2ezEBDbcu3ePtmzZEls8ouxIKFs4W4OwcDacPlwbnDX5FNgTNz\/brqAwJUr7g2bPaOyvJ9J7jZpT\/i9+S1c+\/ZiefX0j4WfnLv0Dnf\/imKuf3XH4BK2d95ukNMbETbihLr9xp3GL1tSp8Fd04x87qOrpg8S4c\/3Jewnb\/\/eqip6e\/gv9buxQ5WcLR06jK5X16cyxv9O1v\/1J2cN+dsrv11L5yfN0au\/H6vc2HfKpRaP31bjjhQds9Bp3THgEHXei5qxefiwFDAZxEEbvSAyiLkz0hkB0LFiwgHr06KEGewiOs2fP0pIlS+jJkyc0evRoNfBi5TfKYAe7a9cuNeg+ffo0Cfc3b96ov7t06UJTpkyh0tJSJVY6deqkhEd5eTlt2LCBcnJyaP78+UrwIA3l4b8BAwbQqVOnksQR2zhkyBBVtm4H7N25cyeNGDGChg8frgbDkydPqrrRroqKCiV00sUkbDxgLzBp3rw57d69mx48eEDz5s1TmLRs2VJ1HKyUnzZtmsIibnhE3amEs8LZsPuwcFb8rPhZcy+InYBBtIQjE3pUBM3AFulRo0aliA3n4IqBlN+IIQIGDhxIEyZMoK1bt1KbNm1o06ZNNHbsWBVFQH1z5syhmzdv0tSpU5X4cUZgjh8\/rsrjiArs4PJev35Nx44do0mTJql79DQM9iiThUd+fn6ibbhXb6eeD\/bDBkQ9Ro4cSfXr10\/Ux21Dftg+ceJEunat+kN9+qWLjTDwcEZgdExgx6JFi1T1Y8aMofPnz9cKHoyJF3+iGhCEs1uVeBXO1viFdPEQztZEusXPVvt8PdItfraOCRiOHCDioE\/tYJCAGGjdurVylLNnz1YRD75056n\/jkjGvn376NKlSyoC06dPnyREkH7gwAEVqTly5IiKKOBtGnmuXr2qhA3qev78eWJQRggV91dVVan7QTJEe\/r27atsxKFsiMLoAozf0Js1a0ZLly5NiQShvm3bttH48eMTERiOcLx8+VJNgTVt2lTVx+WygIFKRzSoa9eukeFRVlamOtaqVauUIHNiAlG4d+9e2r9\/P50+fbrW8AD+eXl56vk0bNgwKs2SVK5wVjiLlxvbPiycFT\/L44742WAuO3YRGDZfX3PCv+nz6iw69EGK52aXL1+etC4EIgVv5Q0aNKCVK1cmhA9+Z5HAoocHI6xr4ekROJiHDx+qv3HxlBGmTjjCgt+xvgYDvL5OxPk4eCoMa16ca2D0v5EP00WtWrVK1KevtdExad++fUKERY0HcGIx6cRk1qxZCYFTm3i4rT0K1i3s7hbOCmdt+7BwttqXip+tHnfEz6bng2MrYNIzX+4SBAQBQUAQEAQEgXcRAREw7+JTfwvabNoJke01MG8BnNKELCAgnM0CyFJFqAjEnbMiYEJ93NktjKdysCbF7XLbIp5dC6OpzbTDCuuSRMBEg3sYpQpnU3dJCmfDYFZ0ZQhn48tZETDR8T60kk0dCLuNmjRpohY8FxQUJOr0y9OoUSO1GNp5obygabzlvF69eqGUh0L87ODt3M4dVthdJgImNOpZF+THP+Fs9XeDsEtSOGtNs1AzCmeT4TTtZI0LZ0XAhNoFoiuMt3vPnDkzafcVatR3UumLeE15wk4Luzy0y6tM044zLO7GxR\/bi+6JSMl+CAhn09slKZz1Y1L20oWzdYuzImCy1zcyrslLqJgKNuUJOy3s8kzizLTjbP369SlfV84YfCnACgHhbM1WfuGsFYWynkk4W3c4KwIm691DKhQEBAFBQBAQBASBTBEQAeOBoH4UAU\/LmMKLmT4Im\/z60QRe5yHZlOvM41w0G0aZUkb4CAhnazAVzobPryhKFM4KZzPhlQiYOixgMnnwQfLKYBAErdq7ty4MBtlCRzibLaQzq0c4KwImEwaJgMlAwDhP6eSDFnFswPXr19Un9fH14Lt371Ljxo3pxIkT6jP8OMSxQ4cO6gRrfL2Xv7zITpdPo+byUA\/+jQu7FjhPu3bt1EGPOBwSERh9Fb3XFuqkE6q1k7idNuPLoHz8Qa9evVTdhYWF6jyoIGVkQk7JGwyBdAYD4WzyCfRuvA+GutydCQLCWfGzmfBHBIylgEE2\/RRqvSPiXKD79+8nzuOBk8QBh9jqXFlZqY5eLyoqouLiYiUScPG\/OR+EDh9z8OjRIyV6IGIgTCBaIGL4pGoIGP6dxQzsgXBCuXyxwBk8eLASIqib63PajPwHDx5Mshl24ywpfPJ7+vTpSjSZysiEmJI3OAJ+g4FwVjgbnFXR5hDOip\/NhGEiYAwCBhES56VHTPQ0feU6xACLEvxfFynOdTQsNHQxAmGgh8Bx1pG+s8YtD85MwkmmpaWlSYdI6jYicqLfg79ZhB0+fDjFZkSJ+IwWtAF\/O23Ry9yxY0dSGZkQU\/IGR8DtLCaUIpxN7j\/C2eDciiqHcLbar4qftWOYCBjLCAwO2+IpFi4CURWcggwB40ZKkNRLwODDQHpkg4UPiwb98Eo3AYP7\/bYPQ2Tp9+i2QMCwzW7rB3QBg2iQfvEAqZdhR0fJlQkCfm+zwtkadIWzmTAtvLzC2QXUo0f1l26dPl\/8rD\/PRMBYChhM6+hiwBmBCSpg\/CIwfgImjAiMl80mMaXDxyKntk7W9af7232H32AgnF1C+ocedV4LZ2unbwhnvSMwus8XP+vOTxEwIQiYnJwctS4Fl20Exm8NjJ+Ayc\/PT4rgeDkGRHm81sDoAkaP1vC6Hbc1MPpaGUwh6WXUjkt8d2sNMhgIZ6vXdwlna7e\/CGdrouLiZ4NzUQSMpYDhHUBIBMTEAAABAklEQVTYMYRwdElJCZWVlak1KE6nqEcmvKaQIGBMu5D8BAzWzXjtDtKbaNpB5BQfPEWGBcK9e\/em3NzclF1I+voKicAE74Bh5vAbDISzk9WideFsmKzLrCzhbPUaSd6gIX42GJ9EwATDS+4WBAQBQUAQEAQEgRggIAImBg9BTBAEBAFBQBAQBASBYAiIgAmGl9wtCAgCgoAgIAgIAjFAQARMDB6CmCAICAKCgCAgCAgCwRAQARMML7lbEBAEBAFBQBAQBGKAgAiYGDwEMUEQEAQEAUFAEBAEgiEgAiYYXnK3ICAICAKCgCAgCMQAAREwMXgIYoIgIAgIAoKAICAIBEPg\/+rvydjAIrvdAAAAAElFTkSuQmCC","height":337.44619799139173,"width":560}}
%---
%[output:0f0272c1]
%   data: {"dataType":"error","outputData":{"errorType":"runtime","text":"Error using <a href=\"matlab:matlab.lang.internal.introspective.errorDocCallback('PhasorArray\/mtimes', 'I:\\Maxime\\OneDrive\\Documents\\GitHub\\phasorArray_Toolbox\\Fonctions\\@PhasorArray\\PhasorArray.m', 710)\" style=\"font-weight:bold\"> * <\/a> (<a href=\"matlab: opentoline('I:\\Maxime\\OneDrive\\Documents\\GitHub\\phasorArray_Toolbox\\Fonctions\\@PhasorArray\\PhasorArray.m',710,0)\">line 710<\/a>)\nThe two arrays must have compatible sizes. Left array o1 is 6 x 3 x 1 and right array is  o2 6 x 3 x 41"}}
%---
