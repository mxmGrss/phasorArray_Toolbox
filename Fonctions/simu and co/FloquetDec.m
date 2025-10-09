function [Wf,Lambda,Phi0T,N] = FloquetDec(Aph,T,varg)
%Output Wf the phasorArray of W(t),

arguments
    Aph
    T=1
    varg.TransSolver {mustBeMember(varg.TransSolver,{'adaptative','forward-euler','RK4'})} ='RK4'
    varg.FixedStepTransPow=18
    varg.InitProbSolver {mustBeMember(varg.InitProbSolver,{'adaptative','forward-euler','RK4'})} ='RK4'
    varg.FixedStepInitProbPow=18
    varg.plot=false
    varg.nT=1
    varg.holdplot=false
    varg.modulo_eig = 0
    varg.precalc_Phi0T=[];
    varg.auto_adjust_precision=true
    varg.jordan = 'false'
    varg.tolSwitch2Jordan = 1e-3;
end

if isa(Aph,'PhasorArray')
    Aph=Aph.Value;
end

Aph=PhasorArray(Aph);
nx=size(Aph,1);

nT=varg.nT;
if isempty(varg.precalc_Phi0T)
    disph('computing transition matrix...')
    Phi0T=TransitionMatrixOverT(Aph,T,"simutype",varg.TransSolver,FSprecpow=varg.FixedStepTransPow);
    disph('computing transition matrix... Done !')
else
    Phi0T=varg.precalc_Phi0T;
end
Sphi=sym(Phi0T);

disph('computing eigen values...')
% [JV,JD] = jordan(Sphi); %extract jordan normal form of PhiOT (JD) and the similarity transform JV
[dV,Q] = eig(Phi0T); %extract eigein vector (dV) and diag matrix of eigenvalues (Q) of PhiOT the transition matrix

nf = get(gcf,'Number');
i=0;
for mui=1:size(Q,1)
    disph("simu ",mui," over ",size(Q,1))
    mu=vpa(Q(mui,mui));
    v0=vpa(dV(:,mui));

    lambda=double(log(vpa(mu))+varg.modulo_eig*2*pi*1i)/T;
    Al=Aph-lambda*eye(nx);

    Al.value;
    %     Al=Aph;
    %     Al(:,:,(end+1)/2)=Al(:,:,(end+1)/2)-eye(size(Al,1))*((log(vpa(mu))+0*2*pi*1i)/T);

    [y_a,t_a]=hmq_sim(Al,nT*T,v0,T,plot=false,solver=varg.InitProbSolver,FSprecpow=varg.FixedStepInitProbPow);
    n=length(t_a)-1; %on retire l'element final pour  avoir un vecteur de 0 à T-Ts
    N.Ny1(:,mui)=(y_a(:,end)-v0); %erreur, normalement à la fin on revient au debut puisque eigen vector
    N.Ny(mui)=norm(y_a(:,end)-v0)/norm(v0);

    Y_a(:,mui,:)=y_a;
    T_a(:,mui)=t_a;
    ty_a=y_a(:,1:end-1);
    W(:,mui,:)=ty_a;
end
fshiftn=((-n/2):(n/2-1))/nT;
fshift=fshiftn/T;

Wf=fftshift(fft((W),[],3),3)/n;

[~,I]=find(fshift==0);
if varg.plot
end


doubleDiag=not(eye(size(Q,1))+diag(ones(size(Q,1)-1,1),1));

Lambda=double(logm(Q))/T;
Lambda((doubleDiag))=0;

[~,I1]=find(mod(fshiftn,1)==0);
Wf=Wf(:,:,I1(2:end));

plotres()
pause(0.1)

if nnz(double(N.Ny>1e-2))>0
    N.Ny
    warning('WARNING in FNC FloquetDec : periodicity deviation >1% ((xsim(T)-xsim(0))/xsim(0)), consider increasing Fs power')
    if varg.auto_adjust_precision
        warning(['Increased Fs power to ',num2str(varg.FixedStepInitProbPow+1)])
        varg.FixedStepInitProbPow=varg.FixedStepInitProbPow+1;
        C=[fieldnames(varg).'; struct2cell(varg).'];
        C=C(:).';
        [Wf,Lambda,Phi0T,N] = FloquetDec(Aph,T,C{:});
        return
    end
end


    function plotres()

        if varg.plot
            for muii=1:size(Lambda)
                i=i+1;
                ty_a=squeeze(W(:,muii,:));
                y_a=squeeze(Y_a(:,muii,:));
                t_a=squeeze(T_a(:,muii));
                figure(nf+i)
                clf
                subplot(221)
                plot(t_a,real(y_a))
                title(['real v_{',num2str(muii),'}(t)'])
                subplot(223)
                plot(t_a,imag(y_a))
                title(['imag v_{',num2str(muii),'}(t)'])
                subplot(2,2,[2 4])
                tt_a=t_a(:,1:end-1);
                stem(fshift,(abs(fftshift(fft((ty_a),[],2),2)))'/n)
                set(gca,'YScale','log')
                xlim([-25 25]/T)
                grid minor
                title(['stem FFT(v_{',num2str(muii),'})'])
            end

            i=i+1;
            figure(nf+i)
            PhasorArray2time(Wf(:,:,I+((-25*nT):nT:(25*nT))),T,0:T/200:nT*T,plot=true,explosed=true,hold=varg.holdplot);
            sgtitle('Matrice W(t)')
        end
    end


end