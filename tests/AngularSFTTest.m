classdef AngularSFTTest < matlab.unittest.TestCase
    methods (TestClassSetup)
        function pathSetup(tc)
            root=fileparts(fileparts(mfilename('fullpath')));
            tc.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(root,'Fonctions'),'IncludingSubfolders',true));
        end
    end
    methods (Test)
        function testFirstCompleteWindowAndPartialWindow(tc)
            th=.73+(0:.1:7); t=th-th(1);
            [p,~,~,~,s]=angularsft(th,t,1,ones(size(t)),0,{},false(1,5));
            expected=min((th-th(1))/(2*pi),1);
            tc.verifyEqual(p{1},expected,'AbsTol',2e-14);
            first=find(th>=th(1)+2*pi,1);
            tc.verifyEqual(s.meta.istart,first);
            [k,f]=find2piAntecedant(th);
            tc.verifyEqual(k(first),1);
            tc.verifyGreaterThan(f(first),0);
        end
        function testVariableSpeedAbsolutePhaseAndMethods(tc)
            t=linspace(0,3,12001); w=2*pi*(4+2*t); th=.73+cumtrapz(t,w);
            x=2+cos(th+.3)+.4*cos(3*th-.2);
            for method={'angle','mixed'}
                p=angularsft(th,t,w,x,[0 1 3],{},false(1,5),method=method{1});
                valid=th>=th(1)+2*pi;
                expected=[2;.5*exp(.3i);.2*exp(-.2i)];
                tc.verifyLessThan(max(abs(p{1}(:,valid)-expected),[],'all'),1e-4);
            end
        end
        function testPlateauShapesAndInputErrors(tc)
            t=linspace(0,2,2001); th=.4+max(t-.2,0)*2*pi;
            x=[ones(size(t));cos(th)];
            [p,~,~,~,s]=angularsft(th,t,[],x,{0,1},{'DC','AC'},false(1,5));
            tc.verifyEqual(numel(p),2);
            tc.verifyTrue(all(isfinite(p{1})));
            tc.verifyEqual(p{1}(s(1).meta.istart:end),ones(1,numel(t)-s(1).meta.istart+1),'AbsTol',1e-12);
            tc.verifyError(@()find2piAntecedant([0 1 .5]),'find2piAntecedant:NotIncreasing');
            tc.verifyError(@()angularsft(0,0,1,1,0,{},false(1,5)),'angularsft:InvalidTime');
            tc.verifyError(@()angularsft(th,t,[],x,{0,1,2},{},false(1,5)),'angularsft:SignalCount');
            [k,f,u]=find2piAntecedant(.4,2,true);
            tc.verifyTrue(isnan(k)&&isnan(f)&&isnan(u));
        end
        function testPlotOptions(tc)
            t=linspace(0,2,1001); th=.7+2*pi*3*t;
            [~,~,~,~,s]=angularsft(th,t,6*pi,[cos(th);sin(th)],[0 1],{'cos','sin'},false(1,5));
            fig=figure('Visible','off'); cleanup=onCleanup(@()close(fig));
            for orientation={'hor','ver'}
                clf(fig); figure(fig);
                plotAngularSFT(s,[1 1 0 0 0],orientation=orientation{1},plotDebut=false,plotOmega=true,xAxes='phase',Hm2plot=1);
                tc.verifyEqual(numel(findall(fig,'Type','axes')),6);
                tc.verifyNotEmpty(findall(fig,'Type','line'));
            end
        end
    end
end
