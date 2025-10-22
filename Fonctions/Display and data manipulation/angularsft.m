function [phasor_cell,theta,omega,IDX, phasorStruct] = angularsft(theta,time,omega,signals,harmonics,NameSignals,PlotTAPRI,optarg)
% ANGULARSFT Compute Sliding Fourier Transform for a non-uniformly sampled periodic signal
% where the fundamental frequency varies over time.
%
%   [phasor_cell,theta,omega,IDX, phasorStruct] = ANGULARSFT(theta, time, omega, signals, harmonics, NameSignals, PlotTAPRI, options)
%   
%   Inputs:
%       theta      - Angular signal (phase) driving the variable frequency of the signal.
%       time       - Time vector corresponding to the signals.
%       omega      - Vector representing the pulsation (frequency) at each time.
%       signals    - Cell array of size p containing the p signals to be analyzed with respect to theta.
%       harmonics  - Cell array of size p specifying the harmonics to analyze for each signal.
%       NameSignals- (Optional) Cell array of names for the signals, default is empty.
%       PlotTAPRI   - 5-element logical vector indicating which plots to display:
%                      1st: Time-domain signal, 2nd: Abs of phasor, 
%                      3rd: Phase of phasor, 4th: Real part of phasor, 5th: Imaginary part of phasor.
%                      Default is [true true false false false].
%       options    - Optional name-value arguments with the following fields:
%           .xAxes       - X-axis type: 'time', 'phase', or 'revolution' (default: 'time').
%           .method      - Analysis method: 'angle' or 'mixed' (default: 'mixed').
%           .orientation - Plot orientation: 'ver' or 'hor' (default: 'hor').
%           .plotDebut   - Logical, whether to plot the starting phase (default: true).
%           .plotOmega   - Logical, whether to plot omega (default: false).
%
%   Outputs:
%       phasor_cell  - Cell array of phasor representations for each signal.
%       theta         - Adjusted angular signal (phase).
%       omega         - The pulsation (frequency) vector.
%       IDX           - Indices corresponding to the angular components of the signal.
%       phasorStruct  - Structure array containing signal information:
%                       .name      - Name of the signal.
%                       .phasors   - The phasor representation.
%                       .signal    - The original signal.
%                       .time      - The time vector.
%                       .harmonics - The harmonics for the signal.
%                       .theta     - The angular signal (phase).
%                       .omega     - The pulsation vector.
%                       .IDX       - Indices corresponding to the angular components.
%
%   Formula used for each harmonic h_k:
%       h_k = integral from (t - Ti) to t exp(i k * theta(tau)) * omega(tau) * signal(tau) d(tau) / 2π
%       where Ti satisfies: theta(Ti) = theta(t) - 2π
%
%   Notes:
%       - If `theta` is empty, it is approximated using cumtrapz over `omega` and `time`.
%       - If `omega` is empty, it is computed using the difference between consecutive `theta` values.
%       - For more precise results, provide real `theta` and `omega` signals.
%
%   Methods:
%       * 'angle' method:
%           - Uses angular representation of the signal, i.e., computes `exp(-i * k * theta)` for each harmonic.
%           - Formula: `h_k = integral from (theta - 2π) to theta exp(-i k * theta) * signal d(theta) / 2π`
%       
%       * 'mixed' method:
%           - Incorporates both angular signal and pulsation, i.e., computes `exp(-i * k * theta) * omega` for each harmonic.
%           - Formula: `h_k = integral(exp(-i k * theta) * omega * signal) / 2π`
%
%   See also: cumtrapz, shift2pi, validate_input, plotAngularSFT
arguments
    theta
    time
    omega
    signals
    harmonics  = 0:5
    NameSignals=[]
    PlotTAPRI {mustBeNumericOrLogical}=[true true false false false]
    optarg.xAxes {mustBeMember(optarg.xAxes,{'time','phase','revolution'})} = 'time'
    optarg.method {mustBeMember(optarg.method,{'angle','mixed'})} = 'angle'
    optarg.orientation {mustBeMember(optarg.orientation,{'ver','hor'})} = 'hor'
    optarg.plotDebut logical = true
    optarg.plotOmega logical = false
    optarg.plotlang = 'fr'
end

[theta,time,omega,signals,harmonics,NameSignals,PlotTAPRI,~,~] = validate_input(theta,time,omega,signals,harmonics,NameSignals,PlotTAPRI,optarg);

%first we check that d theta/dt is not changing sign
% dtth = gradient(theta);
% dt   = gradient(time);
% assert(isempty(find(dtth./dt<0,1)) || isempty(find(dtth./dt>0, 1)))
%
% %First case, increasing theta
% assert(isempty(find(dtth./dt<0,1)))

% Preparing future interpolation

% Call the shift2pi function
[~, istart, ~, ~, ~, L, ~, IDX] = shift2pi(theta, time);

% Use the outputs from shift2pi
% vector of integer part of theta -2pi
% after the loop, L(kk) is the index of the last theta < theta(kk)-2pi
% ie theta(L(kk)) <= theta(kk)-2pi < theta(L(kk)+1)
% if theta(kk)-2pi < theta(1) then L(kk)=1

% vector of fractional part of theta -2pi
% ie F(kk) = (theta(kk)-2pi - theta(L(kk)))/(theta(L(kk)+1)-theta(L(kk)))
% after the loop, F(kk) is the fractional part of theta(kk)-2pi
% ie theta(kk)-2pi = theta(L(kk)) + F(kk)*(theta(L(kk)+1)-theta(L(kk)))

% First step is given theta, find index IDX such that
%           theta(i)-2pi = theta(IDX(i))

% index where theta starts being >2pi+theta(1)
II2pi = find(L == 1, 1, 'last');

% index iimove where theta starts being >=theta(1), ie theta(1:iimove) are almost all equal to theta(1)
II0 = find(theta == 0, 1, 'last');

% find first index istart such that each element of theta(istart:end) is unique
% and theta(istart:end) is increasing
% istart = find(diff(theta) == 0, 1, 'last') + 1;
% istart = iimove; % This is now handled by the shift2pi function

% check if every element of theta(istart:end) is unique
% This is now handled by the shift2pi function

% Interpolation
% IDX = zeros(size(theta));
% IDXi = round(interp1(theta(istart:end), istart:numel(time), theta(istart:end) - 2*pi));
% IDXi(isnan(IDXi)) = 1;
% IDX(istart:end) = IDXi;

% index st theta(IDX(ii)) =theta(ii)-2pi
% but IDX is not an int actually
% theta(IDX_int(ii)) <= theta(ii)-2pi < theta(IDX_int(ii)+1)

% IDX_raw = (interp1(theta, 1:numel(time), theta - 2*pi));
% IDX_raw(isnan(IDX_raw)) = 1;
% IDX_int = round(IDX_raw);
% IDX_frac = IDX_raw - IDX_int;

shift = [theta(IDX == 1) - theta(IDX(IDX == 1)); 2*pi * ones(1, nnz(IDX > 1))'];

phasor_cell=cell(numel(signals),1);
for ii = 1:numel(signals)
    try
        Signal_ii=double(signals{ii});
        if ~isrow(Signal_ii)
            Signal_ii=Signal_ii.';
        end
        Signal_ii(Signal_ii==-Inf)=0;
        Signal_ii(Signal_ii==Inf)=0;

        Hmcs_ii=harmonics{ii};

        if ~isrow(Hmcs_ii)
            Hmcs_ii=Hmcs_ii';
        end


        baseExp = exp(-1i*Hmcs_ii'*theta');
        switch optarg.method
            case 'angle'
                base = baseExp;
                x_int = theta;
            case 'mixed'
                base = baseExp.*omega';
                x_int = time;
        end


        Integral_baseSig             = cumtrapz(x_int,base.*Signal_ii,2)/2/pi;
        Integral_Decal               = zeros(size(Integral_baseSig));
        Integral_Decal(:,istart:end) = interp1(theta(istart:end)',Integral_baseSig(:,istart:end).',theta(istart:end)'-2*pi).';
        Integral_Decal(isnan(Integral_Decal))=0;
        phasor_v2                    = Integral_baseSig-Integral_Decal;


        %base d'analyse
        %BaseAna = baseExp.*omega';

        %integrale et calcul du phaseur
        %IntegraleBanaSig=cumtrapz(time,BaseAna.*Signal_ii,2)/2/pi;
        % IntegraleBanaSig=cumtrapz(time,BaseAna.*Signal_ii,2)./shift';

        %Integrale_Decal=interp1(theta,IntegraleBanaSig.',theta-2*pi);
        %Integrale_Decal(isnan(Integrale_Decal))=0;
        %
        % IntDecal_1=IntegraleBanaSig(:,IDX);
        % IntDecal_2=IntegraleBanaSig(:,IDX+1);
        % phasor_v=IntegraleBanaSig-IntDecal_1.*(1-IDX_frac')-IntDecal_2.*(IDX_frac');


        % IntDecal_1=IntegraleBanaSig(:,L);
        % IntDecal_2=IntegraleBanaSig(:,L+1);
        % phasor_v=IntegraleBanaSig-IntDecal_1.*(1-F')-IntDecal_2.*(F');

        %correction si necessaire, et decalage de 2pi
        % IntDecalee=IntegraleBanaSig(:,max(IDX-cor_i,1));
        % IntDecalee=IntegraleBanaSig(:,(IDX-cor_i));
        % phasor_v=IntegraleBanaSig-IntDecalee;
        %phasor_v2=IntegraleBanaSig-Integrale_Decal.';

        %integral with respect to theta
        %integral using int_theta(t)-2pi^theta(t) exp(i k \phi) x(t(phi))) dphi /2pi
        % IntegraleBanaSig=cumtrapz(theta,baseExp.*Signal_ii,2)/2/pi;

        %IntegralBexpSig_v3 = cumtrapz(theta,baseExp.*Signal_ii,2)/2/pi; %integral with respect to theta
        %interp
        %Integrale_Decal_v3=interp1(theta,IntegralBexpSig_v3.',theta-2*pi);
        %Integrale_Decal_v3(isnan(Integrale_Decal_v3))=0;
        %phasor_v3=IntegralBexpSig_v3-Integrale_Decal_v3.';




        phasor_cell{ii}             = phasor_v2;
        phasorStruct(ii).name       = NameSignals{ii};
        phasorStruct(ii).phasors    = phasor_v2;
        phasorStruct(ii).signal     = signals{ii};
        phasorStruct(ii).time       = time;
        phasorStruct(ii).harmonics  = Hmcs_ii;
        phasorStruct(ii).theta      = theta;
        phasorStruct(ii).omega      = omega;
        phasorStruct(ii).IDX        = IDX;


    catch e
        e.stack
        NameSignals(ii)
        error(e.message)
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Affichage
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if sum(PlotTAPRI)>0
    plotAngularSFT(phasorStruct,PlotTAPRI,"orientation",optarg.orientation,"plotDebut",optarg.plotDebut,"plotOmega",optarg.plotOmega,"xAxes",optarg.xAxes,"lang",optarg.plotlang)
end

    function [theta,time,w,sig,hm,nm_sig,PlotTAPRI,outcell,cor_i]= validate_input(theta,time,w,sig,hm,nm_sig,PlotTAPRI,plotdeb,plotw,options)
        outcell=true;
        if ~iscell(squeeze(sig))
            sig=squeeze(sig);
            if ~isvector(sig)
                if size(sig,1)==numel(time)
                    sig=num2cell(sig,1);

                elseif size(sig,2)==numel(time)
                    sig=num2cell(sig,2);
                else
                    error('error dim of sig')
                end

            else
                sig={squeeze(sig)};
                outcell=false;
            end
        end
        if ~iscell(squeeze(hm))
            hm={squeeze(hm)};
        end

        if ~iscell(nm_sig)
            nm_sig={nm_sig};
        end

        if numel(hm)==1
            myCell=cell(1,numel(sig));
            cell(hm);
            myCell(:) = hm;
            hm=myCell;
        end


        if isempty(nm_sig)
            nm_sig="signal "+string(1:numel(hm))';
        end


        if isempty(theta)
            theta=0;
        end

        switch numel(w)
            case 0
                if numel(theta)~=numel(time)
                    error('\omega cannot be empty with if provided \theta is not a time vector')
                end
                for iter_thetha =1:(numel(theta)-1)
                    w(iter_thetha) = (theta(iter_thetha+1)-theta(iter_thetha))/(time(iter_thetha+1)-time(iter_thetha));
                end
                w(numel(theta))=    w(numel(theta)-1);
            case 1
                w=ones(size(time))*w;
        end
        if numel(theta)==1
            theta=cumtrapz(time,w)+theta;
        end

        theta=unwrap(theta);

        w=squeeze(w);
        if ~iscolumn(w)
            w=w';
        end



        time=squeeze(time);
        if ~isrow(time)
            time=time';
        end
        theta=squeeze(theta);
        if ~iscolumn(theta)
            theta=theta';
        end


        cor_i=0;

    end

end


