function [A] = rand_phasor(nx,ny,h,arg)
    % RAND_PHASOR Generate a random PhasorArray with a specified structure.
    %
    %   A = RAND_PHASOR(NX, NY, H, <name-value arguments>) creates a 3D array 
    %   of complex phasors representing a time-dependent matrix with harmonic components.
    %
    %   Inputs:
    %     NX  - (integer) Number of rows.
    %     NY  - (integer) Number of columns.
    %     H   - (integer) Number of harmonics.
    %
    %   Name-Value Arguments:
    %     'time_structure' (char) - Structure of the generated phasors (default: 'real').
    %                               Options: 'real', 'symmetric', 'antisymmetric', 
    %                               'sdp', 'hurwitz', 'hermitian', 'Q-spec', 'cmplx', 
    %                               'retroHermitian'.
    %     'hurwitzeig' (vector)   - Desired eigenvalues for Hurwitz matrices (default: random).
    %     'T' (scalar)            - Period of the system (default: 1).
    %     'Q' (matrix)            - Matrix used in Q-spectrum methods (default: `[]`).
    %     'output' (char)         - Output format: 'NDarray' (default) or 'PhasorArray'.
    %     'average_power_decay'   - Decay rate of harmonic amplitudes (default: 2).
    %
    %   Outputs:
    %     A - (3D array or PhasorArray) Randomly generated PhasorArray.
    %
    %   Notes:
    %     - Generates phasors with complex Gaussian components.
    %     - The decay factor controls power reduction of higher-order harmonics.
    %     - If `time_structure='hurwitz'`, matrix has Hurwitz eigenvalues in the HArmonic domain (T(A)-N).
    %     - If `time_structure='sdp'`, matrix satisfies semidefinite constraints.
    %
    %   See also: PhasorArray, Sylv_harmonique, randomPhasorArrayWithPole.
arguments
    nx
    ny
    h
    arg.time_structure {mustBeMember(arg.time_structure,{'real','symetric','antisymetric','sdp','hurwitz','Q-spec','hermitian','cmplx','retroHermitian'})} ='real'
    arg.hurwitzeig=-rand(nx,1)'.*(10.^(rand(nx,1)'+(1:nx)))-10
    arg.T=1
    arg.Q=[]
    arg.output {mustBeMember(arg.output,{'NDarray','PhasorArray'})} ='NDarray' %or PhasorArray
    arg.average_power_decay=2;
end
A=zeros(nx,ny,2*h+1);

decay=arg.average_power_decay;

switch arg.time_structure
    case 'real'
        for hh=0:h
            A(:,:,h+1+hh)=((rand(nx,ny)-0.5)*2+1i*(rand(nx,ny)-0.5)*2*(hh~=0)).*(rand(nx,ny)*8)*1/(hh+1)^decay*1/arg.T;
            A(:,:,h+1-hh)=conj(A(:,:,h+1+hh));
        end
    
    case 'symetric'
    for hh=0:h
         A(:,:,h+1+hh)=((rand(nx,ny)-0.5)*2+1i*(rand(nx,ny)-0.5)*2*(hh~=0)).*(rand(nx,ny)*8)*1/(hh+1)^decay*1/arg.T;
         A(:,:,h+1+hh)=(A(:,:,h+1+hh).'+A(:,:,h+1+hh))/2;
         A(:,:,h+1-hh)=conj(A(:,:,h+1+hh));
    end
    case 'antisymetric'
    for hh=0:h
         A(:,:,h+1+hh)=((rand(nx,ny)-0.5)*2+1i*(rand(nx,ny)-0.5)*2*(hh~=0)).*(rand(nx,ny)*8)*1/(hh+1)^decay*1/arg.T;
         A(:,:,h+1+hh)=(A(:,:,h+1+hh).'-A(:,:,h+1+hh))/2;
         A(:,:,h+1-hh)=conj(A(:,:,h+1+hh));
    end

    case 'sdp'
        while(~isempty(find(real(eig(array2TBlocks(A,10*h)))<1e-2,1)))
        for hh=0:h
             A(:,:,h+1+hh)=((rand(nx,nx)-0.5)*2+1i*(rand(nx,nx)-0.5)*2*(hh~=0)).*(rand(nx,nx))*8*1/(hh+1)^decay.*(ones(nx)+eye(nx)*nx*(hh==0))*1/sqrt(arg.T);
    %        A(:,:,h+1+hh)=((A(:,:,h+1+hh).'+A(:,:,h+1+hh))/2)*((A(:,:,h+1+hh).'+A(:,:,h+1+hh))/2)';
             A(:,:,h+1+hh)=((A(:,:,h+1+hh).'+A(:,:,h+1+hh))/2);
             A(:,:,h+1+hh)=A(:,:,h+1+hh)*A(:,:,h+1+hh).';
             A(:,:,h+1-hh)=conj((A(:,:,h+1+hh)));
        end
        end
    case 'hurwitz'
%         arg.hurwitzeig
        hbi=h*4;
       W=rand_phasor(nx,nx,h,"time_structure",'sdp');
       N=N_tb(nx,hbi,arg.T);
%        N=kron(eye(nx),diag(1i*(-hbi:hbi)));
       Qhm=array2TBlocks(diag(arg.hurwitzeig),2*hbi);
       Whm=array2TBlocks(W,2*hbi);
       Amn=Whm*(Qhm-N)*Whm^-1;
       Ahmd=Amn+N;
       A=TB2array(Ahmd,nx);
       A=A(:,:,(end+1)/2+(-h:h));
    case 'Q-spec'
        hbi=h*1;
       W=rand_phasor(nx,nx,h,"time_structure",'sdp');
       N=N_tb(nx,hbi,arg.T);
%        N=kron(eye(nx),diag(1i*(-hbi:hbi)));
       Qhm=array2TBlocks(arg.Q,2*hbi);
       Whm=array2TBlocks(W,2*hbi);
       Amn=Whm*(Qhm-N)*Whm^-1;
       Ahmd=Amn+N;
       A=TB2array(Ahmd,nx);
%        A=A(:,:,(end+1)/2+(-h:h));
    case 'retroHermitian'
        for hh=0:h
            A(:,:,h+1+hh)=((rand(nx,ny)-0.5)*2+1i*(rand(nx,ny)-0.5)*2).*(rand(nx,ny)*8)*1/(hh+1)^decay*1/arg.T;
            A(:,:,h+1+hh)=(A(:,:,h+1+hh)'+A(:,:,h+1+hh))/2;
            A(:,:,h+1-hh)=((rand(nx,ny)-0.5)*2+1i*(rand(nx,ny)-0.5)*2).*(rand(nx,ny)*8)*1/(hh+1)^decay*1/arg.T;
            A(:,:,h+1-hh)=(A(:,:,h+1-hh)'+A(:,:,h+1-hh))/2;
        end
    case 'hermitian'
        for hh=1:h
            A(:,:,h+1+hh)=((rand(nx,ny)-0.5)*2+1i*(rand(nx,ny)-0.5)*2).*(rand(nx,ny)*8)*1/(hh+1)^decay*1/arg.T;
            A(:,:,h+1-hh)=A(:,:,h+1+hh)';
        end
        A(:,:,h+1)=((rand(nx,ny)-0.5)*2+1i*(rand(nx,ny)-0.5)*2).*(rand(nx,ny)*8)*1/(hh+1)^decay*1/arg.T;
        A(:,:,h+1)=(A(:,:,h+1)+A(:,:,h+1)')/2;

    otherwise
        for hh=0:h
            A(:,:,h+1+hh)=((rand(nx,ny)-0.5)*2+1i*(rand(nx,ny)-0.5)*2).*(rand(nx,ny)*2)*1/(hh+1)^decay*1/arg.T;
            A(:,:,h+1-hh)=((rand(nx,ny)-0.5)*2+1i*(rand(nx,ny)-0.5)*2).*(rand(nx,ny)*2)*1/(hh+1)^decay*1/arg.T;
        end
end

switch arg.output
    case 'PhasorArray'
        A=PhasorArray(A);
    otherwise
end

end