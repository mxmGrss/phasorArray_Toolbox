function  TimeArray2plot(At,t,varg)
    %TIMEARRAY2PLOT Plot time-domain data from a 3D array.
    %
    %   TIMEARRAY2PLOT(At, t, varg) generates a plot of the time-domain data
    %   contained in the 3D array `At`. The data is assumed to be organized as
    %   a set of time-series signals along the third dimension of the array.
    %   Depending on the `explosed` option, each individual time-series is plotted
    %   in separate subplots or combined into a single plot.
    %
    %   Inputs:
    %     At       - (3D array) Time-domain data organized as a `N x M x P` matrix,
    %                where `N` and `M` represent the number of rows and columns,
    %                and `P` is the number of time steps.
    %     t        - (vector, optional) Time steps corresponding to the third dimension
    %                of `At`. Default: 1:P if omitted.
    %     varg     - (struct) Optional parameters:
    %                   'explosed' - (logical) If true (default), each component is
    %                                plotted in separate subplots. If false, all components
    %                                are plotted in a single combined plot.
    %                   'hold'     - (logical) If true, holds the current plot and adds 
    %                                the new plot on top of existing ones. Default: false.
    %
    %   Outputs:
    %     None. The function modifies the current figure to display the plot.
    %
    %   Behavior:
    %     - The function plots the values in `At` along the third dimension (time).
    %     - If `explosed` is true, each time-series is plotted in a separate subplot.
    %     - If `explosed` is false, all time-series are plotted in a single combined plot.
    %
    %   Example Usage:
    %     % Plot all time-series in separate subplots with holding enabled
    %     TimeArray2plot(At, t, struct('explosed', true, 'hold', true));
    %
    %     % Plot all time-series in a single combined plot
    %     TimeArray2plot(At, t, struct('explosed', false));
    %
    %     % Plot with custom time steps and holding disabled
    %     TimeArray2plot(At, t, struct('explosed', true, 'hold', false));
    %
    %   See also: plot, subplot, hold.
arguments
    At
    t=[]
    varg.explosed=true
    varg.hold=false
end

s=size(At);

if isempty(t)
    t=1:s(3);
end

if varg.explosed
    for ii=1:s(1)
        for jj=1:s(2)
            subplot(s(1),s(2),(ii-1)*s(2)+jj)
            if varg.hold
                hold on
            end
            plot(t,squeeze(At(ii,jj,:)))
            hold off
        end
    end
else
    if varg.hold
        hold on
    end
    plot(t,reshape(Mt,[],numel(t),1))
    hold off
end


end