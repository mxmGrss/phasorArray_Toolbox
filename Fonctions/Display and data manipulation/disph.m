function disph(inputArg1,opts)
arguments (Repeating)
    inputArg1
end
arguments
    opts.inhib=false
end
if ~opts.inhib
    if isempty(inputArg1)
        disp(strcat(string(datetime(),'HH:mm:ss.SSS')))
    else
        inputArg1=cellfun(@(x) string(x),inputArg1);
        disp(strcat(string(datetime(),'HH:mm:ss.SSS')," : ",[inputArg1{:}]))
    end
end