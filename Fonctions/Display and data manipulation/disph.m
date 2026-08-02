function disph(inputArg1,nvp)
arguments (Repeating)
    inputArg1
end
arguments
    nvp.inhib=false
end
if ~nvp.inhib
    if isempty(inputArg1)
        disp(strcat(string(datetime(),'HH:mm:ss.SSS')))
    else
        inputArg1=cellfun(@(x) string(x),inputArg1);
        disp(strcat(string(datetime(),'HH:mm:ss.SSS')," : ",[inputArg1{:}]))
    end
end