function [src]= getParentList(history, node)

    %source virus list
    
    src = [];
    src = node;
    [s] = getParent(history,node);
    while s ~= 0
        src(end+1)=s;
        [s]=getParent(history,src(end));  
    end
    src(end+1)=s;
    
    function [src epsb] = getParent(history, node)
        src = history(node);
    end
end