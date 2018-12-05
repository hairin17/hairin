subs = [n ... ]

    for s = 1:numel(subs);
    subj = subs(s)
    dirname = [num2str(subj,'%03d')]
    mkdir(dirname);
    end
    
