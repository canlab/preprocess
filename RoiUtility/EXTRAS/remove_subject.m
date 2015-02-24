function EXPT = remove_subject(EXPT,omit)
%EXPT = remove_subject(EXPT,omit number)

sd = EXPT.subjects;

N = fieldnames(EXPT);

for i = 1:length(N)
    
    eval(['f = EXPT.' N{i} ';']);
    
    if length(f) == length(sd)
        f(omit) = [];
    end

    if size(f,1) == length(sd) 
        f(omit,:) = [];      
    end
        
    eval(['EXPT.' N{i} '=f;']);
    
end

if isfield(EXPT,'SNPM')
    if isfield(EXPT.SNPM,'P')
    
    for i = 1:length(EXPT.SNPM.P)
        EXPT.SNPM.P{i}(omit,:) = [];
    end
    end
end

    

if isfield(EXPT,'SUBJECT')
    
    N = fieldnames(EXPT.SUBJECT);

    for i = 1:length(N)
    
        eval(['f = EXPT.SUBJECT.' N{i} ';']);
    
        if length(f) == length(sd)
            f(omit) = [];
        end
        
        if size(f,1) == length(sd)
            f(omit,:) = [];
        end
   
        eval(['EXPT.SUBJECT.' N{i} '=f;']);
    
    end

end


if isfield(EXPT,'FILES')
    
    N = fieldnames(EXPT.FILES);

    for i = 1:length(N)
    
        eval(['f = EXPT.FILES.' N{i} ';']);
    
        if length(f) == length(sd)
            f(omit) = [];
        end
        
        if size(f,1) == length(sd)
            f(omit,:) = [];
        end
   
        eval(['EXPT.FILES.' N{i} '=f;']);
    
    end

end


if isfield(EXPT,'FILT')
    
    N = fieldnames(EXPT.FILT);

    for i = 1:length(N)
    
        eval(['f = EXPT.FILT.' N{i} ';']);
    
        if length(f) == length(sd)
            f(omit) = [];
        end
        
        if size(f,1) == length(sd)
            f(omit,:) = [];
        end
   
        eval(['EXPT.FILT.' N{i} '=f;']);
    
    end

end


if isfield(EXPT,'FIR')
    
    N = fieldnames(EXPT.FIR);

    for i = 1:length(N)
    
        eval(['f = EXPT.FIR.' N{i} ';']);
    
        if length(f) == length(sd)
            f(omit) = [];
        end
        
        if size(f,1) == length(sd)
            f(omit,:) = [];
        end
   
        eval(['EXPT.FIR.' N{i} '=f;']);
    
    end

end