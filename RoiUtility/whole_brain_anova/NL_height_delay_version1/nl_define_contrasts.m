function DX = nl_define_contrasts(DX)
% function DX = nl_define_contrasts(DX)
% Tor Wager
% Sets up matrix of contrasts and descriptions
% Prompting the user for values
%
% DX is structure from EXPT.DX
%
% OR, to use pre-existing ones, try
% mycons = cat(1,EXPT.SUBJECT.contrasts{:});
% EXPT.DX.contrasts = mycons(:,EXPT.DX.regsofinterest);
% EXPT.DX.connames = EXPT.SUBJECT.connames;
%
mycon = 1;
coni = 1;
conlen = length(DX.regsofinterest);


disp(['Define contrasts over these conditions (just return to quit):'])
fprintf(1,'\n')
for i = 1:conlen
    fprintf(1,'%s ',DX.dxnames{i})
end
fprintf(1,'\n')

while ~isempty(mycon) 
    
    mycon = input(['Enter contrast values for contrast ' num2str(coni) ' in brackets []. ']);
    
    if length(mycon) == conlen
        
        if sum(mycon) == 0, 
            
            % enter the contrast
            DX.contrasts(coni,:) = mycon;
        
            % define contrast name
            myname = [];
            for i = 1:conlen
                if mycon(i) > 0, mysign = '+'; elseif mycon(i) < 0, mysign = '-';, end
            
                if ~(mycon(i) == 0)
                    myname = [myname ' ' mysign DX.dxnames{i}];
                end
            end
        
            %DX.connames{coni} = input(['Enter description or contrast name','%s']);
            DX.connames{coni} = myname;
            
            coni = coni + 1;
            fprintf(1,'%s\n',myname)
 
        else
    
            disp(['Contrast does not sum to zero.  Please re-enter. '])
        
        end
        
    elseif ~isempty(mycon)
        
        disp(['Contrast length should be ' num2str(conlen) ' values; Please re-enter.'])
        
    end
    
end

return
