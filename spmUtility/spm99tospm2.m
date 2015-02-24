% SPM2 = spm99tospm2(SPM99)

function SPM2 = spm99tospm2(SPM99)
    SPM2 = SPM99;
    
    assoc_field('xX.RT', 'xY.RT');
    assoc_field('xX.dt', 'xBF.dt');

    function assoc_field(SPM99_field, SPM2_field)
        eval(sprintf('SPM2.%s = SPM99.%s;', SPM2_field, SPM99_field));
        wh_dots = strfind(SPM99_field, '.');
        prior_fields = SPM99_field(1:wh_dots(end)-1);
        last_field = SPM99_field(wh_dots(end)+1:end);
        eval(sprintf('rmfield(SPM2.%s, ''%s'');', prior_fields, last_field));
    end
end