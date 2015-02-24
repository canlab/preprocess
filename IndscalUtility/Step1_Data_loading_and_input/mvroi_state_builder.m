function states = mvroi_state_builder(numimgs,numsubs,varargin)
% states = mvroi_state_builder(numimgs,numsubs,[counterbalance vector])
%
% numimgs: number of time points
% optional counterbalance: vector of 1 / -1 values for each subject



z = zeros(numimgs,1);

disp('This will only work if all subjects have the same design...')
disp('Or if you enter a vector of counterbalanced orders across subjects.')

disp('Enter onsets in TRs for each condition, with 0 indexing the first image.')
disp('When all conditions are entered, press return.')
disp(' ')

go = 1; ind = 1;
while go
    
    st = input(['Enter onsets for condition ' num2str(ind) ' in images (0 is image 1): ' ]);
    
    en = input(['Enter offsets for condition ' num2str(ind) ' in images (0 is image 1): ' ]);
    
    if ~isempty(st)
        for i = 1:length(st)
            z(st(i):en(i)) = ind;
        end
    
    else
        go = 0;
    end
    
    ind = ind + 1;
end

states{1} = z;

for i = 1:numsubs
    states{i} = states{1};
end

if length(varargin) > 0
    ord = varargin{1};
    for i = 1:length(states)
        
        if ord(i) < 0,
            states{i} = -1 * states{i} + (max(states{i}) + 1);
        end
    end
end


return