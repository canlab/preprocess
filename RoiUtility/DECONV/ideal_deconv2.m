function ideal_deconv2(DX,tp,TR)
% function ideal_deconv2(DX,tp,TR)
% tests deconvolution matrix directly
%
% Tor Wager, 10/24/01

mycolors = {'ro-' 'b^--' 'go:' 'y^--' 'mo-' 'r^-' 'g^-' 'b^-' 'k^-' 'm^-' 'y^-'};

hrf = spm_hrf(TR);
hrf = hrf ./ max(hrf);

% -------------------------------------------------------------------
% * set up single response
% -------------------------------------------------------------------

ideal_data = conv(hrf,DX(:,1));
ideal_data = ideal_data(1:size(DX,1));

b = pinv(DX) * ideal_data;

% -------------------------------------------------------------------
% * plot single response
% -------------------------------------------------------------------
figure;

subplot(1,2,1)
plot(hrf)
hold on; 
%plot(b(1:end-1),'rs-')
index = 1;
for i = 1:tp:size(DX,2)-1
    plot(b(i:i+tp-1),mycolors{index})
    myleg{index} = ['Condition' num2str(index)];
    index = index + 1;
end
grid on
title('Single condition response')

% -------------------------------------------------------------------
% * set up multiple response
% -------------------------------------------------------------------

ideal_data = [];
for i = 1:tp:size(DX,2)-1
    mydata = conv(hrf,DX(:,i));
    mydata = mydata(1:size(DX,1));
    ideal_data(:,i) = mydata;
end
ideal_data = sum(ideal_data,2);

b = pinv(DX) * ideal_data;


% -------------------------------------------------------------------
% * plot multiple response
% -------------------------------------------------------------------
subplot(1,2,2)


plot(hrf)
hold on; 
index = 1;
for i = 1:tp:size(DX,2)-1
    plot(b(i:i+tp-1),mycolors{index})
    myleg{index} = ['Condition' num2str(index)];
    index = index + 1;
end

grid on
title('Multiple (all) condition response')