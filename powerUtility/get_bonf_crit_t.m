for df = 5:5:40
	for n = 1:5:30000,
		tc(n,df) = tinv(1-(.05/n),df);
	end
end

figure; hold on; set(gcf,'Color','w')
plot(tc)

df = 5:5:40;
leg = mat2cell(df,1,ones(size(df)));
for i = 1:length(leg), leg{i} = num2str(leg{i});, end
legend(leg)