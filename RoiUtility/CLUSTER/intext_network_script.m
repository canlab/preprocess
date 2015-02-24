

k = 1;
[subcl,out] = cluster_network(union2_cl,k,covs,'EO_network'); cd ..; close all;
EO_network_cl = subcl;
EO_network_out = out;
if mean(out.pc(:,1)) < 0, network_PCs(:,k) = -out.pc(:,1);,else,network_PCs(:,k) = out.pc(:,1);,end
if mean(out.pc(:,1)) < 0, network_scores1(:,k) = -out.score(:,1);,else,network_scores1(:,k) = out.score(:,1);,end
if mean(out.pc(:,1)) < 0, network_scores2(:,k) = -out.score(:,2);,else,network_scores2(:,k) = out.score(:,2);,end

k = 2;
[subcl,out] = cluster_network(union2_cl,k,covs,'EA_network'); cd ..; close all;
EA_network_cl = subcl;
EA_network_out = out;
if mean(out.pc(:,1)) < 0, network_PCs(:,k) = -out.pc(:,1);,else,network_PCs(:,k) = out.pc(:,1);,end
if mean(out.pc(:,1)) < 0, network_scores1(:,k) = -out.score(:,1);,else,network_scores1(:,k) = out.score(:,1);,end
if mean(out.pc(:,1)) < 0, network_scores2(:,k) = -out.score(:,2);,else,network_scores2(:,k) = out.score(:,2);,end

k = 4;
[subcl,out] = cluster_network(union2_cl,k,covs,'IO_network'); cd ..; close all;
IO_network_cl = subcl;
IO_network_out = out;
if mean(out.pc(:,1)) < 0, network_PCs(:,k) = -out.pc(:,1);,else,network_PCs(:,k) = out.pc(:,1);,end
if mean(out.pc(:,1)) < 0, network_scores1(:,k) = -out.score(:,1);,else,network_scores1(:,k) = out.score(:,1);,end
if mean(out.pc(:,1)) < 0, network_scores2(:,k) = -out.score(:,2);,else,network_scores2(:,k) = out.score(:,2);,end

k = 5;
[subcl,out] = cluster_network(union2_cl,k,covs,'IA_network'); cd ..; close all;
IA_network_cl = subcl;
IA_network_out = out;
if mean(out.pc(:,1)) < 0, network_PCs(:,k) = -out.pc(:,1);,else,network_PCs(:,k) = out.pc(:,1);,end
if mean(out.pc(:,1)) < 0, network_scores1(:,k) = -out.score(:,1);,else,network_scores1(:,k) = out.score(:,1);,end
if mean(out.pc(:,1)) < 0, network_scores2(:,k) = -out.score(:,2);,else,network_scores2(:,k) = out.score(:,2);,end

k = 6;
[subcl,out] = cluster_network(union2_cl,k,covs,'IINT_network'); cd ..; close all;
IINT_network_cl = subcl;
IINT_network_out = out;
if mean(out.pc(:,1)) < 0, network_PCs(:,k) = -out.pc(:,1);,else,network_PCs(:,k) = out.pc(:,1);,end
if mean(out.pc(:,1)) < 0, network_scores1(:,k) = -out.score(:,1);,else,network_scores1(:,k) = out.score(:,1);,end
if mean(out.pc(:,1)) < 0, network_scores2(:,k) = -out.score(:,2);,else,network_scores2(:,k) = out.score(:,2);,end

% FLIP COMP2 SCORES FOR IA ONLY, REST ARE ALREADY FLIPPED
% BECAUSE THEY HAVE NEG WEIGHTS ON FIRST PC
network_scores2(:,5) = -network_scores2(:,5);





