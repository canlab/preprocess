
for msub = EXPT.subjects(1:5)

	str = (['[manip_' msub{1} ',test_' msub{1} '] = plac_make_c(msub{1});']);
	eval(str)

end