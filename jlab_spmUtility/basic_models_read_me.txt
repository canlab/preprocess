this folder contains a few functions to implement
random effects analyses (and several forms of basic models)
in SPM without the GUI.

Some functions are modified from the SPM99 distribution
so that instead of asking for user input, the function
reads the values for each question from an input struct
variable, with values stored in INPUT.variable format.

I didn't modify the functions to deal with all possible
basic models - only the ones I've been using.  So there's
no guarantee that things will work perfectly for you without
some work.  But they should at least serve as a starting point.

getfiles2 is a function that lists the file names, much as does
spm_list_files.m, and it works for windows.  I don't know about
unix/linux.  You may have to make your own function if this one
doesn't work.  The output should be filenames with full pathnames
the same as in spm_list_files.m.

Here's how the functions work.  The top level is an example
script, biman5_run_randfx.m, which loops through contrasts in
a study and runs the analysis for each contrast.  You'll have
to create a version of this for your study.

This script calls biman5_random_effects, another script.
Edit this and put in your own values in the user-entered section.
Another example is tor_basic_models, which accomplishes about
the same thing.

Those scripts create a structure called TOR with all the 
input values that you want for the analysis.  That structure
is passed into 
	tor_spm_spm_ui.m, 	which does the rfx analysis, and
	estimateContrasts.m 	a script for estimating the rfx contrast
	tor_spm_results_ui.m	gets the results and prints a table

You'll have to know Matlab scripting, and again - they're 
basically as-is, with no guarantee that they're right or anything.
I just tried to take out the questions, and leave the spm code
as intact as possible.


Tor Wager
University of Michigan