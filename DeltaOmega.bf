//VERBOSITY_LEVEL = 0;
LoadFunctionLibrary("GrabBag");

SetDialogPrompt("Please specify a codon data file:\n");
//fprintf (stdout, "I've got a lovely bunch of " + PROMPT_FOR_FILE);
DataSet _ds_branch = ReadDataFile(PROMPT_FOR_FILE);
DataSetFilter _dsf_branch = CreateFilter(_ds_branch,3,"","",GeneticCodeExclusions);
GetInformation(infotainer, _dsf_branch);
fprintf (stdout, infotainer);
//HarvestFrequencies(_nuc3_branch, dsf, 3, 1, 1);
//nucCF = CF3x4(_nuc3_branch, GeneticCodeExclusions);
fscanf (stdin, "Number", branchNumber);

fprintf (stdout, "loaded and executed!\n");




return 0;
