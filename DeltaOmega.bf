//VERBOSITY_LEVEL = 0;
LoadFunctionLibrary("GrabBag");

SetDialogPrompt("Please specify a codon data file:\n");
//fprintf (stdout, "I've got a lovely bunch of " + PROMPT_FOR_FILE);
DataSet _ds_branch = ReadDataFile(PROMPT_FOR_FILE);
fscanf (stdin, "Number", branchNumber);
DataSetFilter _dsf_branch = CreateFilter(_ds_branch,3,"","",GeneticCodeExclusions);
LoadFunctionLibrary("queryTree");
bNames = BranchName (givenTree, -1);

LikelihoodFunction	branch_base_LF	 = (dsf, givenTree);
fprintf(stdout, "in DeltaOmega, optimizing: \n");
ref = "givenTree." + bNames[branchNumber];
ExecuteCommands ("givenTree." + bNames[branchNumber] + ".omega4 = 1.1;");
ExecuteCommands ("givenTree." + bNames[branchNumber] + ".Paux3 = 0.1;");
Optimize(opt_branch_base_LF, branch_base_LF);

lfOut	= csvFilePath + ".optTree.DO." + branchNumber + ".fit";
LIKELIHOOD_FUNCTION_OUTPUT = 7;
fprintf (lfOut, CLEAR_FILE, branch_base_LF);
LIKELIHOOD_FUNCTION_OUTPUT = 2;

//lfOut	= csvFilePath + ".optTreeRes.DO." + branchNumber + ".fit";
//LIKELIHOOD_FUNCTION_OUTPUT = 7;
//fprintf (lfOut, CLEAR_FILE, opt_branch_base_LF);
//LIKELIHOOD_FUNCTION_OUTPUT = 2;


// Print stuff
GetInformation(infotainer, _dsf_branch);
fprintf (stdout, infotainer);
//HarvestFrequencies(_nuc3_branch, dsf, 3, 1, 1);
//nucCF = CF3x4(_nuc3_branch, GeneticCodeExclusions);

fprintf (stdout, "loaded and executed!\n");




return 0;
