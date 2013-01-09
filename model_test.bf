VERBOSITY_LEVEL				 = 0;

skipCodeSelectionStep 		= 0;
LoadFunctionLibrary("chooseGeneticCode");   // Where are these libraries (for
                                            //reference/use)

LoadFunctionLibrary("GrabBag");
LoadFunctionLibrary("dSdNTreeTools");
LoadFunctionLibrary("CF3x4");
LoadFunctionLibrary("BranchSiteTemplate");
//LoadFunctionLibrary("DeltaOmega.bf");
LoadFunctionLibrary("AddRateClass.bf");


modelList = {};


DataSet 			ds 				= ReadDataFile(PROMPT_FOR_FILE);
DataSetFilter 		dsf 			= CreateFilter(ds,3,"","",GeneticCodeExclusions);

GetInformation(dsf_seq_array, dsf);
taxa_count = Columns(dsf_seq_array);
algn_len = Abs(dsf_seq_array[0]);


// Are these the nucleotide frequences?
// If so, are they harvested but left in dsf? Or are they now in nuc3?
HarvestFrequencies	(nuc3, dsf, 3, 1, 1);

// I don't know what this is yet, but nuc3 is passed in as an argument here.
nucCF						= CF3x4	(nuc3, GeneticCodeExclusions);

// Why is this one separate? Matrix Local with syn and nonsyn? hmmm
PopulateModelMatrix			  ("MGMatrixLocal",  nucCF, "syn", "", "nonsyn");

// So we harvested frequencies earlier (character frequencies), now we're
// building codon frequencies from the nucleotide frequencies?
codon3x4					= BuildCodonFrequencies (nucCF);
// I'm not sure I get this syntax. Is this some sort of constructor, or
// tuple?
Model		MGL				= (MGMatrixLocal, codon3x4, 0);

modelList[1] = "MGL";

LoadFunctionLibrary			  ("queryTree");

// I know a few of these words (they just print thinks, as far as I can tell)
SetDialogPrompt ("Save analysis results to");
//fprintf (PROMPT_FOR_FILE, CLEAR_FILE, KEEP_OPEN,"Branch,Mean_dNdS,Omega1,P1,Omega2,P2,Omega3,P3,LRT,p,p_Holm");
fprintf (PROMPT_FOR_FILE, CLEAR_FILE, KEEP_OPEN,"Branch,Mean_dNdS,LRT,p,p_Holm");
csvFilePath = LAST_FILE_PATH;

fprintf 					  (stdout, "[PHASE 0] Fitting the local MG94 (no site-to-site variation) to obtain initial parameter estimates\n");

// Again, is this syntax for a constructor or a function or a tuple?
// dsf = DataSetFilter
// givenTree doesn't yet exist, does it? Is it a global variable?
// No, it does not. Rather it is a hidden global variable that is populated
// by some other process.
// Yes, givenTree is a hidden global variable that shows up when you
// "queryTree"
// So yeah, it doesn't appear that Optimization ever really happens in the
// code or in execution. Thats chill...
LikelihoodFunction	base_LF	 = (dsf, givenTree);


//Optimize					  (res_base,base_LF);
//GetInformation(treeInformation, givenTree);
//fprintf (csvFilePath + ".givenTree.tree", CLEAR_FILE, treeInformation);

// Not sure what LIKELIHOOD_FUNCTION_OUTPUT is about. I'm assuming it has to
// do with printing, however
//lfOut	= csvFilePath + ".mglocal.fit";
//LIKELIHOOD_FUNCTION_OUTPUT = 7;
//fprintf (lfOut, CLEAR_FILE, base_LF);
//LIKELIHOOD_FUNCTION_OUTPUT = 2;

//localLL						 = res_base[1][0];
//localParams					 = res_base[1][1] + 9;

LoadFunctionLibrary			 ("DescriptiveStatistics");

// Function from DescriptiveStatistics?
totalBranchCount			 = BranchCount(givenTree) + TipCount (givenTree);

//GetInformation	   		   (varNames, "givenTree\\..*\\.omega1");
//localOmegaValues			 = {totalBranchCount,1}["Eval(varNames[_MATRIX_ELEMENT_ROW_])"];

// What now is this syntax?
pValueByBranch				  = {totalBranchCount,10};

// Function from DescriptiveStatistics?
bNames						  = BranchName (givenTree, -1);


// not sure I get why things are called pValueByBranch...
// So this looks to be an important processing point, wherein synonymous and
// nonsynonymous rates are used to determine a a dn/ds ration (omega) for the
// branch.
for (k = 0; k < totalBranchCount; k = k+1)
{
	srate  = Eval ("givenTree." + bNames[k] + ".syn");
	nsrate = Eval ("givenTree." + bNames[k] + ".nonsyn");
	if (srate > 0)
	{
		pValueByBranch [k][0] = Min (10, nsrate/srate);
	}
	else
	{
		pValueByBranch [k][0] = 10;
	}
}

// Why are we negative one indexing here? All values? The Last value?
omegaStats					 = GatherDescriptiveStats (pValueByBranch[-1][0]);

fprintf						 (stdout, "\nLog L = ", localLL, " with ", localParams, " degrees of freedom\n");

//PrintDescriptiveStats		 ("Branch omega values", omegaStats);

Tree						   mixtureTree = treeString;

// Say what? are we saving parameters into temptorary variables or something?
//ReplicateConstraint 		  ("this1.?.t:=this2.?.syn",mixtureTreeG,givenTree);
ReplicateConstraint 		  ("this1.?.t:=this2.?.syn",mixtureTree,givenTree);

// And then clearing them...
ClearConstraints			  (mixtureTree);
//ClearConstraints			  (mixtureTreeG);

// Some options
ASSUME_REVERSIBLE_MODELS	  = 1;

VERBOSITY_LEVEL               = 1;

// Yay new LF
// So this is the likelihood function. But there is no "Optimize()" call...
// yet. There wasn't earlier either, but there is one later that takes this
// as the second parameter (look into the Optimize function)
LikelihoodFunction three_LF  = (dsf,mixtureTree);

//-------------------------------------------------------------------------
// So at this point we have the MG94 model. Now we will go through the
// branches and add&optimize rate classes.


lfOut	= csvFilePath + ".tree.fit";
LIKELIHOOD_FUNCTION_OUTPUT = 7;
fprintf (lfOut, CLEAR_FILE, three_LF);
LIKELIHOOD_FUNCTION_OUTPUT = 2;

Optimize (res_three_LF,three_LF);
fprintf(stdout, "\n");

lfOut	= csvFilePath + ".optTree.fit";
LIKELIHOOD_FUNCTION_OUTPUT = 7;
fprintf (lfOut, CLEAR_FILE, three_LF);
LIKELIHOOD_FUNCTION_OUTPUT = 2;

//for (k = 0; k < totalBranchCount; k = k+1)

iter_likelihood = 0; // Determined after Optimize
//iter_samples = taxa_count * algn_len; // Known, for now the # of codons * the number of sites
iter_samples = algn_len;
iter_parameters = 0; // Known, but can be determined after Optimize
init_parameters = 0;

GetInformation(dsf_seq_array, dsf);
//fprintf(stdout, dsf_seq_array);

Export(three_LF_bak, three_LF);

origRes = res_three_LF[1][0] - 1.0;
orig_likelihood = res_three_LF[1][0];
orig_parameters = res_three_LF[1][1];
orig_bic = calcBIC(orig_likelihood, orig_parameters, iter_samples);

best_models = {};

for (taxI = 0; taxI < totalBranchCount; taxI = taxI + 1)
{
    best_models[taxI] = 1;
}

//fprintf(stdout, "Branch Names:\n");
//fprintf(stdout, bNames);
//fprintf(stdout, "\n");
//fprintf(stdout, "Total Branch Count: ");
//fprintf(stdout, totalBranchCount);
//fprintf(stdout, "\n");

for (branchI = 0; branchI < totalBranchCount; branchI = branchI + 1)
{

    fprintf(stdout, "\n");
    fprintf(stdout, "New branch: \n");
    lastRes = origRes;
    omegaNumber = 1;

    better_bic = 1;
    //last_bic = 100000000;
    last_bic = orig_bic;

    bic_run_count = 0;

    //while (res_three_LF[1][0] > lastRes)
    //for (omegaNumber = 1; omegaNumber < 3; omegaNumber = omegaNumber + 1)
    //while (last_BIC < 0 || iter_bic < last_BIC)
    while (better_bic == 1)
    {
        //fprintf(stdout, "\n");
        addRate2Branch("three_LF", nucCF, bNames[branchI], "MGL", modelList, algn_len, 1);
        //fprintf(stdout, "Current Model list:\n");
        //fprintf(stdout, modelList);
        //fprintf(stdout, "\n");
        omegaNumber = omegaNumber + 1;

        lfOut	= csvFilePath + ".treePlusRate." + bNames[branchI] + "." + omegaNumber + ".fit";
        LIKELIHOOD_FUNCTION_OUTPUT = 7;
        fprintf (lfOut, CLEAR_FILE, three_LF);
        LIKELIHOOD_FUNCTION_OUTPUT = 2;

        VERBOSITY_LEVEL = 10; // 10 prints EVERYTHING
        //VERBOSITY_LEVEL = 0;

        Optimize (res_three_LF,three_LF);
        fprintf(stdout, "\n");

        iter_likelihood = res_three_LF[1][0];
        //iter_parameters = res_three_LF[1][1];
        if (branchI == 0 && omegaNumber == 2)
        {
            init_parameters = res_three_LF[1][1];
        }
        iter_parameters = init_parameters + (2 * omegaNumber) + 2;

        iter_bic = calcBIC(iter_likelihood, iter_parameters, iter_samples);
        fprintf(stdout, "This iterations likelihood: " + iter_likelihood);
        fprintf(stdout, "\n");
        fprintf(stdout, "This iterations parameter count: " + iter_parameters);
        fprintf(stdout, "\n");
        fprintf(stdout, "This iterations sample count: " + iter_samples);
        fprintf(stdout, "\n");
        fprintf(stdout, "This iterations BIC: " + iter_bic);
        fprintf(stdout, "\n");

        lfOut	= csvFilePath + ".optTreePlusRate." + bNames[branchI] + "." + omegaNumber + ".fit";
        LIKELIHOOD_FUNCTION_OUTPUT = 7;
        fprintf (lfOut, CLEAR_FILE, three_LF);
        LIKELIHOOD_FUNCTION_OUTPUT = 2;

        lastRes = res_three_LF[1][0];
        if (iter_bic == last_bic)
        {
            bic_run_count = bic_run_count + 1;
        }
        if ((iter_bic > last_bic) || (bic_run_count > 2))
        {
            fprintf(stdout, "Done with this branch...\n");
            better_bic = 0;
        }
        else
        {
            best_models[branchI] = best_models[branchI] + 1;
        }
        last_bic = iter_bic;
    }
    fprintf(stdout, "\n");
}
fprintf(stdout, "\n");
fprintf(stdout, best_models);
fprintf(stdout, "\n");

REPLACE_TREE_STRUCTURE = 1;

for (branchI = 0; branchI < totalBranchCount; branchI = branchI + 1)
{
    fprintf(stdout, "Adding " + (best_models[branchI]-1) + " to:\n");
    fprintf(stdout, bNames[branchI]);
    fprintf(stdout, "\n");
    for (modelI = 1; modelI < best_models[branchI]; modelI = modelI + 1)
    {
        fprintf(stdout, "A new model!\n");
        if (branchI == 0)
        {
            addRate2Branch("three_LF", nucCF, bNames[branchI], "MGL", modelList, algn_len, 1);
        }
        else
        {
            addRate2BranchAdvanced("three_LF", nucCF, bNames[branchI], "MGL", modelList, algn_len, 1, bNames, best_models);
        }
    }
    fprintf(stdout, "\n");
}

VERBOSITY_LEVEL = 10; // 10 prints EVERYTHING

Optimize (res_three_LF,three_LF);
fprintf(stdout, "\n");

iter_likelihood = res_three_LF[1][0];
iter_parameters = res_three_LF[1][1];

iter_bic = calcBIC(iter_likelihood, iter_parameters, iter_samples);
fprintf(stdout, "This iterations likelihood: " + iter_likelihood);
fprintf(stdout, "\n");
fprintf(stdout, "This iterations parameter count: " + iter_parameters);
fprintf(stdout, "\n");
fprintf(stdout, "This iterations sample count: " + iter_samples);
fprintf(stdout, "\n");
fprintf(stdout, "This iterations BIC: " + iter_bic);
fprintf(stdout, "\n");

lfOut	= csvFilePath + ".finalTree.fit";
LIKELIHOOD_FUNCTION_OUTPUT = 7;
fprintf (lfOut, CLEAR_FILE, three_LF);
LIKELIHOOD_FUNCTION_OUTPUT = 2;

//---- TREE RENDERING -----
// Pretty self explanatory

/*
LoadFunctionLibrary ("ReadDelimitedFiles");
treeString = Format (givenTree, 1,1);
resultTable=ReadCSVTable (csvFilePath, 2);

rows	= Rows (resultTable[2]);
columns = Columns(resultTable[1]);

UseModel (USE_NO_MODEL);
Tree T = treeString;

TREE_OUTPUT_OPTIONS = {"__FONT_SIZE__":14};

tavl = T^0;
for (k = 1; k < Abs (tavl)-1; k+=1)
{
	TREE_OUTPUT_OPTIONS[(tavl[k])["Name"]] = {};
}

for (k = 1; k < Abs (tavl)-1; k+=1)
{
	thisP = (resultTable[1])[k-1][9];

	parentName = (tavl[((tavl[k])["Parent"])])["Name"];

	myRateDistro = {3,2};
	myRateDistro [0][0] = (resultTable[1])[k-1][1];
	myRateDistro [0][1] = (resultTable[1])[k-1][2];
	myRateDistro [1][0] = (resultTable[1])[k-1][3];
	myRateDistro [1][1] = (resultTable[1])[k-1][4];
	myRateDistro [2][0] = (resultTable[1])[k-1][5];
	myRateDistro [2][1] = (resultTable[1])[k-1][6];

	myRateDistro = myRateDistro % 0;

	color1 = makeDNDSColor (myRateDistro[0][0]);
	color2 = makeDNDSColor (myRateDistro[1][0]);
	color3 = makeDNDSColor (myRateDistro[2][0]);

	(TREE_OUTPUT_OPTIONS[(tavl[k])["Name"]])["TREE_OUTPUT_BRANCH_COLOR"] 		= {{color1__[0],color1__[1],color1__[2],myRateDistro__[0][1]}
																				   {color2__[0],color2__[1],color2__[2],myRateDistro__[1][1]}
																				   {color3__[0],color3__[1],color3__[2],myRateDistro__[2][1]}};
	(TREE_OUTPUT_OPTIONS[(tavl[k])["Name"]])["TREE_OUTPUT_BRANCH_LINECAP"]  = 	0;

	if (thisP <= 0.05)
	{
		(TREE_OUTPUT_OPTIONS[(tavl[k])["Name"]])["TREE_OUTPUT_BRANCH_THICKNESS"] = 	14;

	}
	if (Abs((tavl[k])["Children"]) > 0)
	{
		(TREE_OUTPUT_OPTIONS[(tavl[k])["Name"]])["TREE_OUTPUT_BRANCH_TLABEL"] = (tavl[k])["Name"]; 
	}
}

height = TipCount (T) * 36;
psTree = PSTreeString (T,"STRING_SUPPLIED_LENGTHS",{{400,height}});

treePath = csvFilePath + ".ps";

fprintf (treePath, CLEAR_FILE, psTree);
*/
return pValueByBranch;


//------------------------------------------------------------------------------------------------------------------------
function makeDNDSColor (omega)
{
	if (omega < 1)
	{
		return {{0.5*omega__,0.5*omega__,1-0.5*omega__}};
	}
	omega = Min (omega,5);
	return {{0.5+0.125*(omega__-1),0.5-(omega__-1)*0.125,0.5-(omega__-1)*0.125}};
}

//------------------------------------------------------------------------------------------------------------------------
function saveLF (ID)
{
	ExecuteCommands ("GetString(_lfInfo,"+ID+",-1)");
	_stashLF = {};
	for (_k = 0; _k < Columns (_lfInfo["Global Independent"]); _k+=1)
	{
		_stashLF [(_lfInfo["Global Independent"])[_k]] = Eval ((_lfInfo["Global Independent"])[_k]);
	}
	for (_k = 0; _k < Columns (_lfInfo["Local Independent"]); _k+=1)
	{
		_stashLF [(_lfInfo["Local Independent"])[_k]] = Eval ((_lfInfo["Local Independent"])[_k]);
	}

	return _stashLF;
}

function restoreLF (key, value)
{
	ExecuteCommands (key + " = " + value);
	return 0;
}

//------------------------------------------------------------------------------------------------------------------------
function calcBIC(iter_likelihood, iter_num_params, iter_num_samples)
{
    log_lik = iter_likelihood;
    BICtbr = -2 * log_lik;
    log_sam = Log(iter_num_samples);
    BICtbr = BICtbr + (iter_num_params * log_sam);
    return BICtbr;
}
