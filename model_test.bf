VERBOSITY_LEVEL				 = 0;

skipCodeSelectionStep 		= 0;
LoadFunctionLibrary("chooseGeneticCode");   // Where are these libraries (for
                                            //reference/use)

LoadFunctionLibrary("GrabBag");
LoadFunctionLibrary("dSdNTreeTools");
LoadFunctionLibrary("CF3x4");
LoadFunctionLibrary("BranchSiteTemplate");
//LoadFunctionLibrary("DeltaOmega.bf"); XXX no longer needed
LoadFunctionLibrary("AddRateClass.bf");


modelList = {};


DataSet 			ds 				= ReadDataFile(PROMPT_FOR_FILE);
DataSetFilter 		dsf 			= CreateFilter(ds,3,"","",GeneticCodeExclusions);
//GetInformation(infotainer, dsf); XXX no longer needed
//fprintf(stdout, "The original model!"); XXX no longer needed
//fprintf(stdout, infotainer); XXX no longer needed

// Are these the nucleotide frequences?
// If so, are they harvested but left in dsf? Or are they now in nuc3?
HarvestFrequencies	(nuc3, dsf, 3, 1, 1);


// These are the rate class priors for the global fitting
//omega1 = 0.2; XXX no longer needed
//omega2 = 0.5; XXX no longer needed
//omega3 = 1.0; XXX no longer needed

// I don't know what this is yet, but nuc3 is passed in as an argument here.
nucCF						= CF3x4	(nuc3, GeneticCodeExclusions);

// Why is this called multiple times? Should these function names be
// imperative hints at their side effects?
//PopulateModelMatrix			  ("MGMatrix1",  nucCF, "t", "omega1", "");
//PopulateModelMatrix			  ("MGMatrix2",  nucCF, "t", "omega2", "");
//PopulateModelMatrix			  ("MGMatrix3",  nucCF, "t", "omega3", "");

// What is the difference between omega1 and omegaG1?
// What is the significance of declaring this global? I thought all variables
// were in the global namespace?
//global	omegaG1 = 0.2;
// I'm assuming this is a bound such that omegaG1 cannot exceed 1?
//omegaG1 :< 1;
//global	omegaG2 = 0.5;
//omegaG2 :< 1;
//global	omegaG3 = 2.0;
//omegaG3 :> 1;

// So these are the same function calls as earlier but with different inputs.
// Everything has a G in it, for instance
//PopulateModelMatrix			  ("MGMatrix1G",  nucCF, "t", "omegaG1", "");
//PopulateModelMatrix			  ("MGMatrix2G",  nucCF, "t", "omegaG2", "");
//PopulateModelMatrix			  ("MGMatrix3G",  nucCF, "t", "omegaG3", "");


// Why is this one separate? Matrix Local with syn and nonsyn? hmmm
PopulateModelMatrix			  ("MGMatrixLocal",  nucCF, "syn", "", "nonsyn");

// So we harvested frequencies earlier (character frequencies), now we're
// building codon frequencies from the nucleotide frequencies?
codon3x4					= BuildCodonFrequencies (nucCF);
// I'm not sure I get this syntax. Is this some sort of constructor, or
// tuple?
Model		MGL				= (MGMatrixLocal, codon3x4, 0);

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

// We are now translating from the global fitting to the local prior.

// Awesome, so here we are using the omega stats from the global fitting as
// a starting point for the local model fitting
//omegaG1						 = omegaStats["2.5%"];
//omegaG2						 = omegaStats["Median"];
//omegaG3						 = omegaStats["97.5%"];

// oh yay, more variables named aux
//Paux1 						 = 0.3;
//Paux1 						 :< 1;
//Paux2 						 = 0.4;
//Paux2 						 :< 1;
//
//global Paux1G 				  = 0.3;
//global Paux2G 				  = 0.4;


// Because this totally makes sense. Also, is treeString another magical
// global variable?
// I'm not sure if MGG is a transition matrix or a model in the traditional
// sense... but regardless we are generating a model and a tree using our
// Global priors.
// Actually it does in a way, because you are exponentiating the
// instantaneous transition matrix and multiplying it by the proportion of
// the overall model that that sub model has been estimated to contribute.
//Model 		MGG		=		  ("Exp(MGMatrix1G)*Paux1G+Exp(MGMatrix2G)*(1-Paux1G)*Paux2G+Exp(MGMatrix3G)*(1-Paux1G)*(1-Paux2G)",codon3x4,EXPLICIT_FORM_MATRIX_EXPONENTIAL);
//Tree						   mixtureTreeG = treeString;
//
//Model 		MG1		=		  ("Exp(MGMatrix1)*Paux1+Exp(MGMatrix2)*(1-Paux1)*Paux2+Exp(MGMatrix3)*(1-Paux1)*(1-Paux2)",codon3x4,EXPLICIT_FORM_MATRIX_EXPONENTIAL);
Tree						   mixtureTree = treeString;



// Say what? are we saving parameters into temptorary variables or something?
//ReplicateConstraint 		  ("this1.?.t:=this2.?.syn",mixtureTreeG,givenTree);
ReplicateConstraint 		  ("this1.?.t:=this2.?.syn",mixtureTree,givenTree);

// And then clearing them...
ClearConstraints			  (mixtureTree);
//ClearConstraints			  (mixtureTreeG);

// Why are these bound down here? Are they affected by the side effects of
// one of the above functions?
//omegaG1						 :< 1;
//omegaG2						 :< 1;
//Paux1G 						 :< 1;
//Paux2G 						 :< 1;

// Some options
ASSUME_REVERSIBLE_MODELS	  = 1;

VERBOSITY_LEVEL               = 1;


//LikelihoodFunction global_LF  = (dsf, mixtureTreG);

// Yay new LF
// So this is the likelihood function. But there is no "Optimize()" call...
// yet. There wasn't earlier either, but there is one later that takes this
// as the second parameter (XXX look into the Optimize function)
LikelihoodFunction three_LF   = (dsf,mixtureTree);

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


Export(three_LF_bak, three_LF);

origRes = res_three_LF[1][0] - 1.0;

for (branchI = 0; branchI < totalBranchCount; branchI = branchI + 1)
{
    lastRes = origRes;
    omegaNumber = 1;
    //while (res_three_LF[1][0] > lastRes)


    for (omegaNumber = 1; omegaNumber < 3; omegaNumber = omegaNumber + 1)
    {
        fprintf(stdout, "\n");
        addRate2Branch("three_LF", nucCF, bNames[branchI], "MGL", modelList);
        fprintf(stdout, "Current Model list:\n");
        fprintf(stdout, modelList);
        fprintf(stdout, "\n");
        //omegaNumber = omegaNumber + 1;

        lfOut	= csvFilePath + ".treePlusRate." + bNames[branchI] + "." + omegaNumber + ".fit";
        LIKELIHOOD_FUNCTION_OUTPUT = 7;
        fprintf (lfOut, CLEAR_FILE, three_LF);
        LIKELIHOOD_FUNCTION_OUTPUT = 2;

        VERBOSITY_LEVEL = 1;
        //VERBOSITY_LEVEL = 0;

        Optimize (res_three_LF,three_LF);
        fprintf(stdout, "\n");

        lfOut	= csvFilePath + ".optTreePlusRate." + bNames[branchI] + "." + omegaNumber + ".fit";
        LIKELIHOOD_FUNCTION_OUTPUT = 7;
        fprintf (lfOut, CLEAR_FILE, three_LF);
        LIKELIHOOD_FUNCTION_OUTPUT = 2;

        lastRes = res_three_LF[1][0];
    }
}

/*
for (k = 0; k < totalBranchCount; k = k+1)
{
    // I'd like to be able to write out a new tree at every step so that I
    // can inspect the parameters. This requires that the trees are written
    // out in the DeltaOmega.bf, which I haven't implemented yet. Uncomment
    // this when it happens. XXX
    lfOut	= csvFilePath + ".optTree." + k + ".fit";
    LIKELIHOOD_FUNCTION_OUTPUT = 7;
    fprintf (lfOut, CLEAR_FILE, three_LF);
    LIKELIHOOD_FUNCTION_OUTPUT = 2;
    lfOut	= csvFilePath + ".optTree.0.fit";
    inputOptions = {};
    inputOptions["0"] = lfOut;
    inputOptions["1"] = "" + k;
    ExecuteAFile("DeltaOmega.bf", inputOptions);

    // Go go gadget initialization routines!
    if (k == 0)
    {
        expr            = Eval("BranchLength(givenTree,\""+bNames[0]+";EXPECTED_NUMBER_OF_SUBSTITUTIONS\")");
        syn             = 1; nonsyn = 0;
        synM            = Eval(expr);
        syn             = 0; nonsyn = 1;
        nonsynM         = Eval(expr);
    }

 	srate  = Eval ("givenTree." + bNames[k] + ".syn");
	nsrate = Eval ("givenTree." + bNames[k] + ".nonsyn");
    bl = Eval("BranchLength(givenTree,\""+bNames[k]+"\")")*3;

    if (srate > 0)
    {
        baseOmega = nsrate/srate;
    }
    else
    {
        baseOmega = 10000;
    }

    bl = bl / (synM + nonsynM * baseOmega);

    // These look like the same commands that I'm seeing in the rest of this
    // HBL file. Why are they in the ExecuteCommands framework?
    // so here we go. This is the magic. We are adding the omega variables
    // described above to the branches. Awesome. So this is what I need to
    // replicate in the external file.
    ExecuteCommands ("mixtureTree." + bNames[k] + ".t = bl");
    ExecuteCommands ("mixtureTree." + bNames[k] + ".omega1 :< 1;");
	ExecuteCommands ("mixtureTree." + bNames[k] + ".omega2 :< 1;");
    if (baseOmega > 1)
    {
        ExecuteCommands ("mixtureTree." + bNames[k] + ".omega1 = 0.1;");
        ExecuteCommands ("mixtureTree." + bNames[k] + ".omega2 = 1;");
        ExecuteCommands ("mixtureTree." + bNames[k] + ".omega3 = baseOmega;");

        ExecuteCommands ("mixtureTree." + bNames[k] + ".Paux1 = 0.01;");
        ExecuteCommands ("mixtureTree." + bNames[k] + ".Paux2 = 0.01;");
    }
    else
    {
        ExecuteCommands ("mixtureTree." + bNames[k] + ".omega1 = baseOmega;");
        ExecuteCommands ("mixtureTree." + bNames[k] + ".omega2 = 1;");
        ExecuteCommands ("mixtureTree." + bNames[k] + ".omega3 = 2;");

        ExecuteCommands ("mixtureTree." + bNames[k] + ".Paux1 = 0.98;");
        ExecuteCommands ("mixtureTree." + bNames[k] + ".Paux2 = 0.5;");
    }
}


VERBOSITY_LEVEL = 1;
USE_LAST_RESULTS = 1;
OPTIMIZATION_METHOD = 0;

// NOW DO THE ACTUAL LOCAL FITTING!!!! I don't see evidence of this really
// being done with the global model. How come?
fprintf 					  (stdout, "[PHASE 2] Fitting the full LOCAL alternative model (no constraints)\n");
Optimize					  (res_three_LF,three_LF);
fprintf						  (stdout,"\n",three_LF);



lfOut	= csvFilePath + ".fit";
// So LIKELIHOOD_FUNCTION_OUTPUT is a hidden parameter to something?
LIKELIHOOD_FUNCTION_OUTPUT = 7;
fprintf (lfOut, CLEAR_FILE, three_LF);
LIKELIHOOD_FUNCTION_OUTPUT = 2;



// So at this point the optimization stuff is done.
for	(k = 0; k < totalBranchCount; k = k+1)
{
	ref = "mixtureTree."+bNames[k];



    //Dataset ds = ReadDataFile(lfOut);
    //DataSetFilter dsf = CreateFilter(ds,3,"","",GeneticCodeExclusions);
    //HarvestFrequencies (nuc3, dsf, 3, 1, 1);
    //nucCF = CF3x4(nuc3, GeneticCodeExclusions);
    //LikelihoodFunction new_LF = (dsf, mixtureTree);

    //lfOut	= csvFilePath + ".export.fit";
    //Export(lfOut, three_LF);


    //optimize_branch (lfOut);
    //changeLF (lfOut);
    //LIKELIHOOD_FUNCTION_OUTPUT = 7;
    //fprintf (lfOut, CLEAR_FILE, three_LF);
    //LIKELIHOOD_FUNCTION_OUTPUT = 2;
    //_stashLF ["restoreLF"][""];
    //lfOut	= csvFilePath + ".restore.after.fit";
    //LIKELIHOOD_FUNCTION_OUTPUT = 7;
    //fprintf (lfOut, CLEAR_FILE, three_LF);
    //LIKELIHOOD_FUNCTION_OUTPUT = 2;

	thisOmega3 = Eval (ref+".omega3");
	wt3        = Eval ("(1-"+ref+".Paux1)*(1-"+ref+".Paux2)");

	pValueByBranch [k][1] = Eval (ref+".omega1");
	pValueByBranch [k][2] = Eval (ref+".Paux1");
	pValueByBranch [k][3] = Eval (ref+".omega2");
	pValueByBranch [k][4] = Eval ("(1-"+ref+".Paux1)*"+ref+".Paux2");
	pValueByBranch [k][5] = thisOmega3;
	pValueByBranch [k][6] = wt3;

	fprintf (stdout, "\nNode: ", ref,
					 "\n\tClass 1: omega = ", Eval (ref+".omega1"), " weight = ", Eval (ref+".Paux1"),
					 "\n\tClass 2: omega = ", Eval (ref+".omega2"), " weight = ", Eval ("(1-"+ref+".Paux1)*"+ref+".Paux2"),
					 "\n\tClass 3: omega = ", thisOmega3, " weight = ", wt3 , "\n"
					 );

	if (thisOmega3 > 1 && wt3 > 1e-6)
	{
		fprintf (stdout, "...Testing for selection at this branch\n");
        //lfOut	= csvFilePath + ".lf.before." + k + ".fit";
        //LIKELIHOOD_FUNCTION_OUTPUT = 7;
        //fprintf (lfOut, CLEAR_FILE, three_LF);
        //LIKELIHOOD_FUNCTION_OUTPUT = 2;
		_stashLF = saveLF ("three_LF");
        //lfOut	= csvFilePath + ".lf.after." + k + ".fit";
        //LIKELIHOOD_FUNCTION_OUTPUT = 7;
        //fprintf (lfOut, CLEAR_FILE, three_LF);
        //LIKELIHOOD_FUNCTION_OUTPUT = 2;

        // Another round of optimization to test for selection at this branch
		ExecuteCommands ("mixtureTree." + bNames[k] + ".omega3 := 1");
		Optimize					  (res_three_current,three_LF);

        // Ahah, the p-aux variables are weight related
		fprintf (stdout, "\nNode: ", ref,
						 "\n\tClass 1: omega = ", Eval (ref+".omega1"), " weight = ", Eval (ref+".Paux1"),
						 "\n\tClass 2: omega = ", Eval (ref+".omega2"), " weight = ", Eval ("(1-"+ref+".Paux1)*"+ref+".Paux2"),
						 "\n\tClass 3: omega = ", Eval (ref+".omega3"), " weight = ", Eval ("(1-"+ref+".Paux1)*(1-"+ref+".Paux2)"), "\n"
						 );
		pValueByBranch[k][7]			  = 2*(res_three_LF[1][0] - res_three_current[1][0]);
		pValueByBranch[k][8]			  = (1-CChi2 (pValueByBranch[k][7],1))*.5;
		fprintf (stdout, "\np-value = ", pValueByBranch[k][8],"\n\n", three_LF, "\n");

		ExecuteCommands ("mixtureTree." + bNames[k] + ".omega3 :< 1e26");

		if (pValueByBranch[k][7] < (-0.5))
		{
			fprintf 					  (stdout, "[PHASE 2/REPEAT] Detected a convergence problem; refitting the LOCAL alternative model with new starting values\n");
			lfOut	= csvFilePath + ".fit";
            // And another branch, but only if in an error state
			Optimize					  (res_three_LF,three_LF);
			LIKELIHOOD_FUNCTION_OUTPUT = 7;
			fprintf (lfOut, CLEAR_FILE, three_LF);
			LIKELIHOOD_FUNCTION_OUTPUT = 2;
			_stashLF = saveLF ("three_LF");
			k = 0;
		}
		else
		{
			//lfOut	= csvFilePath + ".restore.before.fit";
            //fprintf (stdout, "in this thing because the p-value is > -.5\n");
			//LIKELIHOOD_FUNCTION_OUTPUT = 7;
            //fprintf (lfOut, CLEAR_FILE, three_LF);
			//LIKELIHOOD_FUNCTION_OUTPUT = 2;
			_stashLF ["restoreLF"][""];
			//lfOut	= csvFilePath + ".restore.after.fit";
			//LIKELIHOOD_FUNCTION_OUTPUT = 7;
            //fprintf (lfOut, CLEAR_FILE, three_LF);
			//LIKELIHOOD_FUNCTION_OUTPUT = 2;
		}
	}
	else
	{
		pValueByBranch[k][8] = 1.0;
	}
}

// I'm assuming this is used as a hidden parameter to a later function
OPTIMIZATION_METHOD = 4;


pValueSorter = {totalBranchCount,2};
pValueSorter = pValueSorter["_MATRIX_ELEMENT_ROW_*(_MATRIX_ELEMENT_COLUMN_==0)+pValueByBranch[_MATRIX_ELEMENT_ROW_][8]*(_MATRIX_ELEMENT_COLUMN_==1)"];
pValueSorter = pValueSorter % 1;
pValueSorter = pValueSorter["_MATRIX_ELEMENT_VALUE_*(_MATRIX_ELEMENT_COLUMN_==0)+_MATRIX_ELEMENT_VALUE_*(totalBranchCount-_MATRIX_ELEMENT_ROW_)*(_MATRIX_ELEMENT_COLUMN_==1)"];

fprintf (stdout,"\n\nSummary of branches under episodic selection:\n");
hasBranchesUnderSelection = 0;

pthreshold = 0.05;

for		(k = 0; k < totalBranchCount; k = k+1)
{
	pValueByBranch[pValueSorter[k][0]][9] = Min (1,pValueSorter[k][1]);
	if (pValueSorter[k][1] <= pthreshold)
	{
		fprintf (stdout, "\t", bNames[pValueSorter[k][0]], " p = ", pValueByBranch[pValueSorter[k][0]][9], "\n");
		hasBranchesUnderSelection += 1;
	}
}


if (hasBranchesUnderSelection == 0)
{
	fprintf (stdout, "\tNo branches found to be under selection at p <= ", threshold, "\n");
}


for		(k = 0; k < totalBranchCount; k = k+1)
{
	fprintf (csvFilePath, "\n", bNames[k], ",", Join(",",pValueByBranch[k][-1]));
}
fprintf (csvFilePath, CLOSE_FILE);
*/

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

