// XXX check initialization
// MG94 and BSREL1 don't agree. Initialize BSREL1 properly (check)
// Initialize BSREL2,3 from BSREL1, compare to BSREL1

VERBOSITY_LEVEL				 = 1;
FAST_MODE = 1;

skipCodeSelectionStep 		= 0;
LoadFunctionLibrary("chooseGeneticCode");
LoadFunctionLibrary("GrabBag");
LoadFunctionLibrary("dSdNTreeTools");
LoadFunctionLibrary("CF3x4");
LoadFunctionLibrary("BranchSiteTemplate");
LoadFunctionLibrary("AddRateClass.bf");


// MPI stuff ----------
mpi_mode = 0;
if (MPI_NODE_COUNT > 1)
{
    //fprintf(stdout, "\nWorking in mpimode...\n");
    mpi_mode = 1;
    _MPI_NODE_STATUS = {MPI_NODE_COUNT-1,1}["-1"];
}
//else
//{
    //fprintf(stdout, "\nRunning locally... \n");
//}
// ---------- MPI stuff

// This is the list of available models
// 0 = MGL
// 1 = BSREL1
// 2 = BSREL2 etc
// The list is populated as the models are created.
modelList = {};

// Input from the file
DataSet 			ds 				= ReadDataFile(PROMPT_FOR_FILE);
DataSetFilter 		dsf 			= CreateFilter(ds,3,"","",GeneticCodeExclusions);

GetInformation(dsf_seq_array, dsf);
//taxa_count = Columns(dsf_seq_array); // never used
algn_len = Abs(dsf_seq_array[0]);

HarvestFrequencies	(nuc3, dsf, 3, 1, 1);

nucCF						= CF3x4	(nuc3, GeneticCodeExclusions);

// Make the model matrix for MGL
PopulateModelMatrix			  ("MGMatrixLocal",  nucCF, "syn", "", "nonsyn");

codon3x4					= BuildCodonFrequencies (nucCF);
Model		MGL				= (MGMatrixLocal, codon3x4, 0);
modelList[0] = "MGL";

LoadFunctionLibrary			  ("queryTree");

SetDialogPrompt ("Save analysis results to");
fprintf (PROMPT_FOR_FILE, CLEAR_FILE, KEEP_OPEN,"Branch,Mean_dNdS,LRT,p,p_Holm");
csvFilePath = LAST_FILE_PATH;

fprintf 					  (stdout, "[PHASE 0] Fitting the local MG94 (no site-to-site variation) to obtain initial parameter estimates\n");

LikelihoodFunction mg94_LF  = (dsf,givenTree);

LoadFunctionLibrary			 ("DescriptiveStatistics");
totalBranchCount			 = BranchCount(givenTree) + TipCount (givenTree);
//pValueByBranch = {totalBranchCount,10}; // XXX remove after test
bNames						  = BranchName (givenTree, -1);

// Create the proper MGL parameters on the branches of the tree
/* Unless we go back to basing omega values on descriptive statistics, this
 * can be removed XXX
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
*/

//omegaStats					 = GatherDescriptiveStats (pValueByBranch[-1][0]);
//fprintf						 (stdout, "\nLog L = ", localLL, " with ", localParams, " degrees of freedom\n");


// Print the MGL .fit file
if (VERBOSITY_LEVEL >= 1)
{
    lfOut	= csvFilePath + ".MG94init.fit";
    LIKELIHOOD_FUNCTION_OUTPUT = 7;
    fprintf (lfOut, CLEAR_FILE, mg94_LF);
    LIKELIHOOD_FUNCTION_OUTPUT = 2;
}

// Optimize the MGL likelihood function
Optimize(res_mg94_LF, mg94_LF);

// Print the MGL .fit file post optimization
if (VERBOSITY_LEVEL >= 1)
{
    fprintf(stdout, "\n");
    lfOut	= csvFilePath + ".MG94opt.fit";
    LIKELIHOOD_FUNCTION_OUTPUT = 7;
    fprintf (lfOut, CLEAR_FILE, mg94_LF);
    LIKELIHOOD_FUNCTION_OUTPUT = 2;
}

// Setup variables to be filled at each iteration
iter_likelihood = 0; // Determined after Optimize
iter_samples = algn_len;
iter_parameters = 0; // Known, but can be determined after Optimize
init_parameters = 0;

GetInformation(dsf_seq_array, dsf);

// Calculate the results from the MGL fitting, to see if we can do better
origRes = res_mg94_LF[1][0] - 1.0;
orig_likelihood = res_mg94_LF[1][0];
orig_parameters = res_mg94_LF[1][1];
orig_bic = calcBIC(orig_likelihood, orig_parameters, iter_samples);

if (VERBOSITY_LEVEL >= 1)
{
    fprintf(stdout, "\nMG94 likelihood: " + orig_likelihood + " BIC: " + orig_bic + "\n\n");
}

//-------------------------------------------------------------------------
// LOCAL DONE. Starting the parallel iterative optimizer
//-------------------------------------------------------------------------

Tree						   mixtureTree = treeString;

// These don't effectively copy the t values, but they do create the t
// parameters, which is (for the time being) useful.
ReplicateConstraint 		  ("this1.?.t:=this2.?.syn",mixtureTree,givenTree);

ClearConstraints			  (mixtureTree);
ASSUME_REVERSIBLE_MODELS	  = 1;
VERBOSITY_LEVEL               = 1;

LikelihoodFunction three_LF  = (dsf,mixtureTree);

// In the MPI version of this code we need to keep track of the results from
// all branches at once, so instead of variables to store results from one
// omega to the next, we'll use arrays (ultimately of length =
// totalBranchcount)
model_list = {};
last_bics = {};
best_models = {};
done_branches = {};
working_models = {};
fixed_branches = {};
initial_ts = {};
initial_omegas = {};

/*
ts = {};
temp_omegas = {};

for (mI = 0; mI < totalBranchCount; mI = mI + 1)
{
    temp_omegas[mI] = 10;
    srate = Eval ("givenTree." + bNames[mI] + ".syn");
    nsrate = Eval ("givenTree." + bNames[mI] + ".nonsyn");
    ts[mI] = Eval ("givenTree." + bNames[mI] + ".syn");
    //fprintf(stdout, "" + nsrate + ", " + srate + "\n");
    if (srate > 0)
    {
        temp_omegas[mI] = Min (10, nsrate/srate);
    }
}
*/
//mgl_ts = {};

// Initialize these results to reasonable values
for (initI = 0; initI < totalBranchCount; initI = initI + 1)
{
    initial_ts[initI] = 1;
    initial_omegas[initI] = 0.1;
    best_models[initI] = 1; // So far the optimal omega # = 1 for all branches
    done_branches[initI] = 0; // none of the branches are done yet
    working_models[initI] = 1;  // This will be filled at each step, so these
                                // don't matter as much
    last_bics[initI] = orig_bic;// The first comparison for each additional
                                //omega will be to a version with one omega
                                //per branch, thus MGL's BIC
    fixed_branches[initI] = 0;
}
if (VERBOSITY_LEVEL > 1)
{
    fprintf(stdout, "Fixed Branches afer MG94\n");
    fprintf(stdout, fixed_branches);
    fprintf(stdout, "\n");
}
for (mI = 0; mI < totalBranchCount; mI = mI + 1)
{
    initial_omegas[mI] = 10;
    srate = Eval ("givenTree." + bNames[mI] + ".syn");
    nsrate = Eval ("givenTree." + bNames[mI] + ".nonsyn");
    initial_ts[mI] = Eval ("givenTree." + bNames[mI] + ".syn");
    //fprintf(stdout, "" + nsrate + ", " + srate + "\n");
    if (srate > 0)
    {
        initial_omegas[mI] = Min (10, nsrate/srate);
    }
}

VERBOSITY_LEVEL = 1;

//USE_MG94 = 1;
// XXX Testing interlude: is the comparison to MGL what is screwing us up?
assignModels2Branches(  "three_LF",
                        nucCF,
                        "BSREL1",
                        bNames,
                        working_models,
                        algn_len,
                        model_list,
                        fixed_branches,
                        initial_ts,
                        initial_omegas);
// Optimize the MGL likelihood function
if (VERBOSITY_LEVEL >= 1)
{
    fprintf(stdout, "\n");
    lfOut	= csvFilePath + ".BSREL1init.fit";
    LIKELIHOOD_FUNCTION_OUTPUT = 7;
    fprintf (lfOut, CLEAR_FILE, three_LF);
    LIKELIHOOD_FUNCTION_OUTPUT = 2;
}
Optimize (res_three_LF,three_LF);

// Print the MGL .fit file post optimization
if (VERBOSITY_LEVEL >= 1)
{
    fprintf(stdout, "\n");
    lfOut	= csvFilePath + ".BSREL1opt.fit";
    LIKELIHOOD_FUNCTION_OUTPUT = 7;
    fprintf (lfOut, CLEAR_FILE, three_LF);
    LIKELIHOOD_FUNCTION_OUTPUT = 2;
}

origRes = res_three_LF[1][0] - 1.0;
orig_likelihood = res_three_LF[1][0];
orig_parameters = res_three_LF[1][1];
orig_bic = calcBIC(orig_likelihood, orig_parameters, iter_samples);
if (VERBOSITY_LEVEL >= 1)
{
    fprintf(stdout, "\nBSREL1 likelihood: " + orig_likelihood + " BIC: " + orig_bic + "\n\n");
}

for (mI = 0; mI < totalBranchCount; mI = mI + 1)
{
    initial_omegas[mI] = Eval ("mixtureTree." + bNames[mI] + ".omega1");
    initial_ts[mI] = Eval ("mixtureTree." + bNames[mI] + ".t");
    if (FAST_MODE == 1)
    {
        fixed_branches[mI] = 1;
    }
}
if (VERBOSITY_LEVEL > 1)
{
    fprintf(stdout, "Fixed Branches afer BSREL1\n");
    fprintf(stdout, fixed_branches);
    fprintf(stdout, "\n");
}


branchLengths = {};
// PRINT out the calculated branch lengths
//for (bI = 0; bI < totalBranchCount; bI = bI + 1)
//{
    //calculateBranchLengthByName(modelList, best_models, "mixtureTree", bNames[bI], bI, "branchLengths");
//}
// XXX end testing interlude

/*
for (initI = 0; initI < totalBranchCount; initI = initI + 1)
{
    ExecuteCommands("mgl_ts[" + initI + "] = mixtureTree." + bNames[initI] + ".t;");
    //ExecuteCommands("fprintf(stdout, mixtureTree." + bNames[initI] + ".t);");
    //fprintf(stdout, mgl_ts[initI]);
    //fprintf(stdout, "\n");
}
*/
//if (VERBOSITY_LEVEL >= 2)
//{
    //fprintf(stdout, "\nT results from BSREL1\n");
    //fprintf(stdout, mgl_ts);
    //fprintf(stdout, "\n");
//}

branchesToOptimize = totalBranchCount;

// OPTIMIZE all branches:
while (branchesToOptimize > 0)
{
    for (launchI = 0; launchI < totalBranchCount; launchI = launchI + 1)
    {
        if (done_branches[launchI] == 0)
        {
            // working_models is used to set the model number for each
            // branch in the LF. This can be removed later on and best_models
            // can be used if you want to use the optimization results for
            // all branches to optimize for the next step.
            for (wmI = 0; wmI < totalBranchCount; wmI = wmI + 1)
            {
                working_models[wmI] = 1;
            }
            // One more omega than previously tried!
            working_models[launchI] = best_models[launchI] + 1;
            if (VERBOSITY_LEVEL >= 2)
            {
                fprintf(stdout, "\n");
                fprintf(stdout, working_models);
                fprintf(stdout, "\n");
            }
            if (FAST_MODE == 1)
            {
                fixed_branches[launchI] = 0;
                assignModels2Branches(  "three_LF",
                                        nucCF,
                                        "BSREL1",
                                        bNames,
                                        working_models,
                                        algn_len,
                                        model_list,
                                        fixed_branches,
                                        initial_ts,
                                        initial_omegas);
                fixed_branches[launchI] = 1;
            }
            else
            {
                assignModels2Branches(  "three_LF",
                                        nucCF,
                                        "BSREL1",
                                        bNames,
                                        working_models,
                                        algn_len,
                                        model_list,
                                        fixed_branches,
                                        initial_ts,
                                        initial_omegas);
            }
            // Print
            if (VERBOSITY_LEVEL >= 1)
            {
                lfOut	= csvFilePath + "." + bNames[launchI] + ".omega" + working_models[launchI] + "init.fit";
                LIKELIHOOD_FUNCTION_OUTPUT = 7;
                fprintf (lfOut, CLEAR_FILE, three_LF);
                LIKELIHOOD_FUNCTION_OUTPUT = 2;
            }
            if (mpi_mode)
            {
                // Fork:
                sendAnMPIjob(launchI);
            }
            else
            {
                // Don't fork, optimize:
                OptBranch(launchI);
                if (VERBOSITY_LEVEL >= 1)
                {
                    // Print
                    lfOut	= csvFilePath + "." + bNames[launchI] + ".omega" + working_models[launchI] + "opt.fit";
                    LIKELIHOOD_FUNCTION_OUTPUT = 7;
                    fprintf (lfOut, CLEAR_FILE, three_LF);
                    LIKELIHOOD_FUNCTION_OUTPUT = 2;
                }
            }
        }
    }
    //tempMPIcounter = branchesToOptimize;  // This is more clear than using
                                            //_MPI_NODE_STATUS
    if (mpi_mode)
    {
        // Join:
        tempMPIcounter = +(_MPI_NODE_STATUS["_MATRIX_ELEMENT_VALUE_>=0"]);
        while (tempMPIcounter > 0)
        {
            // This should get the result array back from the node, send it off
            // to get the BIC and test to see if the appropriate element of
            // best_models should be increased (or if done_branches should be
            // set to 1)
            receiveAnMPIjob();
            tempMPIcounter = tempMPIcounter -1;
        }
    }

    branchesToOptimize = 0;
    for (checkI = 0; checkI < totalBranchCount; checkI = checkI + 1)
    {
        if (done_branches[checkI] == 0)
        {
            branchesToOptimize = branchesToOptimize + 1;
        }
    }
}

// DONE: Optimal omegas are now in best_models
if (VERBOSITY_LEVEL >= 1)
{
    fprintf(stdout, "\nBest models: \n");
    fprintf(stdout, best_models);
    fprintf(stdout, "\n");
}

for (bI = 0; bI < totalBranchCount; bI = bI + 1)
{
    fixed_branches[bI] = 0;
}

// Apply best_models to a new likelihood function
assignModels2Branches(  "three_LF",
                        nucCF,
                        "BSREL1",
                        bNames,
                        best_models,
                        algn_len,
                        model_list,
                        fixed_branches,
                        initial_ts,
                        initial_omegas);

if (VERBOSITY_LEVEL >= 1)
{
    lfOut	= csvFilePath + ".finalTreeInit.fit";
    LIKELIHOOD_FUNCTION_OUTPUT = 7;
    fprintf (lfOut, CLEAR_FILE, three_LF);
    LIKELIHOOD_FUNCTION_OUTPUT = 2;
}

Optimize (res_three_LF,three_LF);
fprintf(stdout, "\n");

// Get results from the optimize
iter_likelihood = res_three_LF[1][0];
iter_parameters = res_three_LF[1][1];

iter_bic = calcBIC(iter_likelihood, iter_parameters, iter_samples);
if (VERBOSITY_LEVEL >= 1)
{
    fprintf(stdout, "This iterations likelihood: " + iter_likelihood);
    fprintf(stdout, "\n");
    fprintf(stdout, "This iterations parameter count: " + iter_parameters);
    fprintf(stdout, "\n");
    fprintf(stdout, "This iterations sample count: " + iter_samples);
    fprintf(stdout, "\n");
    fprintf(stdout, "This iterations BIC: " + iter_bic);
    fprintf(stdout, "\n");
}

if (VERBOSITY_LEVEL >= 1)
{
    lfOut	= csvFilePath + ".finalTreeOpt.fit";
    LIKELIHOOD_FUNCTION_OUTPUT = 7;
    fprintf (lfOut, CLEAR_FILE, three_LF);
    LIKELIHOOD_FUNCTION_OUTPUT = 2;
}


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
branchLengths = {};
calculateBranchLengths(modelList, best_models, "mixtureTree", bNames, "branchLengths");
// PRINT out the calculated branch lengths
//for (bI = 0; bI < totalBranchCount; bI = bI + 1)
//{
    //calculateBranchLengthByName(modelList, best_models, "mixtureTree", bNames[bI], bI, "branchLengths");
//}
//fprintf(stdout, "\nBranch lengths:\n");
//fprintf(stdout, branchLengths);
//fprintf(stdout, "\n");
ExecuteCommands ("final_tree_string = Format(" + lfTreeID + ",1,0)");
//fprintf(stdout, "\nFinal Tree before:\n");
//fprintf(stdout, final_tree_string);
//fprintf(stdout, "\n");
fprintf(stdout, "\nFinal Tree:\n");
for (bI = 0; bI < Abs(branchLengths); bI = bI + 1)
{
    final_tree_string = final_tree_string^{{bNames[bI] + ",", bNames[bI] + ":" + branchLengths[bI] + ","}};
    final_tree_string = final_tree_string^{{bNames[bI] + ")", bNames[bI] + ":" + branchLengths[bI] + ")"}};
}
fprintf(stdout, final_tree_string);
fprintf(stdout, "\n");

return pValueByBranch;

//------------------------------------------------------------------------------------------------------------------------
function sendAnMPIjob(branchNumber)
{
    // Find the next open node
    for (mpiNodeI = 0; mpiNodeI < MPI_NODE_COUNT-1; mpiNodeI += 1)
    {
        if (_MPI_NODE_STATUS[mpiNodeI] < 0)
        {
            break;
        }
    }
    // None open, wait on one to return
    if (mpiNodeI == MPI_NODE_COUNT-1)
    {
        Export(three_LF_bak, three_LF);
        mpiNodeI = receiveAnMPIjob ();
        ExecuteCommands(three_LF_bak);
    }
    fprintf(stdout, "\n sending branch number " + branchNumber + " to node number " + mpiNodeI + "plus one\n");
    _MPI_NODE_STATUS[mpiNodeI] = branchNumber;
    // XXX if you've reconstituted three_LF in a receive call before sending
    // you'll be sending an old three_LF
    MPISend(mpiNodeI + 1, three_LF);
    return 0;
}

//------------------------------------------------------------------------------------------------------------------------
function receiveAnMPIjob()
{
    // GET RESULTS
    MPIReceive (-1, fromNode, res_string);
    fromNode += (-1);
    doneID = _MPI_NODE_STATUS[fromNode];
    _MPI_NODE_STATUS[fromNode] = -1;
    // reconstitute likelihood function
    ExecuteCommands(res_string);
    // ******** HIDDEN GLOBAL VARIABLE INSERTED INTO THE NAMESPACE BY HYPHY ********
    res_three_LF = three_LF_MLES;   // it will be your likelihood function
                                    // name followed by _MLES
    // ******** HIDDEN GLOBAL VARIABLE INSERTED INTO THE NAMESPACE BY HYPHY ********
    fprintf(stdout, "\nJob received from " + fromNode + " saving to branch number " + doneID + "\n");
    if (VERBOSITY_LEVEL >= 1)
    {
        // Print
        lfOut	= csvFilePath + "." + bNames[doneID] + ".omega" + (best_models[doneID] + 1) + "opt.fit";
        LIKELIHOOD_FUNCTION_OUTPUT = 7;
        fprintf (lfOut, CLEAR_FILE, three_LF);
        LIKELIHOOD_FUNCTION_OUTPUT = 2;
    }
    //fprintf(stdout, res_three_LF);
    //fprintf(stdout, "\n");
    processResults(res_three_LF, doneID);
    return fromNode;
}

//------------------------------------------------------------------------------------------------------------------------
function OptBranch(nodeI)
{
    Optimize(res_three_LF, three_LF);
    processResults(res_three_LF, nodeI);
    return 0;
}

//------------------------------------------------------------------------------------------------------------------------
function processResults(res_LF, nodeID)
{
    thisRes = res_LF[1][0];
    this_likelihood = res_LF[1][0];
    this_parameters = res_LF[1][1];
    if (FAST_MODE == 1)
    {
        this_parameters = orig_parameters + best_models[nodeID] + 1;
    }
    this_bic = calcBIC(this_likelihood, this_parameters, algn_len); // algn_len is a global variable
    if (VERBOSITY_LEVEL >= 1)
    {
        fprintf(stdout, "\nThis likelihood: " + this_likelihood + " parameter count: " + this_parameters + "\n");
        fprintf(stdout, "This BIC: " + this_bic + " Last bic: " + last_bics[nodeID] + "\n\n");
    }
    if (last_bics[nodeID] <= this_bic) // last_bics is a global variable
    {
        done_branches[nodeID] = 1; // done_branches is a global variable
    }
    else
    {
        last_bics[nodeID] = this_bic;
        best_models[nodeID] = best_models[nodeID] + 1; // best_models is a global variable
    }
    // replaced in the main optimizer loop
    //if (VERBOSITY_LEVEL >= 5)
    //{
        //lfOut	= csvFilePath + ".branch" + nodeID + ".omega" + best_models[nodeID] + ".fit";
        //LIKELIHOOD_FUNCTION_OUTPUT = 7;
        //fprintf (lfOut, CLEAR_FILE, three_LF);
        //LIKELIHOOD_FUNCTION_OUTPUT = 2;
    //}
    return thisRes;
}

//------------------------------------------------------------------------------------------------------------------------
function calculateBranchLengths(modelList, bestModels, treeName, branchNames, branchLengthsID)
{
    for (bI = 0; bI < Abs(bestModels); bI = bI + 1)
    {
        evalstring = "(";
        branchLength = 0;
        ExecuteCommands("t = `treeName`." + branchNames[bI] + ".t;");
        for (mI = 1; mI <= bestModels[bI]; mI = mI + 1)
        {
            // values
            ExecuteCommands("Model M" + mI + "= (MGMatrix" + mI + ", codon3x4, 0);");
            ExecuteCommands("GetString(M" + mI + "L, M" + mI + ", -1);");
            ExecuteCommands("omega" + mI + " = `treeName`."+ branchNames[bI] + ".omega" + mI + ";");
            if (mI != bestModels[bI])
            {
                ExecuteCommands("Paux" + mI + " = `treeName`."+ branchNames[bI] + ".Paux" + mI + ";");
            }

            // OLD
            if (mI > 1)
            {
                evalstring = evalstring + "+";
            }
            //ExecuteCommands("tempML = " + mls[oI] + ";");
            ExecuteCommands("tempML = M" + mI + "L;");
            ExecuteCommands("evalstring = evalstring + \"(" + tempML + ")\";");
            if (mI != bestModels[branchNumber])
            {
                evalstring = evalstring + "*Paux" + mI;
            }
            for (prevMI = mI - 1; prevMI > 0; prevMI = prevMI - 1)
            {
                evalstring = evalstring + "*(1-Paux" + prevMI + ")";
            }

            //omegas[oI] = Eval("`treeName`.`branchName`.omega" + mI);
            //if (oI < bestModels[branchNumber])
            //{
                //pauxs[oI] = Eval("`treeName`.`branchName`.Paux" + oI);
            //}
            // You divide the estimated length by three because canonical branch
            // length calculations use the nucleotide as the unit of evolution. We're
            // using the codon here.
            // END OLD

            // string
            //ExecuteCommands("tempML = M" + mI + "L;");
            //ExecuteCommands("evalString = evalString + (`tempML`)");
            //if (mI != bestModels[bI])
            //{
                //evalstring = evalstring + "*Paux" + mI;
            //}
            //for (prevMI = mI - 1; prevMI > 0; prevMI = prevMI - 1)
            //{
                //evalstring = evalstring + "*(1-Paux" + prevMI + ")";
            //}
        }
        evalstring = evalstring + ")/3";
        //evalstring = evalstring + ")/3";
        //fprintf(stdout, "Evalstring: \n");
        //fprintf(stdout, evalstring);
        //fprintf(stdout, "\n");
        branchLength = Eval(evalstring);
        if (VERBOSITY_LEVEL >= 1)
        {
            fprintf(stdout, "\n");
            fprintf(stdout, branchNames[bI]);
            fprintf(stdout, ":\n");
            fprintf(stdout, "\n");
            fprintf(stdout, evalstring);
            fprintf(stdout, "\n");
            fprintf(stdout, "\n");
            fprintf(stdout, branchLength);
            fprintf(stdout, "\n");
        }
        ExecuteCommands(branchLengthsID + "[" + bI + "] = " + branchLength + ";");
    }

    return 0;
}

//------------------------------------------------------------------------------------------------------------------------
function calculateBranchLengthByName (modelList, bestModels, treeName, branchName, branchNumber, branchLengthsID)
{
    branchLength = 0;

    model_matrices = {};
    for (mI = 1; mI <= bestModels[branchNumber]; mI = mI + 1)
    {
        model_matrices[mI] = "MGMatrix" + mI;
    }

    models = {};
    for (mI = 1; mI <= bestModels[branchNumber]; mI = mI + 1)
    {
        ExecuteCommands ("Model M" + mI + " = (MGMatrix" + mI + ", codon3x4, 0)");
        models[mI] = "M" + mI;
    }

    mls = {};
    for (mI = 1; mI <= bestModels[branchNumber]; mI = mI + 1)
    {
        ExecuteCommands ("GetString (M" + mI + "L, M" + mI + ", -1);");
        mls[mI] = "M" + mI + "L";
    }

    t = Eval("`treeName`.`branchName`.t");
    omegas = {};
    pauxs = {};
    evalstring = "(";
    for (oI = 1; oI <= bestModels[branchNumber]; oI = oI + 1)
    {
        if (oI > 1)
        {
            evalstring = evalstring + "+";
        }
        ExecuteCommands("tempML = " + mls[oI] + ";");
        ExecuteCommands("evalstring = evalstring + \"(" + tempML + ")\";");
        if (oI != bestModels[branchNumber])
        {
            evalstring = evalstring + "*Paux" + oI;
        }
        for (prevOI = oI - 1; prevOI > 0; prevOI = prevOI - 1)
        {
            evalstring = evalstring + "*(1-Paux" + prevOI + ")";
        }

        omegas[oI] = Eval("`treeName`.`branchName`.omega" + oI);
        if (oI < bestModels[branchNumber])
        {
            pauxs[oI] = Eval("`treeName`.`branchName`.Paux" + oI);
        }
    }
    // You divide the estimated length by three because canonical branch
    // length calculations use the nucleotide as the unit of evolution. We're
    // using the codon here.
    evalstring = evalstring + ")/3";
    branchLength = Eval(evalstring);
    if (VERBOSITY_LEVEL >= 1)
    {
        fprintf(stdout, "\n");
        fprintf(stdout, branchName);
        fprintf(stdout, ":\n");
        fprintf(stdout, "\n");
        fprintf(stdout, evalstring);
        fprintf(stdout, "\n");
        fprintf(stdout, "\n");
        fprintf(stdout, branchLength);
        fprintf(stdout, "\n");
    }
    ExecuteCommands(branchLengthsID + "[branchNumber] = " + branchLength + ";");

    return branchLength;
}

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
