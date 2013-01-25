VERBOSITY_LEVEL				 = 0;

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
    fprintf(stdout, "\nWorking in mpimode...\n");
    mpi_mode = 1;
    _MPI_NODE_STATUS = {MPI_NODE_COUNT-1,1}["-1"];
}
else
{
    fprintf(stdout, "\nRunning locally... \n");
}
// ---------- MPI stuff

modelList = {};

DataSet 			ds 				= ReadDataFile(PROMPT_FOR_FILE);
DataSetFilter 		dsf 			= CreateFilter(ds,3,"","",GeneticCodeExclusions);

GetInformation(dsf_seq_array, dsf);
taxa_count = Columns(dsf_seq_array);
algn_len = Abs(dsf_seq_array[0]);

HarvestFrequencies	(nuc3, dsf, 3, 1, 1);

nucCF						= CF3x4	(nuc3, GeneticCodeExclusions);

PopulateModelMatrix			  ("MGMatrixLocal",  nucCF, "syn", "", "nonsyn");

codon3x4					= BuildCodonFrequencies (nucCF);
Model		MGL				= (MGMatrixLocal, codon3x4, 0);
//Model		BSREL1				= (MGMatrixLocal, codon3x4, 0);
//modelList[1] = "MGL";
// Testing
modelList[0] = "MGL";

LoadFunctionLibrary			  ("queryTree");

SetDialogPrompt ("Save analysis results to");
fprintf (PROMPT_FOR_FILE, CLEAR_FILE, KEEP_OPEN,"Branch,Mean_dNdS,LRT,p,p_Holm");
csvFilePath = LAST_FILE_PATH;

fprintf 					  (stdout, "[PHASE 0] Fitting the local MG94 (no site-to-site variation) to obtain initial parameter estimates\n");

LikelihoodFunction	base_LF	 = (dsf, givenTree);

LoadFunctionLibrary			 ("DescriptiveStatistics");
totalBranchCount			 = BranchCount(givenTree) + TipCount (givenTree);
pValueByBranch				  = {totalBranchCount,10};
bNames						  = BranchName (givenTree, -1);

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

omegaStats					 = GatherDescriptiveStats (pValueByBranch[-1][0]);
fprintf						 (stdout, "\nLog L = ", localLL, " with ", localParams, " degrees of freedom\n");

Tree						   mixtureTree = treeString;

ReplicateConstraint 		  ("this1.?.t:=this2.?.syn",mixtureTree,givenTree);

ClearConstraints			  (mixtureTree);
ASSUME_REVERSIBLE_MODELS	  = 1;
VERBOSITY_LEVEL               = 1;

LikelihoodFunction three_LF  = (dsf,mixtureTree);

    /*
    temp_omegas = {};
    for (mI = 0; mI < Abs(model_assignments); mI = mI + 1)
    {
        //temp_omegas = {};
        //if (model_assignments[mI] == 1)
        //{
            temp_omegas[mI] = 10;
            srate = Eval (lfTreeID + "." + branchName + ".syn");
            nsrate = Eval (lfTreeID + "." + branchName + ".nonsyn");
            fprintf(stdout, "" + nsrate + ", " + srate + "\n");
            if (srate > 0)
            {
                temp_omegas[mI] = Min (10, nsrate/srate);
            }
        //}
        //else
        //{
            //for (oI = 1; oI < model_assignments[mI]; oI = oI + 1)
            //{
                //ExecuteCommands( "tempOmega = " + lfTreeID + "." + branch_names[mI] + ".omega" + oI + ";");
                //temp_omegas[oI] = tempOmega;
            //}
        //}
        //temp_branches[mI] = tempOmega;
    }
    */
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


iter_likelihood = 0; // Determined after Optimize
iter_samples = algn_len;
iter_parameters = 0; // Known, but can be determined after Optimize
init_parameters = 0;

GetInformation(dsf_seq_array, dsf);

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

for (branchI = 0; branchI < totalBranchCount; branchI = branchI + 1)
{
    if (mpi_mode)
    {
        fprintf(stdout, "\nSending an MPI job\n");
        sendAnMPIjob ("three_LF", nucCF, bNames[branchI], branchI, "MGL", modelList, algn_len, 1, origRes, orig_bic);
    }
    else
    {
        best_models[branchI] = optimizeBranchOmegas("three_LF", nucCF, bNames[branchI], branchI, "MGL", modelList, algn_len, 1, origRes, orig_bic)
    /*
        fprintf(stdout, "\n");
        fprintf(stdout, "New branch: \n");
        lastRes = origRes;
        omegaNumber = 1;

        better_bic = 1;
        last_bic = orig_bic;

        bic_run_count = 0;

        while (better_bic == 1)
        {
            addRate2Branch("three_LF", nucCF, bNames[branchI], "MGL", modelList, algn_len, 1);
            //addRate2Branch("three_LF", nucCF, bNames[branchI], "BSREL1", modelList, algn_len, 1);
            omegaNumber = omegaNumber + 1;

            lfOut	= csvFilePath + ".treePlusRate." + bNames[branchI] + "." + omegaNumber + ".fit";
            LIKELIHOOD_FUNCTION_OUTPUT = 7;
            fprintf (lfOut, CLEAR_FILE, three_LF);
            LIKELIHOOD_FUNCTION_OUTPUT = 2;

            VERBOSITY_LEVEL = 10;   // 10 prints EVERYTHING
                                    // 0 prints nothing

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
        */
    } // endif (mpi_mode)
}
if (mpi_mode)
{
    tempMPIcounter = totalBranchCount;
    tempMPIcounter = +(_MPI_NODE_STATUS["_MATRIX_ELEMENT_VALUE_>=0"]);
    while (tempMPIcounter > 0)
    {
        receiveAnMPIjob ();
        fprintf(stdout, "\nMPI job received\n");
        //best_models[results[0]] = results[1];
        tempMPIcounter = tempMPIcounter -1;
    }
}



fprintf(stdout, "\n");
fprintf(stdout, best_models);
fprintf(stdout, "\n");

//assignModels2Branches("three_LF", nucCF, "MGL", bNames, best_models, algn_len, model_list);
// Testing
/*
for (bI = 0; bI < totalBranchCount; bI = bI + 1)
{
    best_models[bI] = 1;
}
*/
// XXX If the optimal number of omegas is discovered using MPI, BSREL1 might
// not exist, in fact the entire model list may be empty on the host node.
assignModels2Branches("three_LF", nucCF, "BSREL1", bNames, best_models, algn_len, model_list);

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
for (bI = 0; bI < totalBranchCount; bI = bI + 1)
{
    calculateBranchLengthByName(modelList, best_models, "mixtureTree", bNames[bI], bI);
}
//calculateBranchLengthByName(modelList, best_models, "mixtureTree", bNames[0], 0);
//fprintf(stdout, extractBranchInformation ("mixtureTree", "omega", "MGMatrix", "Paux", "codon3x4", 0));

return pValueByBranch;

//------------------------------------------------------------------------------------------------------------------------
function optimizeBranchOmegas(lfID, nucCF, branchName, branchNumber, defaultModel, modelList, algn_len, replace_tree, origRes, orig_bic)
{
    tempOmegaNumber = 1;
    lastRes = origRes;
    omegaNumber = 1;
    better_bic = 1;
    last_bic = orig_bic;

    bic_run_count = 0;

    while (better_bic == 1)
    {
        addRate2Branch("three_LF", nucCF, branchName, "MGL", modelList, algn_len, 1);
        omegaNumber = omegaNumber + 1;

        /*
        lfOut	= csvFilePath + ".treePlusRate." + bNames[branchI] + "." + omegaNumber + ".fit";
        LIKELIHOOD_FUNCTION_OUTPUT = 7;
        fprintf (lfOut, CLEAR_FILE, three_LF);
        LIKELIHOOD_FUNCTION_OUTPUT = 2;

        VERBOSITY_LEVEL = 10;   // 10 prints EVERYTHING
                                // 0 prints nothing
        */

        Optimize (res_three_LF,three_LF);
        //fprintf(stdout, "\n");

        iter_likelihood = res_three_LF[1][0];
        iter_parameters = res_three_LF[1][1];

        iter_bic = calcBIC(iter_likelihood, iter_parameters, iter_samples);
        /*
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
        */

        lastRes = res_three_LF[1][0];
        if (iter_bic == last_bic)
        {
            bic_run_count = bic_run_count + 1;
        }
        if ((iter_bic > last_bic) || (bic_run_count > 2))
        {
            //fprintf(stdout, "Done with this branch...\n");
            better_bic = 0;
        }
        else
        {
            tempOmegaNumber = tempOmegaNumber + 1;
        }
        last_bic = iter_bic;
    }
    // fprintf(stdout, "\n");
    return tempOmegaNumber;
}

//------------------------------------------------------------------------------------------------------------------------
function sendAnMPIjob (lfID, nucCF, branchName, branchNumber, defaultModel, modelList, algn_len, replace_tree, origRes, orig_bic)
{
    fprintf(stdout, "\nLooking for MPI node...\n");
    for (mpiNodeI = 0; mpiNodeI < MPI_NODE_COUNT-1; mpiNodeI += 1)
    {
        if (_MPI_NODE_STATUS[mpiNodeI] < 0)
        {
            break;
        }
    }
    fprintf(stdout, "\nMPI node found\n");
    if (mpiNodeI == MPI_NODE_COUNT-1)
    {
        mpiNodeI = receiveAnMPIjob ();
    }
    _MPI_NODE_STATUS[mpiNodeI] = branchNumber;
    // Reference best_models[branchI] = optimizeBranchOmegas("three_LF", nucCF, bNames[branchI], branchI, "MGL", modelList, algn_len, 1, origRes, orig_bic)
    mpiCommandString = "optimizeBranchOmegas(\"`lfID`\", nucCF, \"`branchName`\", branchNumber, \"`defaultModel`\", modelList, algn_len, replace_tree, origRes, orig_bic)";
    //fprintf(stdout, "\n");
    //fprintf(stdout, mpiCommandString);
    //fprintf(stdout, "\n");
    //ExecuteCommands(mpiCommandString);
    //MPISend(mpiNodeI + 1, "optimizeBranchOmegas(`lfID`, `nucCF`, `branchName`, `branchNumber`, `defaultModel`, `modelList`, `algn_len`, `replace_tree`, `origRes`, `orig_bic`)");
    // XXX So the problem is that variables don't transfer to the MPI nodes.
    // Nor do objects. So the modificaiton needs to be done locally and the
    // optimization distally.
    omg_best_omega_ever = 3;
    mpiCommandString = "return "+ omg_best_omega_ever + ";";
    MPISend(mpiNodeI + 1, mpiCommandString);
    return 0;
}

//------------------------------------------------------------------------------------------------------------------------
function receiveAnMPIjob ()
{
    MPIReceive (-1, fromNode, result);
    fromNode += (-1);
    doneID = _MPI_NODE_STATUS[fromNode]-1;
    _MPI_NODE_STATUS[fromNode] = -1;
    //tbr = {};
    //tbr[0] = doneID;
    //tbr[1] = result;
    best_models[doneID] = result;
    return fromNode;
}

//------------------------------------------------------------------------------------------------------------------------
function calculateBranchLengthByName (modelList, bestModels, treeName, branchName, branchNumber)
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
        ExecuteCommands ("Model M" + mI + " = (MGMatrix" + 1 + ", codon3x4, 0)");
        models[mI] = "M" + mI;
    }

    mls = {};
    for (mI = 1; mI <= bestModels[branchNumber]; mI = mI + 1)
    {
        ExecuteCommands ("GetString (M" + mI + "L, M" + mI + ", -1);");
        //GetString(tempML, modelList[mI], -1);
        //mls[mI] = tempML;
        mls[mI] = "M" + mI + "L";
        //fprintf(stdout, tempML);
    }

    t = Eval("`treeName`.`branchName`.t");
    omegas = {};
    pauxs = {};
    evalstring = "(";
    for (oI = 1; oI <= bestModels[branchNumber]; oI = oI + 1)
    {
        // String
        if (oI > 1)
        {
            evalstring = evalstring + "+";
        }
        //evalstring = evalstring + "(" + mls[oI] + ")";
        //evalstring = evalstring + "(" + mls[oI] + ")";
        //evalstring = evalstring + "(" + mls[oI] + ")";
        //ExecuteCommands("evalstring = evalstring + \"(" + mls[oI] + ")\";");
        ExecuteCommands("tempML = " + mls[oI] + ";");
        ExecuteCommands("evalstring = evalstring + \"(" + tempML + ")\";");
        //evalstring = evalstring + "(" + mls[oI] + ")";
        if (oI != bestModels[branchNumber])
        {
            evalstring = evalstring + "*Paux" + oI;
        }
        for (prevOI = oI - 1; prevOI > 0; prevOI = prevOI - 1)
        {
            evalstring = evalstring + "*(1-Paux" + prevOI + ")";
        }

        // Values
        omegas[oI] = Eval("`treeName`.`branchName`.omega" + oI);
        if (oI < bestModels[branchNumber])
        {
            pauxs[oI] = Eval("`treeName`.`branchName`.Paux" + oI);
        }
    }
    //evalstring = evalstring + ")/" + bestModels[branchNumber];
    evalstring = evalstring + ")/3";


    fprintf(stdout, "\n");
    fprintf(stdout, branchName);
    fprintf(stdout, ":\n");

    fprintf(stdout, "\n");
    fprintf(stdout, evalstring);
    fprintf(stdout, "\n");
    branchLength = Eval(evalstring);
    fprintf(stdout, "\n");
    fprintf(stdout, branchLength);
    fprintf(stdout, "\n");

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
