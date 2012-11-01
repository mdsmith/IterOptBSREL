//VERBOSITY_LEVEL = 0;

//chickenWings = 0.3;
//GetInformation(numberInfo, chickenWings);
//fprintf(stdout, numberInfo[0] + "\n");

function addRate2BranchNumber(lfID, branchName)
{
    ExecuteCommands ("GetString(lfInfoAA, " + lfID + ", -1)");
    //fprintf(stdout, lfInfoAA);
    lfTree = lfInfoAA["Trees"];
    lfTreeID = lfTree[0];
    ExecuteCommands ("orig_tree_string = Format(" + lfTreeID + ",1,1)");

    lfdsf = lfInfoAA["Datafilters"];
    lfdsfID = lfdsf[0];
    ExecuteCommands ("HarvestFrequencies (lfnuc3," + dsfID + ", 3, 1, 1)");
    lfnucCF = CF3x4 (lfnuc3, GeneticCodeExclusions);


    //fprintf(stdout, "Branch: \n");
    //ExecuteCommands ( "fprintf(stdout, " + lfTreeID + "." + branchName + ".nonsyn)");
    //ExecuteCommands ( "branch_nsrate = Eval(\"" + lfTreeID + "." + branchName + ".nonsyn)");

    currentParams = {};
    paramProportions = {};


    presOmegas = 1;
    nextOmega = 1;
    // Check and store previous Omegas
    for (presOmegas = 1; presOmegas < 6; presOmegas = presOmegas + 1)
    {
        // Grab the previous omega. If it exists omegaInfo will be a matrix of length three.
        // Otherwise it will be a matrix of length one (I think).
        ExecuteCommands("GetInformation(omegaInfo, " + lfTreeID + "." + branchName + ".omega" + presOmegas + ")");
        // Check to see if this previous omega exists.
        if (Abs(omegaInfo) == 3)
        {
            // Alright, it does. Now lets get information regarding its proportion (remember,
            // the number of proportions is (the number of rate classes) - 1.
            prevOmega = presOmegas - 1;
            if (prevOmega >= 1)
            {
                ExecuteCommands("GetInformation(PauxInfo, " + lfTreeID + "." + branchName + ".Paux" + prevOmega + ")");
                // Save the proportion (part)
                ParamProportions[prevOmega] = PauxInfo[0];
            }
            // Save the rate class rate
            currentParams[presOmegas] = omegaInfo[0];
            // Iterate the future rate class identity counter
            nextOmega = nextOmega + 1;
        }
    }
    newProportion = paramProportions[Abs(paramProportions)-1] * .9;
    ExecuteCommands(lfTreeID + "." + branchName + ".omega" + nextOmega + " = .9"); // XXX fix value
    ExecuteCommands(lfTreeID + "." + branchName + ".Paux" + nextOmega + " = " + newProportion); // XXX fix propotion

    for (matrixI = 1; matrixI < nextOmega; matrixI = matrixI + 1)
    {
        ExecuteCommands("PopulateModelMatrix(\"MGMatrix" + matrixI + "\", lfnucCF, \"t\", \"omega" + matrixI + "\", \"\")");
    }

    //branch_nsrate = Eval(lfTreeID + "." + branchName + ".omega1");
    //fprintf(stdout, branch_nsrate);
    //fprintf(stdout, "\n");

    lfModel = lfInfoAA["Models"];
    lfModelID = lfModel[0];
    //Parameter(s) missing in Model definition. Must have a matrix and a
    //compatible eqiulibrium frequencies vector.

    lfbaseFreqs = lfInfoAA["Base frequencies"];
    lfbaseFreqsID = lfdsf[0];
    //ExecuteCommands("lfdsf = CreateFilter(" + lfdfID + ",3,\"\",\"\",GeneticCodeExclusions)");

    // Define the new model
    //ExecuteCommands ("Model BSREL = " + lfModelID);
    matrixString = "";
    for (paramI = 1; paramI < Abs(currentParams); paramI = paramI + 1)
    {
        if (paramI > 1)
        {
            matrixString = matrixString + "+";
        }
        // We need the matrix from the rate class
        matrixString = matrixString + "Exp(MGMatrix" + paramI + ")"; // XXX source MGMatrix's
        // However we don't store the proportion for the last rate class (DOF = # rate classes - 1)
        if (paramI != Abs(currentParams) - 1)
        {
            matrixString = matrixString + "*Paux" + paramI;
        }
        // As we don't store the proportion for the last rate class, we need to be able to determine
        // that proportion from the other proportions. This relationship is used below:
        for (prevParamI = paramI - 1; paramI > 0; paramI = paramI - 1)
        {
            matrixString = matrixString + "*(1-pPaux" + prevParamI + ")";
        }
    }
    ExecuteCommands ("Model BSREL = (matrixString, " + lfbaseFreqsID + ", EXPLICIT_FORM_MATRIX_EXPONENTIAL)");

    new_tree_string = orig_tree_string;
    //list_of_models {};
    // correct import?:
    //totalBranchCount = BranchCount(orig_tree) + TipCount(orig_tree);
    //bNames = BranchName (orig_tree, -1);
    fprintf(stdout, new_tree_string);
    fprintf(stdout, "\n");
    //UseModel(MG94); // test to make sure this works...
    ExecuteCommands ("UseModel(" + lfModelID + ")");
    fprintf(stdout, new_tree_string);
    fprintf(stdout, "\n");

    new_tree_string = new_tree_string^{{branchName, branchName + "{BSREL}"}};
    fprintf(stdout, new_tree_string);
    fprintf(stdout, "\n");

    //for (k = 0; k < Abs(bNames); k = k+1)
    //{
        //if (bNames[k] == branchName)
        //{
            //// insert BSREL
        //}
    //}
    Tree new_branch_tree = new_tree_string;
    // should this be a modelID string or the model?
    //list_of_models[branchNumber] = BSREL
    // implement:
    //Tree new_branch_tree = setModels(orig_treeString, list_of_models);

    ExecuteCommands("LikelihoodFunction newLF = (" + lfdsfID + ", new_branch_tree)");

    return "newLF";
}

//function setModels(orig_treeString, list_of_models)
//{
    //Tree tbrT = orig_treeString;
    //for (k = 0; k < Abs(list_of_models);
    //{
        //
    //}
//
    //return tbrT
//}



return 0;
