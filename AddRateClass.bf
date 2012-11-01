//VERBOSITY_LEVEL = 0;

//chickenWings = 0.3;
//GetInformation(numberInfo, chickenWings);
//fprintf(stdout, numberInfo[0] + "\n");

function addRate2BranchNumber(lfID, branchName)
{
    ExecuteCommands ("GetString(lfInfoAA, " + lfID + ", -1)");
    lfTree = lfInfoAA["Trees"];
    lfTreeID = lfTree[0];
    ExecuteCommands ("orig_tree_string = Format(" + lfTreeID + ",1,1)");

    fprintf(stdout, "Branch: \n");
    //ExecuteCommands ( "fprintf(stdout, " + lfTreeID + "." + branchName + ".nonsyn)");
    //ExecuteCommands ( "branch_nsrate = Eval(\"" + lfTreeID + "." + branchName + ".nonsyn)");
    presOmegas = 1;
    nextOmega = 1;
    for (presOmegas = 1; presOmegas < 6; presOmegas = presOmegas + 1)
    {
        ExecuteCommands("GetInformation(omegaInfo, " + lfTreeID + "." + branchName + ".omega" + presOmegas + ")");
        if (Abs(omegaInfo) == 3)
        {
            nextOmega = nextOmega + 1;
        }
    }
    ExecuteCommands(lfTreeID + "." + branchName + ".omega" + nextOmega + " = .9");
    ExecuteCommands(lfTreeID + "." + branchName + ".Paux" + nextOmega + " = 1");
    // XXX OK, excellent, I can find out which omega to add. Now I
    // need to keep track of and harvest the omegas that are already
    // there so that I may build a new model out of them.
    //branch_nsrate = Eval(lfTreeID + "." + branchName + ".omega1");
    //fprintf(stdout, branch_nsrate);
    //fprintf(stdout, "\n");

    lfModel = lfInfoAA["Models"];
    lfModelID = lfModel[0];
    //Parameter(s) missing in Model definition. Must have a matrix and a
    //compatible eqiulibrium frequencies vector.

    // XXX Define the new model
    //ExecuteCommands ("Model BSREL = " + lfModelID);

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
    // XXX Tree new_branch_tree = new_tree_string;
    // should this be a modelID string or the model?
    //list_of_models[branchNumber] = BSREL
    // implement:
    //Tree new_branch_tree = setModels(orig_treeString, list_of_models);

    // XXX LikelihoodFunction newLF = (dsf, new_branch_tree);

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
