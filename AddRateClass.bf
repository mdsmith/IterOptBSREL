LoadFunctionLibrary("DistributionFunctions.bf");

function addRate2Branch(lfID, nucCF, branchName, defaultModel, modelList, algn_len, replace_tree)
{
    model_assignments = {};
    branch_names = {};
    addRate2BranchAdvanced(lfID, nucCF, branchName, defaultModel, modelList, algn_len, replace_tree, branch_names, model_assignments);
    return 0;
}

// Testing with nucCF passed in
function addRate2BranchAdvanced(lfID, nucCF, branchName, defaultModel, modelList, algn_len, replace_tree, branch_names, model_assignments)
{
    ExecuteCommands ("GetString(lfInfoAA, " + lfID + ", -1)");
    lfTree = lfInfoAA["Trees"];
    lfTreeID = lfTree[0];
    ExecuteCommands ("orig_tree_string = Format(" + lfTreeID + ",1,1)");

    lfdsf = lfInfoAA["Datafilters"];
    lfdsfID = lfdsf[0];
    ExecuteCommands ("HarvestFrequencies (lfnuc3," + lfdsfID + ", 3, 1, 1)");
    lfnucCF = nucCF;

    currentParams = {};
    paramProportions = {};

    // Check and store previous Omegas
    numOmegas = 1;
    nextOmega = 1;
    openOmegaFound = 0;
    while (openOmegaFound != 1)
    {
        // Grab the previous omega. If it exists omegaInfo will be a matrix of length three.
        // Otherwise it will be a matrix of length one (I think).
        ExecuteCommands("GetInformation(omegaInfo, " + lfTreeID + "." + branchName + ".omega" + numOmegas + ")");
        // Check to see if this previous omega exists.
        if (Columns(omegaInfo) == 3)
        {
            // Alright, it does. Now lets get information regarding its proportion (remember,
            // the number of proportions is (the number of rate classes) - 1.
            prevOmega = numOmegas - 1;
            if (prevOmega >= 1)
            {
                ExecuteCommands("GetInformation(PauxInfo, " + lfTreeID + "." + branchName + ".Paux" + prevOmega + ")");
                // Save the proportion (part)
                paramProportions[prevOmega] = PauxInfo[0];
            }
            // Save the rate class rate
            currentParams[numOmegas] = omegaInfo[0];
            // Iterate the future rate class identity counter
            nextOmega = nextOmega + 1;
        }
        else
        {
            openOmegaFound = 1;
        }
        numOmegas = numOmegas + 1;
    }
    fprintf(stdout, "\n");
    fprintf(stdout, "Last set's values:\n");
    fprintf(stdout, currentParams);
    fprintf(stdout, "\n");
    fprintf(stdout, paramProportions);
    fprintf(stdout, "\n");

    // So I guess if nextOmega = 1 then the original n/ns rate needs to
    // become omega1 and there needs to be an additional omega2 rate class.
    // This will of course mean nextOmega = 2, and this need to make sense
    // going forward
    if (nextOmega == 1)
    {

        currentParams[1] = 10;
        currentParams[2] = currentParams[1]*2;
        paramProportions[1] = (algn_len-1)/algn_len;
        srate = Eval (lfTreeID + "." + branchName + ".syn");
        nsrate = Eval (lfTreeID + "." + branchName + ".nonsyn");
        if (srate > 0)
        {
            //currentParams[0] = Min (10, nsrate/srate);
            currentParams[1] = Min (10, nsrate/srate);
            currentParams[2] = currentParams[1]*2;
        }
        nextOmega = 2;
    }
    else
    {
        // in the event that this is the first Paux, make this at least
        // one site smaller than all of the sites. For later sites this
        // will diminish the proportion by less than a site, so it isn't
        // a big deal.
        paramProportions[nextOmega - 1] = ((algn_len -1)/algn_len);

        avail_site_found = 0;

        // I believe we can search back to nextOmega - 1 as well, instead of having to start with
        // nextOmega - 2. This is because the proportion of nextOmega - 2 was optimized to take
        // into account the number of sites that had to go to nextOmega - 1 when its Paux was
        // essentially 1
        for (site_source_I = nextOmega - 1; site_source_I > 0; site_source_I = site_source_I - 1)
        {
            if (avail_site_found == 0)
            {
                this_prop = paramProportions[site_source_I];
                for (remI = site_source_I - 1; remI > 0; remI = remI - 1)
                {
                    this_prop = this_prop * (1 - paramProportions[remI]);
                }
                // this_prop now contains the proportion of site_source_I omega.
                // If this proportion is greater than two sites, take one and
                // break the loop. Otherwise keep looking.
                if (this_prop > (2 * (1/algn_len)))
                {
                    avail_site_found = 1;
                    paramProportions[site_source_I] = (this_prop - (1/algn_len))/this_prop;
                }
            }
        }
        if (avail_site_found == 0)
        {
            fprintf(stdout, "Impossibly, no omega has a proportion greater than one site...");
        }

        currentParams[nextOmega] = currentParams[nextOmega - 1] * 2;  // XXX this should be selected
                                        // more intelligently
    }

    // Build a model matrix for each rate class to be used (in proportion)
    // in the LF.
    for (matrixI = 1; matrixI <= nextOmega; matrixI = matrixI + 1)
    {
        ExecuteCommands("PopulateModelMatrix(\"MGMatrix" + matrixI + "\", lfnucCF, \"t\", \"omega" + matrixI + "\", \"\");");
    }

    lfModel = lfInfoAA["Models"];
    lfModelID = lfModel[0];

    lfbaseFreqs = lfInfoAA["Base frequencies"];
    lfbaseFreqsID = lfbaseFreqs[0];

    matrixString = "";
    for (paramI = 1; paramI <= nextOmega; paramI = paramI + 1)
    {
        if (paramI > 1)
        {
            matrixString = matrixString + "+";
        }
        // We need the matrix from the rate class
        matrixString = matrixString + "Exp(MGMatrix" + paramI + ")";
        // However we don't store the proportion for the last rate class (DOF = # rate classes - 1)
        if (paramI != nextOmega)
        {
            matrixString = matrixString + "*Paux" + paramI;
        }
        // As we don't store the proportion for the last rate class, we need to be able to determine
        // that proportion from the other proportions. This relationship is used below:
        for (prevParamI = paramI - 1; prevParamI > 0; prevParamI = prevParamI - 1)
        {
            matrixString = matrixString + "*(1-Paux" + prevParamI + ")";
        }
    }
    fprintf(stdout, "matrixString:\n");
    fprintf(stdout, "\n");
    fprintf(stdout, matrixString);
    fprintf(stdout, "\n");
    //ExecuteCommands ("Model BSREL = (matrixString, " + lfbaseFreqsID + ", EXPLICIT_FORM_MATRIX_EXPONENTIAL);");


    // XXX There is a pretty good chance that much of the code above is
    // similarly conditioned, but sorting it from what must be computed to
    // set branch parameters is for another time.
    if (modelList[nextOmega] == 0)
    {
        ExecuteCommands ("Model BSREL" + nextOmega + " = (matrixString, " + lfbaseFreqsID + ", EXPLICIT_FORM_MATRIX_EXPONENTIAL);");
        modelList[nextOmega] = "BSREL" + nextOmega;
    }

    // But it would be nice to have a BSREL1.
    // Testing
    /*
    */
    if (modelList[1] == 0)
    {
        // All we need is a matrix string
        matrixString = "Exp(MGMatrix1)";
        ExecuteCommands ("Model BSREL1 = (matrixString, " + lfbaseFreqsID + ", EXPLICIT_FORM_MATRIX_EXPONENTIAL);");
        modelList[1] = "BSREL1";
    }

    new_tree_string = orig_tree_string;
    ExecuteCommands ("UseModel(" + defaultModel + ")");

    //new_tree_string = new_tree_string^{{branchName, branchName + "{BSREL}"}};

    // Yeah, model_assignments should be of length zero or of length equal to
    // the number of branches. I may implement a check for that later.
    for (mod_assgn_I = 0; mod_assgn_I < Abs(model_assignments); mod_assgn_I = mod_assgn_I + 1)
    {
        if (model_assignments[mod_assgn_I] > 1)
        {
            new_tree_string = new_tree_string^{{branch_names[mod_assgn_I] + ":", branch_names[mod_assgn_I] + "{BSREL" + model_assignments[mod_assgn_I] + "}:"}};
        }
    }

    if (Abs(model_assignments) == 0)
    {
        new_tree_string = new_tree_string^{{branchName + ":", branchName + "{BSREL" + nextOmega + "}:"}};
    }
    fprintf(stdout, new_tree_string);
    fprintf(stdout, "\n");

    if (replace_tree == 1)
    {
        REPLACE_TREE_STRUCTURE = 1;
        fprintf(stdout, "\nThe tree should be replaced\n");
    }
    else
    {
        REPLACE_TREE_STRUCTURE = 0;
    }

    ExecuteCommands ("Tree `lfTreeID` = `new_tree_string`");
    ExecuteCommands ("LikelihoodFunction `lfID` = (`lfdsfID`, `lfTreeID`);");


    fprintf(stdout, "\n");
    fprintf(stdout, "Initial values:\n");
    fprintf(stdout, currentParams);
    fprintf(stdout, "\n");
    fprintf(stdout, paramProportions);
    fprintf(stdout, "\n");

    for (omegaI = 1; omegaI <= nextOmega; omegaI = omegaI + 1)
    {
        ExecuteCommands(lfTreeID + "." + branchName + ".omega" + omegaI + " = " + currentParams[omegaI] + ";");
        if (omegaI != nextOmega)
        {
            ExecuteCommands(lfTreeID + "." + branchName + ".Paux" + omegaI + " = " + paramProportions[omegaI] + ";");
            ExecuteCommands(lfTreeID + "." + branchName + ".Paux" + omegaI + " :< 1;");
        }
    }

    return 0;
}

// XXX this initialization regimen is producing subpar results
function assignModels2Branches(lfID, nucCF, defaultModel, branch_names, model_assignments, algn_len, model_list)
{

    ExecuteCommands ("UseModel(" + defaultModel + ")");

    ExecuteCommands ("GetString(lfInfoAA, " + lfID + ", -1)");
    lfTree = lfInfoAA["Trees"];
    lfTreeID = lfTree[0];
    lfdsf = lfInfoAA["Datafilters"];
    lfdsfID = lfdsf[0];
    ExecuteCommands ("HarvestFrequencies (lfnuc3," + lfdsfID + ", 3, 1, 1)");
    lfnucCF = nucCF;
    lfbaseFreqs = lfInfoAA["Base frequencies"];
    lfbaseFreqsID = lfbaseFreqs[0];
    //lfModel = lfInfoAA["Models"];
    //lfModelID = lfModel[0];

    // Check to make sure all of the models exist. If not, create them and
    // add them to model_list. This should eventually be abstracted into a
    // separate function
    for (mI = 0; mI < Abs(model_assignments); mI = mI + 1)
    {
        //fprintf(stdout, model_list);
        //fprintf(stdout, "\n");
        //fprintf(stdout, "Model list model: " + model_assignments[mI] + "\n");
        if (model_list[model_assignments[mI]] == 0)
        {
            for (matrixI = 1; matrixI <= model_assignments[mI]; matrixI = matrixI + 1)
            {
                ExecuteCommands("PopulateModelMatrix(\"MGMatrix" + matrixI + "\", lfnucCF, \"t\", \"omega" + matrixI + "\", \"\");");
            }
            matrixString = "";
            for (paramI = 1; paramI <= model_assignments[mI]; paramI = paramI + 1)
            {
                if (paramI > 1)
                {
                    matrixString = matrixString + "+";
                }
                // We need the matrix from the rate class
                matrixString = matrixString + "Exp(MGMatrix" + paramI + ")";
                // However we don't store the proportion for the last rate class (DOF = # rate classes - 1)
                if (paramI != model_assignments[mI])
                {
                    matrixString = matrixString + "*Paux" + paramI;
                }
                // As we don't store the proportion for the last rate class, we need to be able to determine
                // that proportion from the other proportions. This relationship is used below:
                for (prevParamI = paramI - 1; prevParamI > 0; prevParamI = prevParamI - 1)
                {
                    matrixString = matrixString + "*(1-Paux" + prevParamI + ")";
                }
            }
            ExecuteCommands ("Model BSREL" + model_assignments[mI] + " = (matrixString, " + lfbaseFreqsID + ", EXPLICIT_FORM_MATRIX_EXPONENTIAL);");
            modelList[nextOmega] = "BSREL" + model_assignments[mI];
        }
    }


    if (Abs(temp_omegas) == 0)
    {
        temp_omegas = {};
        for (mI = 0; mI < Abs(model_assignments); mI = mI + 1)
        {
            temp_omegas[mI] = 10;
            srate = Eval (lfTreeID + "." + branch_names[mI] + ".syn");
            nsrate = Eval (lfTreeID + "." + branch_names[mI] + ".nonsyn");
            fprintf(stdout, "" + nsrate + ", " + srate + "\n");
            if (srate > 0)
            {
                temp_omegas[mI] = Min (10, nsrate/srate);
            }
        }
    }

    fprintf(stdout, "\nomega results: \n");
    fprintf(stdout, temp_omegas);
    fprintf(stdout, "\n");

    ExecuteCommands ("orig_tree_string = Format(" + lfTreeID + ",1,1)");
    final_tree_string = orig_tree_string;
    for (mod_assgn_I = 0; mod_assgn_I < Abs(model_assignments); mod_assgn_I = mod_assgn_I + 1)
    {
        if (model_assignments[mod_assgn_I] > 1)
        {
            final_tree_string = final_tree_string^{{branch_names[mod_assgn_I] + ":", branch_names[mod_assgn_I] + "{BSREL" + model_assignments[mod_assgn_I] + "}:"}};
        }
    }
    fprintf(stdout, final_tree_string);
    fprintf(stdout, "\n");

    REPLACE_TREE_STRUCTURE = 1;

    ExecuteCommands ("Tree `lfTreeID` = `final_tree_string`");
    ExecuteCommands ("LikelihoodFunction `lfID` = (`lfdsfID`, `lfTreeID`);");

    for (mod_assgn_I = 0; mod_assgn_I < Abs(model_assignments); mod_assgn_I = mod_assgn_I + 1)
    {
        if (model_assignments[mod_assgn_I] > 0)
        {
            for (omegaI = 1; omegaI <= model_assignments[omegaI]; omegaI = omegaI + 1)
            {
                ExecuteCommands(lfTreeID + "." + branch_names[mod_assgn_I] + ".omega" + omegaI + " = " + (temp_omegas[mod_assgn_I] * (2^omegaI)) + ";");
                if (omegaI != model_assignments[mod_assgn_I])
                {
                    //fprintf(stdout, "\nomega: " + omegaI + " opt: " + model_assignments[mod_assgn_I] + " branch: " + branch_names[mod_assgn_I] + "\n");
                    ExecuteCommands(lfTreeID + "." + branch_names[mod_assgn_I] + ".Paux" + omegaI + " = " + ((algn_len - 1)/algn_len) + ";");
                    ExecuteCommands(lfTreeID + "." + branch_names[mod_assgn_I] + ".Paux" + omegaI + " :< 1;");
                }
            }
            if (model_assignments[mod_assgn_I] > 1)
            {
                ExecuteCommands(lfTreeID + "." + branch_names[mod_assgn_I] + ".Paux1 = " + (1 - (model_assignments[mod_assgn_I] * (1/algn_len))) + ";");
            }
        }
    }

    return 0;
}

return 0;
