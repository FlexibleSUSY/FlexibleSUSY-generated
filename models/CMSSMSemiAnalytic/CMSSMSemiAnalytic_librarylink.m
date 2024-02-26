Print["================================"];
Print["FlexibleSUSY 2.8.0"];
Print["CMSSMSemiAnalytic"];
Print["http://flexiblesusy.hepforge.org"];
Print["================================"];

libCMSSMSemiAnalytic = FileNameJoin[{Directory[], "models", "CMSSMSemiAnalytic", "CMSSMSemiAnalytic_librarylink.so"}];

FSCMSSMSemiAnalyticGetSettings = LibraryFunctionLoad[libCMSSMSemiAnalytic, "FSCMSSMSemiAnalyticGetSettings", LinkObject, LinkObject];
FSCMSSMSemiAnalyticGetSMInputParameters = LibraryFunctionLoad[libCMSSMSemiAnalytic, "FSCMSSMSemiAnalyticGetSMInputParameters", LinkObject, LinkObject];
FSCMSSMSemiAnalyticGetInputParameters = LibraryFunctionLoad[libCMSSMSemiAnalytic, "FSCMSSMSemiAnalyticGetInputParameters", LinkObject, LinkObject];
FSCMSSMSemiAnalyticGetProblems = LibraryFunctionLoad[libCMSSMSemiAnalytic, "FSCMSSMSemiAnalyticGetProblems", LinkObject, LinkObject];
FSCMSSMSemiAnalyticGetWarnings = LibraryFunctionLoad[libCMSSMSemiAnalytic, "FSCMSSMSemiAnalyticGetWarnings", LinkObject, LinkObject];
FSCMSSMSemiAnalyticToSLHA = LibraryFunctionLoad[libCMSSMSemiAnalytic, "FSCMSSMSemiAnalyticToSLHA", LinkObject, LinkObject];

FSCMSSMSemiAnalyticOpenHandleLib = LibraryFunctionLoad[libCMSSMSemiAnalytic, "FSCMSSMSemiAnalyticOpenHandle", {{Real,1}}, Integer];
FSCMSSMSemiAnalyticCloseHandle = LibraryFunctionLoad[libCMSSMSemiAnalytic, "FSCMSSMSemiAnalyticCloseHandle", {Integer}, Void];

FSCMSSMSemiAnalyticSetLib = LibraryFunctionLoad[libCMSSMSemiAnalytic, "FSCMSSMSemiAnalyticSet", {Integer, {Real,1}}, Void];

FSCMSSMSemiAnalyticCalculateSpectrum = LibraryFunctionLoad[libCMSSMSemiAnalytic, "FSCMSSMSemiAnalyticCalculateSpectrum", LinkObject, LinkObject];
FSCMSSMSemiAnalyticCalculateObservables = LibraryFunctionLoad[libCMSSMSemiAnalytic, "FSCMSSMSemiAnalyticCalculateObservables", LinkObject, LinkObject];


FSCMSSMSemiAnalyticCalculateSpectrum::error = "`1`";
FSCMSSMSemiAnalyticCalculateSpectrum::warning = "`1`";

FSCMSSMSemiAnalyticCalculateObservables::error = "`1`";
FSCMSSMSemiAnalyticCalculateObservables::warning = "`1`";


FSCMSSMSemiAnalytic::info = "`1`";
FSCMSSMSemiAnalytic::nonum = "Error: `1` is not a numeric input value!";
FSCMSSMSemiAnalyticMessage[s_] := Message[FSCMSSMSemiAnalytic::info, s];

FSCMSSMSemiAnalyticCheckIsNumeric[a_?NumericQ] := a;
FSCMSSMSemiAnalyticCheckIsNumeric[a_] := (Message[FSCMSSMSemiAnalytic::nonum, a]; Abort[]);

fsDefaultSettings = {
      precisionGoal -> 1.*^-4,           (* FlexibleSUSY[0] *)
      maxIterations -> 0,                (* FlexibleSUSY[1] *)
      solver -> 2,     (* FlexibleSUSY[2] *)
      calculateStandardModelMasses -> 0, (* FlexibleSUSY[3] *)
      poleMassLoopOrder -> 4,            (* FlexibleSUSY[4] *)
      ewsbLoopOrder -> 4,                (* FlexibleSUSY[5] *)
      betaFunctionLoopOrder -> 4,        (* FlexibleSUSY[6] *)
      thresholdCorrectionsLoopOrder -> 4,(* FlexibleSUSY[7] *)
      higgs2loopCorrectionAtAs -> 1,     (* FlexibleSUSY[8] *)
      higgs2loopCorrectionAbAs -> 1,     (* FlexibleSUSY[9] *)
      higgs2loopCorrectionAtAt -> 1,     (* FlexibleSUSY[10] *)
      higgs2loopCorrectionAtauAtau -> 1, (* FlexibleSUSY[11] *)
      forceOutput -> 0,                  (* FlexibleSUSY[12] *)
      topPoleQCDCorrections -> 3,        (* FlexibleSUSY[13] *)
      betaZeroThreshold -> 1.*^-11,      (* FlexibleSUSY[14] *)
      forcePositiveMasses -> 0,          (* FlexibleSUSY[16] *)
      poleMassScale -> 0,                (* FlexibleSUSY[17] *)
      eftPoleMassScale -> 0,             (* FlexibleSUSY[18] *)
      eftMatchingScale -> 0,             (* FlexibleSUSY[19] *)
      eftMatchingLoopOrderUp -> 2,       (* FlexibleSUSY[20] *)
      eftMatchingLoopOrderDown -> 1,     (* FlexibleSUSY[21] *)
      eftHiggsIndex -> 0,                (* FlexibleSUSY[22] *)
      calculateBSMMasses -> 1,           (* FlexibleSUSY[23] *)
      thresholdCorrections -> 124111421, (* FlexibleSUSY[24] *)
      higgs3loopCorrectionRenScheme -> 0,(* FlexibleSUSY[25] *)
      higgs3loopCorrectionAtAsAs -> 1,   (* FlexibleSUSY[26] *)
      higgs3loopCorrectionAbAsAs -> 1,   (* FlexibleSUSY[27] *)
      higgs3loopCorrectionAtAtAs -> 1,   (* FlexibleSUSY[28] *)
      higgs3loopCorrectionAtAtAt -> 1,   (* FlexibleSUSY[29] *)
      higgs4loopCorrectionAtAsAsAs -> 1, (* FlexibleSUSY[30] *)
      loopLibrary -> 0,                  (* FlexibleSUSY[31] *)
      calculateAMM -> 2.0,               (* FlexibleSUSY[32] *)
      parameterOutputScale -> 0          (* MODSEL[12] *)
};

fsDefaultSMParameters = {
    alphaEmMZ -> 1/127.916, (* SMINPUTS[1] *)
    GF -> 1.1663787*^-5,    (* SMINPUTS[2] *)
    alphaSMZ -> 0.1184,     (* SMINPUTS[3] *)
    MZ -> 91.1876,          (* SMINPUTS[4] *)
    mbmb -> 4.18,           (* SMINPUTS[5] *)
    Mt -> 173.34,           (* SMINPUTS[6] *)
    Mtau -> 1.777,          (* SMINPUTS[7] *)
    Mv3 -> 0,               (* SMINPUTS[8] *)
    MW -> 80.385,           (* SMINPUTS[9] *)
    Me -> 0.000510998902,   (* SMINPUTS[11] *)
    Mv1 -> 0,               (* SMINPUTS[12] *)
    Mm -> 0.1056583715,     (* SMINPUTS[13] *)
    Mv2 -> 0,               (* SMINPUTS[14] *)
    md2GeV -> 0.00475,      (* SMINPUTS[21] *)
    mu2GeV -> 0.0024,       (* SMINPUTS[22] *)
    ms2GeV -> 0.104,        (* SMINPUTS[23] *)
    mcmc -> 1.27,           (* SMINPUTS[24] *)
    CKMTheta12 -> 0,
    CKMTheta13 -> 0,
    CKMTheta23 -> 0,
    CKMDelta -> 0,
    PMNSTheta12 -> 0,
    PMNSTheta13 -> 0,
    PMNSTheta23 -> 0,
    PMNSDelta -> 0,
    PMNSAlpha1 -> 0,
    PMNSAlpha2 -> 0,
    alphaEm0 -> 1/137.035999074,
    Mh -> 125.09
};

fsCMSSMSemiAnalyticDefaultInputParameters = {
   m12 -> 0,
   TanBeta -> 0,
   Azero -> 0,
   MuInput -> 0
};

Options[FSCMSSMSemiAnalyticOpenHandle] = {
    Sequence @@ fsDefaultSettings,
    Sequence @@ fsDefaultSMParameters,
    Sequence @@ fsCMSSMSemiAnalyticDefaultInputParameters

};

FSCMSSMSemiAnalyticOpenHandle[a___, (fsSettings | fsSMParameters | fsModelParameters ) -> s_List, r___] :=
    FSCMSSMSemiAnalyticOpenHandle[a, Sequence @@ s, r];

FSCMSSMSemiAnalyticOpenHandle[OptionsPattern[]] :=
    FSCMSSMSemiAnalyticOpenHandleLib[
        FSCMSSMSemiAnalyticCheckIsNumeric /@ {
            (* spectrum generator settings *)
            OptionValue[precisionGoal],
            OptionValue[maxIterations],
            OptionValue[solver],
            OptionValue[calculateStandardModelMasses],
            OptionValue[poleMassLoopOrder],
            OptionValue[ewsbLoopOrder],
            OptionValue[betaFunctionLoopOrder],
            OptionValue[thresholdCorrectionsLoopOrder],
            OptionValue[higgs2loopCorrectionAtAs],
            OptionValue[higgs2loopCorrectionAbAs],
            OptionValue[higgs2loopCorrectionAtAt],
            OptionValue[higgs2loopCorrectionAtauAtau],
            OptionValue[forceOutput],
            OptionValue[topPoleQCDCorrections],
            OptionValue[betaZeroThreshold],
            OptionValue[forcePositiveMasses],
            OptionValue[poleMassScale],
            OptionValue[eftPoleMassScale],
            OptionValue[eftMatchingScale],
            OptionValue[eftMatchingLoopOrderUp],
            OptionValue[eftMatchingLoopOrderDown],
            OptionValue[eftHiggsIndex],
            OptionValue[calculateBSMMasses],
            OptionValue[thresholdCorrections],
            OptionValue[higgs3loopCorrectionRenScheme],
            OptionValue[higgs3loopCorrectionAtAsAs],
            OptionValue[higgs3loopCorrectionAbAsAs],
            OptionValue[higgs3loopCorrectionAtAtAs],
            OptionValue[higgs3loopCorrectionAtAtAt],
            OptionValue[higgs4loopCorrectionAtAsAsAs],
            OptionValue[loopLibrary],
            OptionValue[calculateAMM],
            OptionValue[parameterOutputScale],

            (* Standard Model input parameters *)
            OptionValue[alphaEmMZ],
            OptionValue[GF],
            OptionValue[alphaSMZ],
            OptionValue[MZ],
            OptionValue[mbmb],
            OptionValue[Mt],
            OptionValue[Mtau],
            OptionValue[Mv3],
            OptionValue[MW],
            OptionValue[Me],
            OptionValue[Mv1],
            OptionValue[Mm],
            OptionValue[Mv2],
            OptionValue[md2GeV],
            OptionValue[mu2GeV],
            OptionValue[ms2GeV],
            OptionValue[mcmc],
            OptionValue[CKMTheta12],
            OptionValue[CKMTheta13],
            OptionValue[CKMTheta23],
            OptionValue[CKMDelta],
            OptionValue[PMNSTheta12],
            OptionValue[PMNSTheta13],
            OptionValue[PMNSTheta23],
            OptionValue[PMNSDelta],
            OptionValue[PMNSAlpha1],
            OptionValue[PMNSAlpha2],
            OptionValue[alphaEm0],
            OptionValue[Mh]

            (* CMSSMSemiAnalytic input parameters *)
            ,
            OptionValue[m12],
            OptionValue[TanBeta],
            OptionValue[Azero],
            OptionValue[MuInput]

        }
];

Options[FSCMSSMSemiAnalyticSet] = Options[FSCMSSMSemiAnalyticOpenHandle];

FSCMSSMSemiAnalyticSet[handle_Integer, a___, (fsSettings | fsSMParameters | fsModelParameters) -> s_List, r___] :=
    FSCMSSMSemiAnalyticSet[handle, a, Sequence @@ s, r];

FSCMSSMSemiAnalyticSet[handle_Integer, p:OptionsPattern[]] :=
    FSCMSSMSemiAnalyticSetLib[
        handle,
        ReleaseHold[Hold[FSCMSSMSemiAnalyticCheckIsNumeric /@ {
            (* spectrum generator settings *)
            OptionValue[precisionGoal],
            OptionValue[maxIterations],
            OptionValue[solver],
            OptionValue[calculateStandardModelMasses],
            OptionValue[poleMassLoopOrder],
            OptionValue[ewsbLoopOrder],
            OptionValue[betaFunctionLoopOrder],
            OptionValue[thresholdCorrectionsLoopOrder],
            OptionValue[higgs2loopCorrectionAtAs],
            OptionValue[higgs2loopCorrectionAbAs],
            OptionValue[higgs2loopCorrectionAtAt],
            OptionValue[higgs2loopCorrectionAtauAtau],
            OptionValue[forceOutput],
            OptionValue[topPoleQCDCorrections],
            OptionValue[betaZeroThreshold],
            OptionValue[forcePositiveMasses],
            OptionValue[poleMassScale],
            OptionValue[eftPoleMassScale],
            OptionValue[eftMatchingScale],
            OptionValue[eftMatchingLoopOrderUp],
            OptionValue[eftMatchingLoopOrderDown],
            OptionValue[eftHiggsIndex],
            OptionValue[calculateBSMMasses],
            OptionValue[thresholdCorrections],
            OptionValue[higgs3loopCorrectionRenScheme],
            OptionValue[higgs3loopCorrectionAtAsAs],
            OptionValue[higgs3loopCorrectionAbAsAs],
            OptionValue[higgs3loopCorrectionAtAtAs],
            OptionValue[higgs3loopCorrectionAtAtAt],
            OptionValue[higgs4loopCorrectionAtAsAsAs],
            OptionValue[loopLibrary],
            OptionValue[calculateAMM],
            OptionValue[parameterOutputScale],

            (* Standard Model input parameters *)
            OptionValue[alphaEmMZ],
            OptionValue[GF],
            OptionValue[alphaSMZ],
            OptionValue[MZ],
            OptionValue[mbmb],
            OptionValue[Mt],
            OptionValue[Mtau],
            OptionValue[Mv3],
            OptionValue[MW],
            OptionValue[Me],
            OptionValue[Mv1],
            OptionValue[Mm],
            OptionValue[Mv2],
            OptionValue[md2GeV],
            OptionValue[mu2GeV],
            OptionValue[ms2GeV],
            OptionValue[mcmc],
            OptionValue[CKMTheta12],
            OptionValue[CKMTheta13],
            OptionValue[CKMTheta23],
            OptionValue[CKMDelta],
            OptionValue[PMNSTheta12],
            OptionValue[PMNSTheta13],
            OptionValue[PMNSTheta23],
            OptionValue[PMNSDelta],
            OptionValue[PMNSAlpha1],
            OptionValue[PMNSAlpha2],
            OptionValue[alphaEm0],
            OptionValue[Mh]

            (* CMSSMSemiAnalytic input parameters *)
            ,
            OptionValue[m12],
            OptionValue[TanBeta],
            OptionValue[Azero],
            OptionValue[MuInput]

        }] /. HoldPattern[OptionValue[param_]] :> param /.
        { p } /.
        FSCMSSMSemiAnalyticGetSettings[handle] /.
        FSCMSSMSemiAnalyticGetSMInputParameters[handle] /.
        FSCMSSMSemiAnalyticGetInputParameters[handle]]];
