Print["================================"];
Print["FlexibleSUSY 2.6.1"];
Print["CE6SSM"];
Print["http://flexiblesusy.hepforge.org"];
Print["================================"];

libCE6SSM = FileNameJoin[{Directory[], "models", "CE6SSM", "CE6SSM_librarylink.so"}];

FSCE6SSMGetSettings = LibraryFunctionLoad[libCE6SSM, "FSCE6SSMGetSettings", LinkObject, LinkObject];
FSCE6SSMGetSMInputParameters = LibraryFunctionLoad[libCE6SSM, "FSCE6SSMGetSMInputParameters", LinkObject, LinkObject];
FSCE6SSMGetInputParameters = LibraryFunctionLoad[libCE6SSM, "FSCE6SSMGetInputParameters", LinkObject, LinkObject];
FSCE6SSMGetProblems = LibraryFunctionLoad[libCE6SSM, "FSCE6SSMGetProblems", LinkObject, LinkObject];
FSCE6SSMGetWarnings = LibraryFunctionLoad[libCE6SSM, "FSCE6SSMGetWarnings", LinkObject, LinkObject];
FSCE6SSMToSLHA = LibraryFunctionLoad[libCE6SSM, "FSCE6SSMToSLHA", LinkObject, LinkObject];

FSCE6SSMOpenHandleLib = LibraryFunctionLoad[libCE6SSM, "FSCE6SSMOpenHandle", {{Real,1}}, Integer];
FSCE6SSMCloseHandle = LibraryFunctionLoad[libCE6SSM, "FSCE6SSMCloseHandle", {Integer}, Void];

FSCE6SSMSetLib = LibraryFunctionLoad[libCE6SSM, "FSCE6SSMSet", {Integer, {Real,1}}, Void];

FSCE6SSMCalculateSpectrum = LibraryFunctionLoad[libCE6SSM, "FSCE6SSMCalculateSpectrum", LinkObject, LinkObject];
FSCE6SSMCalculateObservables = LibraryFunctionLoad[libCE6SSM, "FSCE6SSMCalculateObservables", LinkObject, LinkObject];
FSCE6SSMCalculateDecays = LibraryFunctionLoad[libCE6SSM, "FSCE6SSMCalculateDecays", LinkObject, LinkObject];

FSCE6SSMCalculateSpectrum::error = "`1`";
FSCE6SSMCalculateSpectrum::warning = "`1`";

FSCE6SSMCalculateObservables::error = "`1`";
FSCE6SSMCalculateObservables::warning = "`1`";

FSCE6SSMCalculateDecays::error = "`1`";
FSCE6SSMCalculateDecays::warning = "`1`";

FSCE6SSM::info = "`1`";
FSCE6SSM::nonum = "Error: `1` is not a numeric input value!";
FSCE6SSMMessage[s_] := Message[FSCE6SSM::info, s];

FSCE6SSMCheckIsNumeric[a_?NumericQ] := a;
FSCE6SSMCheckIsNumeric[a_] := (Message[FSCE6SSM::nonum, a]; Abort[]);

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

fdDefaultSettings = {
   minBRtoPrint -> 1*^-5,
   maxHigherOrderCorrections -> 4,
   alphaThomson -> 1,
   offShellVV -> 2
};

fsCE6SSMDefaultInputParameters = {
   TanBeta -> 0,
   m0SqGuess -> 0,
   m12Guess -> 0,
   AzeroGuess -> 0,
   LambdaInput -> 0,
   KappaInput -> 0,
   MuPrimeInput -> 0,
   BMuPrimeInput -> 0,
   vsInput -> 0,
   Lambda12Input -> 0
};

Options[FSCE6SSMOpenHandle] = {
    Sequence @@ fsDefaultSettings,
    Sequence @@ fsDefaultSMParameters,
    Sequence @@ fsCE6SSMDefaultInputParameters
   , Sequence @@ fdDefaultSettings
};

FSCE6SSMOpenHandle[a___, (fsSettings | fsSMParameters | fsModelParameters | fdSettings) -> s_List, r___] :=
    FSCE6SSMOpenHandle[a, Sequence @@ s, r];

FSCE6SSMOpenHandle[OptionsPattern[]] :=
    FSCE6SSMOpenHandleLib[
        FSCE6SSMCheckIsNumeric /@ {
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

            (* CE6SSM input parameters *)
            ,
            OptionValue[TanBeta],
            OptionValue[m0SqGuess],
            OptionValue[m12Guess],
            OptionValue[AzeroGuess],
            OptionValue[LambdaInput],
            OptionValue[KappaInput],
            OptionValue[MuPrimeInput],
            OptionValue[BMuPrimeInput],
            OptionValue[vsInput],
            OptionValue[Lambda12Input]
            ,
            OptionValue[minBRtoPrint],
            OptionValue[maxHigherOrderCorrections],
            OptionValue[alphaThomson],
            OptionValue[offShellVV]
        }
];

Options[FSCE6SSMSet] = Options[FSCE6SSMOpenHandle];

FSCE6SSMSet[handle_Integer, a___, (fsSettings | fsSMParameters | fsModelParameters) -> s_List, r___] :=
    FSCE6SSMSet[handle, a, Sequence @@ s, r];

FSCE6SSMSet[handle_Integer, p:OptionsPattern[]] :=
    FSCE6SSMSetLib[
        handle,
        ReleaseHold[Hold[FSCE6SSMCheckIsNumeric /@ {
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

            (* CE6SSM input parameters *)
            ,
            OptionValue[TanBeta],
            OptionValue[m0SqGuess],
            OptionValue[m12Guess],
            OptionValue[AzeroGuess],
            OptionValue[LambdaInput],
            OptionValue[KappaInput],
            OptionValue[MuPrimeInput],
            OptionValue[BMuPrimeInput],
            OptionValue[vsInput],
            OptionValue[Lambda12Input]
            ,
            OptionValue[minBRtoPrint],
            OptionValue[maxHigherOrderCorrections],
            OptionValue[alphaThomson],
            OptionValue[offShellVV]
        }] /. HoldPattern[OptionValue[param_]] :> param /.
        { p } /.
        FSCE6SSMGetSettings[handle] /.
        FSCE6SSMGetSMInputParameters[handle] /.
        FSCE6SSMGetInputParameters[handle]]];
