Print["================================"];
Print["FlexibleSUSY 2.6.1"];
Print["NMSSM"];
Print["http://flexiblesusy.hepforge.org"];
Print["================================"];

libNMSSM = FileNameJoin[{Directory[], "models", "NMSSM", "NMSSM_librarylink.so"}];

FSNMSSMGetSettings = LibraryFunctionLoad[libNMSSM, "FSNMSSMGetSettings", LinkObject, LinkObject];
FSNMSSMGetSMInputParameters = LibraryFunctionLoad[libNMSSM, "FSNMSSMGetSMInputParameters", LinkObject, LinkObject];
FSNMSSMGetInputParameters = LibraryFunctionLoad[libNMSSM, "FSNMSSMGetInputParameters", LinkObject, LinkObject];
FSNMSSMGetProblems = LibraryFunctionLoad[libNMSSM, "FSNMSSMGetProblems", LinkObject, LinkObject];
FSNMSSMGetWarnings = LibraryFunctionLoad[libNMSSM, "FSNMSSMGetWarnings", LinkObject, LinkObject];
FSNMSSMToSLHA = LibraryFunctionLoad[libNMSSM, "FSNMSSMToSLHA", LinkObject, LinkObject];

FSNMSSMOpenHandleLib = LibraryFunctionLoad[libNMSSM, "FSNMSSMOpenHandle", {{Real,1}}, Integer];
FSNMSSMCloseHandle = LibraryFunctionLoad[libNMSSM, "FSNMSSMCloseHandle", {Integer}, Void];

FSNMSSMSetLib = LibraryFunctionLoad[libNMSSM, "FSNMSSMSet", {Integer, {Real,1}}, Void];

FSNMSSMCalculateSpectrum = LibraryFunctionLoad[libNMSSM, "FSNMSSMCalculateSpectrum", LinkObject, LinkObject];
FSNMSSMCalculateObservables = LibraryFunctionLoad[libNMSSM, "FSNMSSMCalculateObservables", LinkObject, LinkObject];
FSNMSSMCalculateDecays = LibraryFunctionLoad[libNMSSM, "FSNMSSMCalculateDecays", LinkObject, LinkObject];

FSNMSSMCalculateSpectrum::error = "`1`";
FSNMSSMCalculateSpectrum::warning = "`1`";

FSNMSSMCalculateObservables::error = "`1`";
FSNMSSMCalculateObservables::warning = "`1`";

FSNMSSMCalculateDecays::error = "`1`";
FSNMSSMCalculateDecays::warning = "`1`";

FSNMSSM::info = "`1`";
FSNMSSM::nonum = "Error: `1` is not a numeric input value!";
FSNMSSMMessage[s_] := Message[FSNMSSM::info, s];

FSNMSSMCheckIsNumeric[a_?NumericQ] := a;
FSNMSSMCheckIsNumeric[a_] := (Message[FSNMSSM::nonum, a]; Abort[]);

fsDefaultSettings = {
      precisionGoal -> 1.*^-4,           (* FlexibleSUSY[0] *)
      maxIterations -> 0,                (* FlexibleSUSY[1] *)
      solver -> 1,     (* FlexibleSUSY[2] *)
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

fsNMSSMDefaultInputParameters = {
   m0 -> 0,
   m12 -> 0,
   TanBeta -> 0,
   SignvS -> 0,
   Azero -> 0,
   LambdaInput -> 0
};

Options[FSNMSSMOpenHandle] = {
    Sequence @@ fsDefaultSettings,
    Sequence @@ fsDefaultSMParameters,
    Sequence @@ fsNMSSMDefaultInputParameters
   , Sequence @@ fdDefaultSettings
};

FSNMSSMOpenHandle[a___, (fsSettings | fsSMParameters | fsModelParameters | fdSettings) -> s_List, r___] :=
    FSNMSSMOpenHandle[a, Sequence @@ s, r];

FSNMSSMOpenHandle[OptionsPattern[]] :=
    FSNMSSMOpenHandleLib[
        FSNMSSMCheckIsNumeric /@ {
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

            (* NMSSM input parameters *)
            ,
            OptionValue[m0],
            OptionValue[m12],
            OptionValue[TanBeta],
            OptionValue[SignvS],
            OptionValue[Azero],
            OptionValue[LambdaInput]
            ,
            OptionValue[minBRtoPrint],
            OptionValue[maxHigherOrderCorrections],
            OptionValue[alphaThomson],
            OptionValue[offShellVV]
        }
];

Options[FSNMSSMSet] = Options[FSNMSSMOpenHandle];

FSNMSSMSet[handle_Integer, a___, (fsSettings | fsSMParameters | fsModelParameters) -> s_List, r___] :=
    FSNMSSMSet[handle, a, Sequence @@ s, r];

FSNMSSMSet[handle_Integer, p:OptionsPattern[]] :=
    FSNMSSMSetLib[
        handle,
        ReleaseHold[Hold[FSNMSSMCheckIsNumeric /@ {
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

            (* NMSSM input parameters *)
            ,
            OptionValue[m0],
            OptionValue[m12],
            OptionValue[TanBeta],
            OptionValue[SignvS],
            OptionValue[Azero],
            OptionValue[LambdaInput]
            ,
            OptionValue[minBRtoPrint],
            OptionValue[maxHigherOrderCorrections],
            OptionValue[alphaThomson],
            OptionValue[offShellVV]
        }] /. HoldPattern[OptionValue[param_]] :> param /.
        { p } /.
        FSNMSSMGetSettings[handle] /.
        FSNMSSMGetSMInputParameters[handle] /.
        FSNMSSMGetInputParameters[handle]]];
