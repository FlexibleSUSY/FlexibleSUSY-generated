libSMSSM = FileNameJoin[{Directory[], "models", "SMSSM", "SMSSM_librarylink.so"}];

FSSMSSMGetSettings = LibraryFunctionLoad[libSMSSM, "FSSMSSMGetSettings", LinkObject, LinkObject];
FSSMSSMGetSMInputParameters = LibraryFunctionLoad[libSMSSM, "FSSMSSMGetSMInputParameters", LinkObject, LinkObject];
FSSMSSMGetInputParameters = LibraryFunctionLoad[libSMSSM, "FSSMSSMGetInputParameters", LinkObject, LinkObject];

FSSMSSMOpenHandleLib = LibraryFunctionLoad[libSMSSM, "FSSMSSMOpenHandle", {{Real,1}}, Integer];
FSSMSSMCloseHandle = LibraryFunctionLoad[libSMSSM, "FSSMSSMCloseHandle", {Integer}, Void];

FSSMSSMSetLib = LibraryFunctionLoad[libSMSSM, "FSSMSSMSet", {Integer, {Real,1}}, Void];

FSSMSSMCalculateSpectrum = LibraryFunctionLoad[libSMSSM, "FSSMSSMCalculateSpectrum", LinkObject, LinkObject];
FSSMSSMCalculateObservables = LibraryFunctionLoad[libSMSSM, "FSSMSSMCalculateObservables", LinkObject, LinkObject];

FSSMSSMCalculateSpectrum::error = "`1`";
FSSMSSMCalculateSpectrum::warning = "`1`";

FSSMSSMCalculateObservables::error = "`1`";
FSSMSSMCalculateObservables::warning = "`1`";

FSSMSSM::info = "`1`";
FSSMSSMMessage[s_] := Message[FSSMSSM::info, s];

fsDefaultSettings = {
      precisionGoal -> 1.*^-4,           (* FlexibleSUSY[0] *)
      maxIterations -> 0,                (* FlexibleSUSY[1] *)
      calculateStandardModelMasses -> 0, (* FlexibleSUSY[3] *)
      poleMassLoopOrder -> 2,            (* FlexibleSUSY[4] *)
      ewsbLoopOrder -> 2,                (* FlexibleSUSY[5] *)
      betaFunctionLoopOrder -> 3,        (* FlexibleSUSY[6] *)
      thresholdCorrectionsLoopOrder -> 2,(* FlexibleSUSY[7] *)
      higgs2loopCorrectionAtAs -> 1,     (* FlexibleSUSY[8] *)
      higgs2loopCorrectionAbAs -> 1,     (* FlexibleSUSY[9] *)
      higgs2loopCorrectionAtAt -> 1,     (* FlexibleSUSY[10] *)
      higgs2loopCorrectionAtauAtau -> 1, (* FlexibleSUSY[11] *)
      forceOutput -> 0,                  (* FlexibleSUSY[12] *)
      topPoleQCDCorrections -> 1,        (* FlexibleSUSY[13] *)
      betaZeroThreshold -> 1.*^-11,      (* FlexibleSUSY[14] *)
      forcePositiveMasses -> 0,          (* FlexibleSUSY[16] *)
      poleMassScale -> 0,                (* FlexibleSUSY[17] *)
      eftPoleMassScale -> 0,             (* FlexibleSUSY[18] *)
      eftMatchingScale -> 0,             (* FlexibleSUSY[19] *)
      eftMatchingLoopOrderUp -> 2,       (* FlexibleSUSY[20] *)
      eftMatchingLoopOrderDown -> 1,     (* FlexibleSUSY[21] *)
      eftHiggsIndex -> 0,                (* FlexibleSUSY[22] *)
      calculateBSMMasses -> 1,           (* FlexibleSUSY[23] *)
      parameterOutputScale -> 0          (* MODSEL[12] *)
};

fsDefaultSMParameters = {
    alphaEmMZ -> 1/127.916, (* SMINPUTS[1] *)
    GF -> 1.16637*^-5,      (* SMINPUTS[2] *)
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

fsSMSSMDefaultInputParameters = {
   m0 -> 0,
   m12 -> 0,
   TanBeta -> 0,
   SignMu -> 0,
   Azero -> 0,
   LambdaInput -> 0,
   KappaInput -> 0,
   LambdaSInput -> 0,
   L1Input -> 0,
   MSInput -> 0,
   BMSInput -> 0
};

Options[FSSMSSMOpenHandle] = {
    Sequence @@ fsDefaultSettings,
    Sequence @@ fsDefaultSMParameters,
    Sequence @@ fsSMSSMDefaultInputParameters
};

FSSMSSMOpenHandle[a___, (fsSettings | fsSMParameters | fsModelParameters) -> s_List, r___] :=
    FSSMSSMOpenHandle[a, Sequence @@ s, r];

FSSMSSMOpenHandle[OptionsPattern[]] :=
    FSSMSSMOpenHandleLib[
        {
            (* spectrum generator settings *)
            OptionValue[precisionGoal],
            OptionValue[maxIterations],
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

            (* SMSSM input parameters *)
            ,
            OptionValue[m0],
            OptionValue[m12],
            OptionValue[TanBeta],
            OptionValue[SignMu],
            OptionValue[Azero],
            OptionValue[LambdaInput],
            OptionValue[KappaInput],
            OptionValue[LambdaSInput],
            OptionValue[L1Input],
            OptionValue[MSInput],
            OptionValue[BMSInput]
        }
];

Options[FSSMSSMSet] = Options[FSSMSSMOpenHandle];

FSSMSSMSet[handle_Integer, a___, (fsSettings | fsSMParameters | fsModelParameters) -> s_List, r___] :=
    FSSMSSMSet[handle, a, Sequence @@ s, r];

FSSMSSMSet[handle_Integer, p:OptionsPattern[]] :=
    FSSMSSMSetLib[
        handle,
        ReleaseHold[Hold[{
            (* spectrum generator settings *)
            OptionValue[precisionGoal],
            OptionValue[maxIterations],
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

            (* SMSSM input parameters *)
            ,
            OptionValue[m0],
            OptionValue[m12],
            OptionValue[TanBeta],
            OptionValue[SignMu],
            OptionValue[Azero],
            OptionValue[LambdaInput],
            OptionValue[KappaInput],
            OptionValue[LambdaSInput],
            OptionValue[L1Input],
            OptionValue[MSInput],
            OptionValue[BMSInput]
        }] /. HoldPattern[OptionValue[param_]] :> param /.
        { p } /.
        FSSMSSMGetSettings[handle] /.
        FSSMSSMGetSMInputParameters[handle] /.
        FSSMSSMGetInputParameters[handle]]];