libE6SSMtower = FileNameJoin[{Directory[], "models", "E6SSMtower", "E6SSMtower_librarylink.so"}];

FSE6SSMtowerGetSettings = LibraryFunctionLoad[libE6SSMtower, "FSE6SSMtowerGetSettings", LinkObject, LinkObject];
FSE6SSMtowerGetSMInputParameters = LibraryFunctionLoad[libE6SSMtower, "FSE6SSMtowerGetSMInputParameters", LinkObject, LinkObject];
FSE6SSMtowerGetInputParameters = LibraryFunctionLoad[libE6SSMtower, "FSE6SSMtowerGetInputParameters", LinkObject, LinkObject];

FSE6SSMtowerOpenHandleLib = LibraryFunctionLoad[libE6SSMtower, "FSE6SSMtowerOpenHandle", {{Real,1}}, Integer];
FSE6SSMtowerCloseHandle = LibraryFunctionLoad[libE6SSMtower, "FSE6SSMtowerCloseHandle", {Integer}, Void];

FSE6SSMtowerSetLib = LibraryFunctionLoad[libE6SSMtower, "FSE6SSMtowerSet", {Integer, {Real,1}}, Void];

FSE6SSMtowerCalculateSpectrum = LibraryFunctionLoad[libE6SSMtower, "FSE6SSMtowerCalculateSpectrum", LinkObject, LinkObject];
FSE6SSMtowerCalculateObservables = LibraryFunctionLoad[libE6SSMtower, "FSE6SSMtowerCalculateObservables", LinkObject, LinkObject];

FSE6SSMtowerCalculateSpectrum::error = "`1`";
FSE6SSMtowerCalculateSpectrum::warning = "`1`";

FSE6SSMtowerCalculateObservables::error = "`1`";
FSE6SSMtowerCalculateObservables::warning = "`1`";

FSE6SSMtower::info = "`1`";
FSE6SSMtower::nonum = "Error: `1` is not a numeric input value!";
FSE6SSMtowerMessage[s_] := Message[FSE6SSMtower::info, s];

FSE6SSMtowerCheckIsNumeric[a_?NumericQ] := a;
FSE6SSMtowerCheckIsNumeric[a_] := (Message[FSE6SSMtower::nonum, a]; Abort[]);

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

fsE6SSMtowerDefaultInputParameters = {
   MSUSY -> 0,
   M1Input -> 0,
   M2Input -> 0,
   M3Input -> 0,
   MuInput -> 0,
   mAInput -> 0,
   TanBeta -> 0,
   LambdaInput -> 0,
   gNInput -> 0,
   M4Input -> 0,
   mHp2Input -> 0,
   mHpbar2Input -> 0,
   MuPrInput -> 0,
   BMuPrInput -> 0,
   mq2Input -> {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}},
   mu2Input -> {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}},
   md2Input -> {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}},
   ml2Input -> {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}},
   me2Input -> {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}},
   AuInput -> {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}},
   AdInput -> {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}},
   AeInput -> {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}},
   Lambda12Input -> {{0, 0}, {0, 0}},
   ALambda12Input -> {{0, 0}, {0, 0}},
   KappaInput -> {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}},
   AKappaInput -> {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}},
   mDx2Input -> {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}},
   mDxbar2Input -> {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}},
   mH1I2Input -> {{0, 0}, {0, 0}},
   mH2I2Input -> {{0, 0}, {0, 0}},
   msI2Input -> {{0, 0}, {0, 0}}
};

Options[FSE6SSMtowerOpenHandle] = {
    Sequence @@ fsDefaultSettings,
    Sequence @@ fsDefaultSMParameters,
    Sequence @@ fsE6SSMtowerDefaultInputParameters
};

FSE6SSMtowerOpenHandle[a___, (fsSettings | fsSMParameters | fsModelParameters) -> s_List, r___] :=
    FSE6SSMtowerOpenHandle[a, Sequence @@ s, r];

FSE6SSMtowerOpenHandle[OptionsPattern[]] :=
    FSE6SSMtowerOpenHandleLib[
        FSE6SSMtowerCheckIsNumeric /@ {
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

            (* E6SSMtower input parameters *)
            ,
            OptionValue[MSUSY],
            OptionValue[M1Input],
            OptionValue[M2Input],
            OptionValue[M3Input],
            OptionValue[MuInput],
            OptionValue[mAInput],
            OptionValue[TanBeta],
            OptionValue[LambdaInput],
            OptionValue[gNInput],
            OptionValue[M4Input],
            OptionValue[mHp2Input],
            OptionValue[mHpbar2Input],
            OptionValue[MuPrInput],
            OptionValue[BMuPrInput],
            OptionValue[mq2Input][[1,1]],
            OptionValue[mq2Input][[1,2]],
            OptionValue[mq2Input][[1,3]],
            OptionValue[mq2Input][[2,1]],
            OptionValue[mq2Input][[2,2]],
            OptionValue[mq2Input][[2,3]],
            OptionValue[mq2Input][[3,1]],
            OptionValue[mq2Input][[3,2]],
            OptionValue[mq2Input][[3,3]],
            OptionValue[mu2Input][[1,1]],
            OptionValue[mu2Input][[1,2]],
            OptionValue[mu2Input][[1,3]],
            OptionValue[mu2Input][[2,1]],
            OptionValue[mu2Input][[2,2]],
            OptionValue[mu2Input][[2,3]],
            OptionValue[mu2Input][[3,1]],
            OptionValue[mu2Input][[3,2]],
            OptionValue[mu2Input][[3,3]],
            OptionValue[md2Input][[1,1]],
            OptionValue[md2Input][[1,2]],
            OptionValue[md2Input][[1,3]],
            OptionValue[md2Input][[2,1]],
            OptionValue[md2Input][[2,2]],
            OptionValue[md2Input][[2,3]],
            OptionValue[md2Input][[3,1]],
            OptionValue[md2Input][[3,2]],
            OptionValue[md2Input][[3,3]],
            OptionValue[ml2Input][[1,1]],
            OptionValue[ml2Input][[1,2]],
            OptionValue[ml2Input][[1,3]],
            OptionValue[ml2Input][[2,1]],
            OptionValue[ml2Input][[2,2]],
            OptionValue[ml2Input][[2,3]],
            OptionValue[ml2Input][[3,1]],
            OptionValue[ml2Input][[3,2]],
            OptionValue[ml2Input][[3,3]],
            OptionValue[me2Input][[1,1]],
            OptionValue[me2Input][[1,2]],
            OptionValue[me2Input][[1,3]],
            OptionValue[me2Input][[2,1]],
            OptionValue[me2Input][[2,2]],
            OptionValue[me2Input][[2,3]],
            OptionValue[me2Input][[3,1]],
            OptionValue[me2Input][[3,2]],
            OptionValue[me2Input][[3,3]],
            OptionValue[AuInput][[1,1]],
            OptionValue[AuInput][[1,2]],
            OptionValue[AuInput][[1,3]],
            OptionValue[AuInput][[2,1]],
            OptionValue[AuInput][[2,2]],
            OptionValue[AuInput][[2,3]],
            OptionValue[AuInput][[3,1]],
            OptionValue[AuInput][[3,2]],
            OptionValue[AuInput][[3,3]],
            OptionValue[AdInput][[1,1]],
            OptionValue[AdInput][[1,2]],
            OptionValue[AdInput][[1,3]],
            OptionValue[AdInput][[2,1]],
            OptionValue[AdInput][[2,2]],
            OptionValue[AdInput][[2,3]],
            OptionValue[AdInput][[3,1]],
            OptionValue[AdInput][[3,2]],
            OptionValue[AdInput][[3,3]],
            OptionValue[AeInput][[1,1]],
            OptionValue[AeInput][[1,2]],
            OptionValue[AeInput][[1,3]],
            OptionValue[AeInput][[2,1]],
            OptionValue[AeInput][[2,2]],
            OptionValue[AeInput][[2,3]],
            OptionValue[AeInput][[3,1]],
            OptionValue[AeInput][[3,2]],
            OptionValue[AeInput][[3,3]],
            OptionValue[Lambda12Input][[1,1]],
            OptionValue[Lambda12Input][[1,2]],
            OptionValue[Lambda12Input][[2,1]],
            OptionValue[Lambda12Input][[2,2]],
            OptionValue[ALambda12Input][[1,1]],
            OptionValue[ALambda12Input][[1,2]],
            OptionValue[ALambda12Input][[2,1]],
            OptionValue[ALambda12Input][[2,2]],
            OptionValue[KappaInput][[1,1]],
            OptionValue[KappaInput][[1,2]],
            OptionValue[KappaInput][[1,3]],
            OptionValue[KappaInput][[2,1]],
            OptionValue[KappaInput][[2,2]],
            OptionValue[KappaInput][[2,3]],
            OptionValue[KappaInput][[3,1]],
            OptionValue[KappaInput][[3,2]],
            OptionValue[KappaInput][[3,3]],
            OptionValue[AKappaInput][[1,1]],
            OptionValue[AKappaInput][[1,2]],
            OptionValue[AKappaInput][[1,3]],
            OptionValue[AKappaInput][[2,1]],
            OptionValue[AKappaInput][[2,2]],
            OptionValue[AKappaInput][[2,3]],
            OptionValue[AKappaInput][[3,1]],
            OptionValue[AKappaInput][[3,2]],
            OptionValue[AKappaInput][[3,3]],
            OptionValue[mDx2Input][[1,1]],
            OptionValue[mDx2Input][[1,2]],
            OptionValue[mDx2Input][[1,3]],
            OptionValue[mDx2Input][[2,1]],
            OptionValue[mDx2Input][[2,2]],
            OptionValue[mDx2Input][[2,3]],
            OptionValue[mDx2Input][[3,1]],
            OptionValue[mDx2Input][[3,2]],
            OptionValue[mDx2Input][[3,3]],
            OptionValue[mDxbar2Input][[1,1]],
            OptionValue[mDxbar2Input][[1,2]],
            OptionValue[mDxbar2Input][[1,3]],
            OptionValue[mDxbar2Input][[2,1]],
            OptionValue[mDxbar2Input][[2,2]],
            OptionValue[mDxbar2Input][[2,3]],
            OptionValue[mDxbar2Input][[3,1]],
            OptionValue[mDxbar2Input][[3,2]],
            OptionValue[mDxbar2Input][[3,3]],
            OptionValue[mH1I2Input][[1,1]],
            OptionValue[mH1I2Input][[1,2]],
            OptionValue[mH1I2Input][[2,1]],
            OptionValue[mH1I2Input][[2,2]],
            OptionValue[mH2I2Input][[1,1]],
            OptionValue[mH2I2Input][[1,2]],
            OptionValue[mH2I2Input][[2,1]],
            OptionValue[mH2I2Input][[2,2]],
            OptionValue[msI2Input][[1,1]],
            OptionValue[msI2Input][[1,2]],
            OptionValue[msI2Input][[2,1]],
            OptionValue[msI2Input][[2,2]]
        }
];

Options[FSE6SSMtowerSet] = Options[FSE6SSMtowerOpenHandle];

FSE6SSMtowerSet[handle_Integer, a___, (fsSettings | fsSMParameters | fsModelParameters) -> s_List, r___] :=
    FSE6SSMtowerSet[handle, a, Sequence @@ s, r];

FSE6SSMtowerSet[handle_Integer, p:OptionsPattern[]] :=
    FSE6SSMtowerSetLib[
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

            (* E6SSMtower input parameters *)
            ,
            OptionValue[MSUSY],
            OptionValue[M1Input],
            OptionValue[M2Input],
            OptionValue[M3Input],
            OptionValue[MuInput],
            OptionValue[mAInput],
            OptionValue[TanBeta],
            OptionValue[LambdaInput],
            OptionValue[gNInput],
            OptionValue[M4Input],
            OptionValue[mHp2Input],
            OptionValue[mHpbar2Input],
            OptionValue[MuPrInput],
            OptionValue[BMuPrInput],
            OptionValue[mq2Input][[1,1]],
            OptionValue[mq2Input][[1,2]],
            OptionValue[mq2Input][[1,3]],
            OptionValue[mq2Input][[2,1]],
            OptionValue[mq2Input][[2,2]],
            OptionValue[mq2Input][[2,3]],
            OptionValue[mq2Input][[3,1]],
            OptionValue[mq2Input][[3,2]],
            OptionValue[mq2Input][[3,3]],
            OptionValue[mu2Input][[1,1]],
            OptionValue[mu2Input][[1,2]],
            OptionValue[mu2Input][[1,3]],
            OptionValue[mu2Input][[2,1]],
            OptionValue[mu2Input][[2,2]],
            OptionValue[mu2Input][[2,3]],
            OptionValue[mu2Input][[3,1]],
            OptionValue[mu2Input][[3,2]],
            OptionValue[mu2Input][[3,3]],
            OptionValue[md2Input][[1,1]],
            OptionValue[md2Input][[1,2]],
            OptionValue[md2Input][[1,3]],
            OptionValue[md2Input][[2,1]],
            OptionValue[md2Input][[2,2]],
            OptionValue[md2Input][[2,3]],
            OptionValue[md2Input][[3,1]],
            OptionValue[md2Input][[3,2]],
            OptionValue[md2Input][[3,3]],
            OptionValue[ml2Input][[1,1]],
            OptionValue[ml2Input][[1,2]],
            OptionValue[ml2Input][[1,3]],
            OptionValue[ml2Input][[2,1]],
            OptionValue[ml2Input][[2,2]],
            OptionValue[ml2Input][[2,3]],
            OptionValue[ml2Input][[3,1]],
            OptionValue[ml2Input][[3,2]],
            OptionValue[ml2Input][[3,3]],
            OptionValue[me2Input][[1,1]],
            OptionValue[me2Input][[1,2]],
            OptionValue[me2Input][[1,3]],
            OptionValue[me2Input][[2,1]],
            OptionValue[me2Input][[2,2]],
            OptionValue[me2Input][[2,3]],
            OptionValue[me2Input][[3,1]],
            OptionValue[me2Input][[3,2]],
            OptionValue[me2Input][[3,3]],
            OptionValue[AuInput][[1,1]],
            OptionValue[AuInput][[1,2]],
            OptionValue[AuInput][[1,3]],
            OptionValue[AuInput][[2,1]],
            OptionValue[AuInput][[2,2]],
            OptionValue[AuInput][[2,3]],
            OptionValue[AuInput][[3,1]],
            OptionValue[AuInput][[3,2]],
            OptionValue[AuInput][[3,3]],
            OptionValue[AdInput][[1,1]],
            OptionValue[AdInput][[1,2]],
            OptionValue[AdInput][[1,3]],
            OptionValue[AdInput][[2,1]],
            OptionValue[AdInput][[2,2]],
            OptionValue[AdInput][[2,3]],
            OptionValue[AdInput][[3,1]],
            OptionValue[AdInput][[3,2]],
            OptionValue[AdInput][[3,3]],
            OptionValue[AeInput][[1,1]],
            OptionValue[AeInput][[1,2]],
            OptionValue[AeInput][[1,3]],
            OptionValue[AeInput][[2,1]],
            OptionValue[AeInput][[2,2]],
            OptionValue[AeInput][[2,3]],
            OptionValue[AeInput][[3,1]],
            OptionValue[AeInput][[3,2]],
            OptionValue[AeInput][[3,3]],
            OptionValue[Lambda12Input][[1,1]],
            OptionValue[Lambda12Input][[1,2]],
            OptionValue[Lambda12Input][[2,1]],
            OptionValue[Lambda12Input][[2,2]],
            OptionValue[ALambda12Input][[1,1]],
            OptionValue[ALambda12Input][[1,2]],
            OptionValue[ALambda12Input][[2,1]],
            OptionValue[ALambda12Input][[2,2]],
            OptionValue[KappaInput][[1,1]],
            OptionValue[KappaInput][[1,2]],
            OptionValue[KappaInput][[1,3]],
            OptionValue[KappaInput][[2,1]],
            OptionValue[KappaInput][[2,2]],
            OptionValue[KappaInput][[2,3]],
            OptionValue[KappaInput][[3,1]],
            OptionValue[KappaInput][[3,2]],
            OptionValue[KappaInput][[3,3]],
            OptionValue[AKappaInput][[1,1]],
            OptionValue[AKappaInput][[1,2]],
            OptionValue[AKappaInput][[1,3]],
            OptionValue[AKappaInput][[2,1]],
            OptionValue[AKappaInput][[2,2]],
            OptionValue[AKappaInput][[2,3]],
            OptionValue[AKappaInput][[3,1]],
            OptionValue[AKappaInput][[3,2]],
            OptionValue[AKappaInput][[3,3]],
            OptionValue[mDx2Input][[1,1]],
            OptionValue[mDx2Input][[1,2]],
            OptionValue[mDx2Input][[1,3]],
            OptionValue[mDx2Input][[2,1]],
            OptionValue[mDx2Input][[2,2]],
            OptionValue[mDx2Input][[2,3]],
            OptionValue[mDx2Input][[3,1]],
            OptionValue[mDx2Input][[3,2]],
            OptionValue[mDx2Input][[3,3]],
            OptionValue[mDxbar2Input][[1,1]],
            OptionValue[mDxbar2Input][[1,2]],
            OptionValue[mDxbar2Input][[1,3]],
            OptionValue[mDxbar2Input][[2,1]],
            OptionValue[mDxbar2Input][[2,2]],
            OptionValue[mDxbar2Input][[2,3]],
            OptionValue[mDxbar2Input][[3,1]],
            OptionValue[mDxbar2Input][[3,2]],
            OptionValue[mDxbar2Input][[3,3]],
            OptionValue[mH1I2Input][[1,1]],
            OptionValue[mH1I2Input][[1,2]],
            OptionValue[mH1I2Input][[2,1]],
            OptionValue[mH1I2Input][[2,2]],
            OptionValue[mH2I2Input][[1,1]],
            OptionValue[mH2I2Input][[1,2]],
            OptionValue[mH2I2Input][[2,1]],
            OptionValue[mH2I2Input][[2,2]],
            OptionValue[msI2Input][[1,1]],
            OptionValue[msI2Input][[1,2]],
            OptionValue[msI2Input][[2,1]],
            OptionValue[msI2Input][[2,2]]
        }] /. HoldPattern[OptionValue[param_]] :> FSE6SSMtowerCheckIsNumeric[param] /.
        { p } /.
        FSE6SSMtowerGetSettings[handle] /.
        FSE6SSMtowerGetSMInputParameters[handle] /.
        FSE6SSMtowerGetInputParameters[handle]]];
