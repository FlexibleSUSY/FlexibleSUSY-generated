libMRSSM = FileNameJoin[{Directory[], "models", "MRSSM", "MRSSM_librarylink.so"}];

FSMRSSMGetSettings = LibraryFunctionLoad[libMRSSM, "FSMRSSMGetSettings", LinkObject, LinkObject];
FSMRSSMGetSMInputParameters = LibraryFunctionLoad[libMRSSM, "FSMRSSMGetSMInputParameters", LinkObject, LinkObject];
FSMRSSMGetInputParameters = LibraryFunctionLoad[libMRSSM, "FSMRSSMGetInputParameters", LinkObject, LinkObject];

FSMRSSMOpenHandleLib = LibraryFunctionLoad[libMRSSM, "FSMRSSMOpenHandle", {{Real,1}}, Integer];
FSMRSSMCloseHandle = LibraryFunctionLoad[libMRSSM, "FSMRSSMCloseHandle", {Integer}, Void];

FSMRSSMSetLib = LibraryFunctionLoad[libMRSSM, "FSMRSSMSet", {Integer, {Real,1}}, Void];

FSMRSSMCalculateSpectrum = LibraryFunctionLoad[libMRSSM, "FSMRSSMCalculateSpectrum", LinkObject, LinkObject];
FSMRSSMCalculateObservables = LibraryFunctionLoad[libMRSSM, "FSMRSSMCalculateObservables", LinkObject, LinkObject];

FSMRSSMCalculateSpectrum::error = "`1`";
FSMRSSMCalculateSpectrum::warning = "`1`";

FSMRSSMCalculateObservables::error = "`1`";
FSMRSSMCalculateObservables::warning = "`1`";

FSMRSSM::info = "`1`";
FSMRSSMMessage[s_] := Message[FSMRSSM::info, s];

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

fsMRSSMDefaultInputParameters = {
   TanBeta -> 0,
   LamTDInput -> 0,
   LamTUInput -> 0,
   LamSDInput -> 0,
   LamSUInput -> 0,
   MuInput -> 0,
   MuDInput -> 0,
   MuUInput -> 0,
   vTInput -> 0,
   vSInput -> 0,
   BMuInput -> 0,
   BMuDInput -> 0,
   BMuUInput -> 0,
   mq2Input -> {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}},
   ml2Input -> {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}},
   md2Input -> {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}},
   mu2Input -> {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}},
   me2Input -> {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}},
   moc2Input -> 0,
   mRd2Input -> 0,
   mRu2Input -> 0,
   MDBSInput -> 0,
   MDWBTInput -> 0,
   MDGocInput -> 0
};

Options[FSMRSSMOpenHandle] = {
    Sequence @@ fsDefaultSettings,
    Sequence @@ fsDefaultSMParameters,
    Sequence @@ fsMRSSMDefaultInputParameters
};

FSMRSSMOpenHandle[a___, (fsSettings | fsSMParameters | fsModelParameters) -> s_List, r___] :=
    FSMRSSMOpenHandle[a, Sequence @@ s, r];

FSMRSSMOpenHandle[OptionsPattern[]] :=
    FSMRSSMOpenHandleLib[
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

            (* MRSSM input parameters *)
            ,
            OptionValue[TanBeta],
            OptionValue[LamTDInput],
            OptionValue[LamTUInput],
            OptionValue[LamSDInput],
            OptionValue[LamSUInput],
            OptionValue[MuInput],
            OptionValue[MuDInput],
            OptionValue[MuUInput],
            OptionValue[vTInput],
            OptionValue[vSInput],
            OptionValue[BMuInput],
            OptionValue[BMuDInput],
            OptionValue[BMuUInput],
            OptionValue[mq2Input][[1,1]],
            OptionValue[mq2Input][[1,2]],
            OptionValue[mq2Input][[1,3]],
            OptionValue[mq2Input][[2,1]],
            OptionValue[mq2Input][[2,2]],
            OptionValue[mq2Input][[2,3]],
            OptionValue[mq2Input][[3,1]],
            OptionValue[mq2Input][[3,2]],
            OptionValue[mq2Input][[3,3]],
            OptionValue[ml2Input][[1,1]],
            OptionValue[ml2Input][[1,2]],
            OptionValue[ml2Input][[1,3]],
            OptionValue[ml2Input][[2,1]],
            OptionValue[ml2Input][[2,2]],
            OptionValue[ml2Input][[2,3]],
            OptionValue[ml2Input][[3,1]],
            OptionValue[ml2Input][[3,2]],
            OptionValue[ml2Input][[3,3]],
            OptionValue[md2Input][[1,1]],
            OptionValue[md2Input][[1,2]],
            OptionValue[md2Input][[1,3]],
            OptionValue[md2Input][[2,1]],
            OptionValue[md2Input][[2,2]],
            OptionValue[md2Input][[2,3]],
            OptionValue[md2Input][[3,1]],
            OptionValue[md2Input][[3,2]],
            OptionValue[md2Input][[3,3]],
            OptionValue[mu2Input][[1,1]],
            OptionValue[mu2Input][[1,2]],
            OptionValue[mu2Input][[1,3]],
            OptionValue[mu2Input][[2,1]],
            OptionValue[mu2Input][[2,2]],
            OptionValue[mu2Input][[2,3]],
            OptionValue[mu2Input][[3,1]],
            OptionValue[mu2Input][[3,2]],
            OptionValue[mu2Input][[3,3]],
            OptionValue[me2Input][[1,1]],
            OptionValue[me2Input][[1,2]],
            OptionValue[me2Input][[1,3]],
            OptionValue[me2Input][[2,1]],
            OptionValue[me2Input][[2,2]],
            OptionValue[me2Input][[2,3]],
            OptionValue[me2Input][[3,1]],
            OptionValue[me2Input][[3,2]],
            OptionValue[me2Input][[3,3]],
            OptionValue[moc2Input],
            OptionValue[mRd2Input],
            OptionValue[mRu2Input],
            OptionValue[MDBSInput],
            OptionValue[MDWBTInput],
            OptionValue[MDGocInput]
        }
];

Options[FSMRSSMSet] = Options[FSMRSSMOpenHandle];

FSMRSSMSet[handle_Integer, a___, (fsSettings | fsSMParameters | fsModelParameters) -> s_List, r___] :=
    FSMRSSMSet[handle, a, Sequence @@ s, r];

FSMRSSMSet[handle_Integer, p:OptionsPattern[]] :=
    FSMRSSMSetLib[
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

            (* MRSSM input parameters *)
            ,
            OptionValue[TanBeta],
            OptionValue[LamTDInput],
            OptionValue[LamTUInput],
            OptionValue[LamSDInput],
            OptionValue[LamSUInput],
            OptionValue[MuInput],
            OptionValue[MuDInput],
            OptionValue[MuUInput],
            OptionValue[vTInput],
            OptionValue[vSInput],
            OptionValue[BMuInput],
            OptionValue[BMuDInput],
            OptionValue[BMuUInput],
            OptionValue[mq2Input][[1,1]],
            OptionValue[mq2Input][[1,2]],
            OptionValue[mq2Input][[1,3]],
            OptionValue[mq2Input][[2,1]],
            OptionValue[mq2Input][[2,2]],
            OptionValue[mq2Input][[2,3]],
            OptionValue[mq2Input][[3,1]],
            OptionValue[mq2Input][[3,2]],
            OptionValue[mq2Input][[3,3]],
            OptionValue[ml2Input][[1,1]],
            OptionValue[ml2Input][[1,2]],
            OptionValue[ml2Input][[1,3]],
            OptionValue[ml2Input][[2,1]],
            OptionValue[ml2Input][[2,2]],
            OptionValue[ml2Input][[2,3]],
            OptionValue[ml2Input][[3,1]],
            OptionValue[ml2Input][[3,2]],
            OptionValue[ml2Input][[3,3]],
            OptionValue[md2Input][[1,1]],
            OptionValue[md2Input][[1,2]],
            OptionValue[md2Input][[1,3]],
            OptionValue[md2Input][[2,1]],
            OptionValue[md2Input][[2,2]],
            OptionValue[md2Input][[2,3]],
            OptionValue[md2Input][[3,1]],
            OptionValue[md2Input][[3,2]],
            OptionValue[md2Input][[3,3]],
            OptionValue[mu2Input][[1,1]],
            OptionValue[mu2Input][[1,2]],
            OptionValue[mu2Input][[1,3]],
            OptionValue[mu2Input][[2,1]],
            OptionValue[mu2Input][[2,2]],
            OptionValue[mu2Input][[2,3]],
            OptionValue[mu2Input][[3,1]],
            OptionValue[mu2Input][[3,2]],
            OptionValue[mu2Input][[3,3]],
            OptionValue[me2Input][[1,1]],
            OptionValue[me2Input][[1,2]],
            OptionValue[me2Input][[1,3]],
            OptionValue[me2Input][[2,1]],
            OptionValue[me2Input][[2,2]],
            OptionValue[me2Input][[2,3]],
            OptionValue[me2Input][[3,1]],
            OptionValue[me2Input][[3,2]],
            OptionValue[me2Input][[3,3]],
            OptionValue[moc2Input],
            OptionValue[mRd2Input],
            OptionValue[mRu2Input],
            OptionValue[MDBSInput],
            OptionValue[MDWBTInput],
            OptionValue[MDGocInput]
        }] /. HoldPattern[OptionValue[param_]] :> param /.
        { p } /.
        FSMRSSMGetSettings[handle] /.
        FSMRSSMGetSMInputParameters[handle] /.
        FSMRSSMGetInputParameters[handle]]];
