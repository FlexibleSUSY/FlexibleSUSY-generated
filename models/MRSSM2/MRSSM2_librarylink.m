Print["================================"];
Print["FlexibleSUSY 2.8.0"];
Print["MRSSM2"];
Print["http://flexiblesusy.hepforge.org"];
Print["================================"];

libMRSSM2 = FileNameJoin[{Directory[], "models", "MRSSM2", "MRSSM2_librarylink.so"}];

FSMRSSM2GetSettings = LibraryFunctionLoad[libMRSSM2, "FSMRSSM2GetSettings", LinkObject, LinkObject];
FSMRSSM2GetSMInputParameters = LibraryFunctionLoad[libMRSSM2, "FSMRSSM2GetSMInputParameters", LinkObject, LinkObject];
FSMRSSM2GetInputParameters = LibraryFunctionLoad[libMRSSM2, "FSMRSSM2GetInputParameters", LinkObject, LinkObject];
FSMRSSM2GetProblems = LibraryFunctionLoad[libMRSSM2, "FSMRSSM2GetProblems", LinkObject, LinkObject];
FSMRSSM2GetWarnings = LibraryFunctionLoad[libMRSSM2, "FSMRSSM2GetWarnings", LinkObject, LinkObject];
FSMRSSM2ToSLHA = LibraryFunctionLoad[libMRSSM2, "FSMRSSM2ToSLHA", LinkObject, LinkObject];

FSMRSSM2OpenHandleLib = LibraryFunctionLoad[libMRSSM2, "FSMRSSM2OpenHandle", {{Real,1}}, Integer];
FSMRSSM2CloseHandle = LibraryFunctionLoad[libMRSSM2, "FSMRSSM2CloseHandle", {Integer}, Void];

FSMRSSM2SetLib = LibraryFunctionLoad[libMRSSM2, "FSMRSSM2Set", {Integer, {Real,1}}, Void];

FSMRSSM2CalculateSpectrum = LibraryFunctionLoad[libMRSSM2, "FSMRSSM2CalculateSpectrum", LinkObject, LinkObject];
FSMRSSM2CalculateObservables = LibraryFunctionLoad[libMRSSM2, "FSMRSSM2CalculateObservables", LinkObject, LinkObject];


FSMRSSM2CalculateSpectrum::error = "`1`";
FSMRSSM2CalculateSpectrum::warning = "`1`";

FSMRSSM2CalculateObservables::error = "`1`";
FSMRSSM2CalculateObservables::warning = "`1`";


FSMRSSM2::info = "`1`";
FSMRSSM2::nonum = "Error: `1` is not a numeric input value!";
FSMRSSM2Message[s_] := Message[FSMRSSM2::info, s];

FSMRSSM2CheckIsNumeric[a_?NumericQ] := a;
FSMRSSM2CheckIsNumeric[a_] := (Message[FSMRSSM2::nonum, a]; Abort[]);

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

fsMRSSM2DefaultInputParameters = {
   TanBeta -> 0,
   Ms -> 0,
   LamTDInput -> 0,
   LamTUInput -> 0,
   LamSDInput -> 0,
   LamSUInput -> 0,
   MuDInput -> 0,
   MuUInput -> 0,
   BMuInput -> 0,
   mq2Input -> {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}},
   ml2Input -> {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}},
   md2Input -> {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}},
   mu2Input -> {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}},
   me2Input -> {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}},
   mS2Input -> 0,
   mT2Input -> 0,
   moc2Input -> 0,
   mRd2Input -> 0,
   mRu2Input -> 0,
   MDBSInput -> 0,
   MDWBTInput -> 0,
   MDGocInput -> 0
};

Options[FSMRSSM2OpenHandle] = {
    Sequence @@ fsDefaultSettings,
    Sequence @@ fsDefaultSMParameters,
    Sequence @@ fsMRSSM2DefaultInputParameters

};

FSMRSSM2OpenHandle[a___, (fsSettings | fsSMParameters | fsModelParameters ) -> s_List, r___] :=
    FSMRSSM2OpenHandle[a, Sequence @@ s, r];

FSMRSSM2OpenHandle[OptionsPattern[]] :=
    FSMRSSM2OpenHandleLib[
        FSMRSSM2CheckIsNumeric /@ {
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

            (* MRSSM2 input parameters *)
            ,
            OptionValue[TanBeta],
            OptionValue[Ms],
            OptionValue[LamTDInput],
            OptionValue[LamTUInput],
            OptionValue[LamSDInput],
            OptionValue[LamSUInput],
            OptionValue[MuDInput],
            OptionValue[MuUInput],
            OptionValue[BMuInput],
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
            OptionValue[mS2Input],
            OptionValue[mT2Input],
            OptionValue[moc2Input],
            OptionValue[mRd2Input],
            OptionValue[mRu2Input],
            OptionValue[MDBSInput],
            OptionValue[MDWBTInput],
            OptionValue[MDGocInput]

        }
];

Options[FSMRSSM2Set] = Options[FSMRSSM2OpenHandle];

FSMRSSM2Set[handle_Integer, a___, (fsSettings | fsSMParameters | fsModelParameters) -> s_List, r___] :=
    FSMRSSM2Set[handle, a, Sequence @@ s, r];

FSMRSSM2Set[handle_Integer, p:OptionsPattern[]] :=
    FSMRSSM2SetLib[
        handle,
        ReleaseHold[Hold[FSMRSSM2CheckIsNumeric /@ {
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

            (* MRSSM2 input parameters *)
            ,
            OptionValue[TanBeta],
            OptionValue[Ms],
            OptionValue[LamTDInput],
            OptionValue[LamTUInput],
            OptionValue[LamSDInput],
            OptionValue[LamSUInput],
            OptionValue[MuDInput],
            OptionValue[MuUInput],
            OptionValue[BMuInput],
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
            OptionValue[mS2Input],
            OptionValue[mT2Input],
            OptionValue[moc2Input],
            OptionValue[mRd2Input],
            OptionValue[mRu2Input],
            OptionValue[MDBSInput],
            OptionValue[MDWBTInput],
            OptionValue[MDGocInput]

        }] /. HoldPattern[OptionValue[param_]] :> param /.
        { p } /.
        FSMRSSM2GetSettings[handle] /.
        FSMRSSM2GetSMInputParameters[handle] /.
        FSMRSSM2GetInputParameters[handle]]];