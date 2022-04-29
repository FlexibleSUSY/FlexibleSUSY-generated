Print["================================"];
Print["FlexibleSUSY 2.6.2"];
Print["E6SSMEFTHiggs"];
Print["http://flexiblesusy.hepforge.org"];
Print["================================"];

libE6SSMEFTHiggs = FileNameJoin[{Directory[], "models", "E6SSMEFTHiggs", "E6SSMEFTHiggs_librarylink.so"}];

FSE6SSMEFTHiggsGetSettings = LibraryFunctionLoad[libE6SSMEFTHiggs, "FSE6SSMEFTHiggsGetSettings", LinkObject, LinkObject];
FSE6SSMEFTHiggsGetSMInputParameters = LibraryFunctionLoad[libE6SSMEFTHiggs, "FSE6SSMEFTHiggsGetSMInputParameters", LinkObject, LinkObject];
FSE6SSMEFTHiggsGetInputParameters = LibraryFunctionLoad[libE6SSMEFTHiggs, "FSE6SSMEFTHiggsGetInputParameters", LinkObject, LinkObject];
FSE6SSMEFTHiggsGetProblems = LibraryFunctionLoad[libE6SSMEFTHiggs, "FSE6SSMEFTHiggsGetProblems", LinkObject, LinkObject];
FSE6SSMEFTHiggsGetWarnings = LibraryFunctionLoad[libE6SSMEFTHiggs, "FSE6SSMEFTHiggsGetWarnings", LinkObject, LinkObject];
FSE6SSMEFTHiggsToSLHA = LibraryFunctionLoad[libE6SSMEFTHiggs, "FSE6SSMEFTHiggsToSLHA", LinkObject, LinkObject];

FSE6SSMEFTHiggsOpenHandleLib = LibraryFunctionLoad[libE6SSMEFTHiggs, "FSE6SSMEFTHiggsOpenHandle", {{Real,1}}, Integer];
FSE6SSMEFTHiggsCloseHandle = LibraryFunctionLoad[libE6SSMEFTHiggs, "FSE6SSMEFTHiggsCloseHandle", {Integer}, Void];

FSE6SSMEFTHiggsSetLib = LibraryFunctionLoad[libE6SSMEFTHiggs, "FSE6SSMEFTHiggsSet", {Integer, {Real,1}}, Void];

FSE6SSMEFTHiggsCalculateSpectrum = LibraryFunctionLoad[libE6SSMEFTHiggs, "FSE6SSMEFTHiggsCalculateSpectrum", LinkObject, LinkObject];
FSE6SSMEFTHiggsCalculateObservables = LibraryFunctionLoad[libE6SSMEFTHiggs, "FSE6SSMEFTHiggsCalculateObservables", LinkObject, LinkObject];
FSE6SSMEFTHiggsCalculateDecays = LibraryFunctionLoad[libE6SSMEFTHiggs, "FSE6SSMEFTHiggsCalculateDecays", LinkObject, LinkObject];

FSE6SSMEFTHiggsCalculateSpectrum::error = "`1`";
FSE6SSMEFTHiggsCalculateSpectrum::warning = "`1`";

FSE6SSMEFTHiggsCalculateObservables::error = "`1`";
FSE6SSMEFTHiggsCalculateObservables::warning = "`1`";

FSE6SSMEFTHiggsCalculateDecays::error = "`1`";
FSE6SSMEFTHiggsCalculateDecays::warning = "`1`";

FSE6SSMEFTHiggs::info = "`1`";
FSE6SSMEFTHiggs::nonum = "Error: `1` is not a numeric input value!";
FSE6SSMEFTHiggsMessage[s_] := Message[FSE6SSMEFTHiggs::info, s];

FSE6SSMEFTHiggsCheckIsNumeric[a_?NumericQ] := a;
FSE6SSMEFTHiggsCheckIsNumeric[a_] := (Message[FSE6SSMEFTHiggs::nonum, a]; Abort[]);

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

fsE6SSMEFTHiggsDefaultInputParameters = {
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
   AeInput -> {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}},
   AdInput -> {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}},
   AuInput -> {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}},
   ml2Input -> {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}},
   me2Input -> {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}},
   mq2Input -> {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}},
   md2Input -> {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}},
   mu2Input -> {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}},
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

Options[FSE6SSMEFTHiggsOpenHandle] = {
    Sequence @@ fsDefaultSettings,
    Sequence @@ fsDefaultSMParameters,
    Sequence @@ fsE6SSMEFTHiggsDefaultInputParameters
   , Sequence @@ fdDefaultSettings
};

FSE6SSMEFTHiggsOpenHandle[a___, (fsSettings | fsSMParameters | fsModelParameters | fdSettings) -> s_List, r___] :=
    FSE6SSMEFTHiggsOpenHandle[a, Sequence @@ s, r];

FSE6SSMEFTHiggsOpenHandle[OptionsPattern[]] :=
    FSE6SSMEFTHiggsOpenHandleLib[
        FSE6SSMEFTHiggsCheckIsNumeric /@ {
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

            (* E6SSMEFTHiggs input parameters *)
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
            OptionValue[AeInput][[1,1]],
            OptionValue[AeInput][[1,2]],
            OptionValue[AeInput][[1,3]],
            OptionValue[AeInput][[2,1]],
            OptionValue[AeInput][[2,2]],
            OptionValue[AeInput][[2,3]],
            OptionValue[AeInput][[3,1]],
            OptionValue[AeInput][[3,2]],
            OptionValue[AeInput][[3,3]],
            OptionValue[AdInput][[1,1]],
            OptionValue[AdInput][[1,2]],
            OptionValue[AdInput][[1,3]],
            OptionValue[AdInput][[2,1]],
            OptionValue[AdInput][[2,2]],
            OptionValue[AdInput][[2,3]],
            OptionValue[AdInput][[3,1]],
            OptionValue[AdInput][[3,2]],
            OptionValue[AdInput][[3,3]],
            OptionValue[AuInput][[1,1]],
            OptionValue[AuInput][[1,2]],
            OptionValue[AuInput][[1,3]],
            OptionValue[AuInput][[2,1]],
            OptionValue[AuInput][[2,2]],
            OptionValue[AuInput][[2,3]],
            OptionValue[AuInput][[3,1]],
            OptionValue[AuInput][[3,2]],
            OptionValue[AuInput][[3,3]],
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
            OptionValue[mq2Input][[1,1]],
            OptionValue[mq2Input][[1,2]],
            OptionValue[mq2Input][[1,3]],
            OptionValue[mq2Input][[2,1]],
            OptionValue[mq2Input][[2,2]],
            OptionValue[mq2Input][[2,3]],
            OptionValue[mq2Input][[3,1]],
            OptionValue[mq2Input][[3,2]],
            OptionValue[mq2Input][[3,3]],
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
            ,
            OptionValue[minBRtoPrint],
            OptionValue[maxHigherOrderCorrections],
            OptionValue[alphaThomson],
            OptionValue[offShellVV]
        }
];

Options[FSE6SSMEFTHiggsSet] = Options[FSE6SSMEFTHiggsOpenHandle];

FSE6SSMEFTHiggsSet[handle_Integer, a___, (fsSettings | fsSMParameters | fsModelParameters) -> s_List, r___] :=
    FSE6SSMEFTHiggsSet[handle, a, Sequence @@ s, r];

FSE6SSMEFTHiggsSet[handle_Integer, p:OptionsPattern[]] :=
    FSE6SSMEFTHiggsSetLib[
        handle,
        ReleaseHold[Hold[FSE6SSMEFTHiggsCheckIsNumeric /@ {
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

            (* E6SSMEFTHiggs input parameters *)
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
            OptionValue[AeInput][[1,1]],
            OptionValue[AeInput][[1,2]],
            OptionValue[AeInput][[1,3]],
            OptionValue[AeInput][[2,1]],
            OptionValue[AeInput][[2,2]],
            OptionValue[AeInput][[2,3]],
            OptionValue[AeInput][[3,1]],
            OptionValue[AeInput][[3,2]],
            OptionValue[AeInput][[3,3]],
            OptionValue[AdInput][[1,1]],
            OptionValue[AdInput][[1,2]],
            OptionValue[AdInput][[1,3]],
            OptionValue[AdInput][[2,1]],
            OptionValue[AdInput][[2,2]],
            OptionValue[AdInput][[2,3]],
            OptionValue[AdInput][[3,1]],
            OptionValue[AdInput][[3,2]],
            OptionValue[AdInput][[3,3]],
            OptionValue[AuInput][[1,1]],
            OptionValue[AuInput][[1,2]],
            OptionValue[AuInput][[1,3]],
            OptionValue[AuInput][[2,1]],
            OptionValue[AuInput][[2,2]],
            OptionValue[AuInput][[2,3]],
            OptionValue[AuInput][[3,1]],
            OptionValue[AuInput][[3,2]],
            OptionValue[AuInput][[3,3]],
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
            OptionValue[mq2Input][[1,1]],
            OptionValue[mq2Input][[1,2]],
            OptionValue[mq2Input][[1,3]],
            OptionValue[mq2Input][[2,1]],
            OptionValue[mq2Input][[2,2]],
            OptionValue[mq2Input][[2,3]],
            OptionValue[mq2Input][[3,1]],
            OptionValue[mq2Input][[3,2]],
            OptionValue[mq2Input][[3,3]],
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
            ,
            OptionValue[minBRtoPrint],
            OptionValue[maxHigherOrderCorrections],
            OptionValue[alphaThomson],
            OptionValue[offShellVV]
        }] /. HoldPattern[OptionValue[param_]] :> param /.
        { p } /.
        FSE6SSMEFTHiggsGetSettings[handle] /.
        FSE6SSMEFTHiggsGetSMInputParameters[handle] /.
        FSE6SSMEFTHiggsGetInputParameters[handle]]];
