Print["================================"];
Print["FlexibleSUSY 2.8.0"];
Print["THDMIIMSSMBC"];
Print["http://flexiblesusy.hepforge.org"];
Print["================================"];

libTHDMIIMSSMBC = FileNameJoin[{Directory[], "models", "THDMIIMSSMBC", "THDMIIMSSMBC_librarylink.so"}];

FSTHDMIIMSSMBCGetSettings = LibraryFunctionLoad[libTHDMIIMSSMBC, "FSTHDMIIMSSMBCGetSettings", LinkObject, LinkObject];
FSTHDMIIMSSMBCGetSMInputParameters = LibraryFunctionLoad[libTHDMIIMSSMBC, "FSTHDMIIMSSMBCGetSMInputParameters", LinkObject, LinkObject];
FSTHDMIIMSSMBCGetInputParameters = LibraryFunctionLoad[libTHDMIIMSSMBC, "FSTHDMIIMSSMBCGetInputParameters", LinkObject, LinkObject];
FSTHDMIIMSSMBCGetProblems = LibraryFunctionLoad[libTHDMIIMSSMBC, "FSTHDMIIMSSMBCGetProblems", LinkObject, LinkObject];
FSTHDMIIMSSMBCGetWarnings = LibraryFunctionLoad[libTHDMIIMSSMBC, "FSTHDMIIMSSMBCGetWarnings", LinkObject, LinkObject];
FSTHDMIIMSSMBCToSLHA = LibraryFunctionLoad[libTHDMIIMSSMBC, "FSTHDMIIMSSMBCToSLHA", LinkObject, LinkObject];

FSTHDMIIMSSMBCOpenHandleLib = LibraryFunctionLoad[libTHDMIIMSSMBC, "FSTHDMIIMSSMBCOpenHandle", {{Real,1}}, Integer];
FSTHDMIIMSSMBCCloseHandle = LibraryFunctionLoad[libTHDMIIMSSMBC, "FSTHDMIIMSSMBCCloseHandle", {Integer}, Void];

FSTHDMIIMSSMBCSetLib = LibraryFunctionLoad[libTHDMIIMSSMBC, "FSTHDMIIMSSMBCSet", {Integer, {Real,1}}, Void];

FSTHDMIIMSSMBCCalculateSpectrum = LibraryFunctionLoad[libTHDMIIMSSMBC, "FSTHDMIIMSSMBCCalculateSpectrum", LinkObject, LinkObject];
FSTHDMIIMSSMBCCalculateObservables = LibraryFunctionLoad[libTHDMIIMSSMBC, "FSTHDMIIMSSMBCCalculateObservables", LinkObject, LinkObject];


FSTHDMIIMSSMBCCalculateSpectrum::error = "`1`";
FSTHDMIIMSSMBCCalculateSpectrum::warning = "`1`";

FSTHDMIIMSSMBCCalculateObservables::error = "`1`";
FSTHDMIIMSSMBCCalculateObservables::warning = "`1`";


FSTHDMIIMSSMBC::info = "`1`";
FSTHDMIIMSSMBC::nonum = "Error: `1` is not a numeric input value!";
FSTHDMIIMSSMBCMessage[s_] := Message[FSTHDMIIMSSMBC::info, s];

FSTHDMIIMSSMBCCheckIsNumeric[a_?NumericQ] := a;
FSTHDMIIMSSMBCCheckIsNumeric[a_] := (Message[FSTHDMIIMSSMBC::nonum, a]; Abort[]);

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

fsTHDMIIMSSMBCDefaultInputParameters = {
   TanBeta -> 0,
   MSUSY -> 0,
   MEWSB -> 0,
   MuInput -> 0,
   M1Input -> 0,
   M2Input -> 0,
   MAInput -> 0,
   LambdaLoopOrder -> 0,
   AeInput -> {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}},
   AdInput -> {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}},
   AuInput -> {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}},
   mslInput -> {0, 0, 0},
   mseInput -> {0, 0, 0},
   msqInput -> {0, 0, 0},
   msdInput -> {0, 0, 0},
   msuInput -> {0, 0, 0}
};

Options[FSTHDMIIMSSMBCOpenHandle] = {
    Sequence @@ fsDefaultSettings,
    Sequence @@ fsDefaultSMParameters,
    Sequence @@ fsTHDMIIMSSMBCDefaultInputParameters

};

FSTHDMIIMSSMBCOpenHandle[a___, (fsSettings | fsSMParameters | fsModelParameters ) -> s_List, r___] :=
    FSTHDMIIMSSMBCOpenHandle[a, Sequence @@ s, r];

FSTHDMIIMSSMBCOpenHandle[OptionsPattern[]] :=
    FSTHDMIIMSSMBCOpenHandleLib[
        FSTHDMIIMSSMBCCheckIsNumeric /@ {
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

            (* THDMIIMSSMBC input parameters *)
            ,
            OptionValue[TanBeta],
            OptionValue[MSUSY],
            OptionValue[MEWSB],
            OptionValue[MuInput],
            OptionValue[M1Input],
            OptionValue[M2Input],
            OptionValue[MAInput],
            OptionValue[LambdaLoopOrder],
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
            OptionValue[mslInput][[1]],
            OptionValue[mslInput][[2]],
            OptionValue[mslInput][[3]],
            OptionValue[mseInput][[1]],
            OptionValue[mseInput][[2]],
            OptionValue[mseInput][[3]],
            OptionValue[msqInput][[1]],
            OptionValue[msqInput][[2]],
            OptionValue[msqInput][[3]],
            OptionValue[msdInput][[1]],
            OptionValue[msdInput][[2]],
            OptionValue[msdInput][[3]],
            OptionValue[msuInput][[1]],
            OptionValue[msuInput][[2]],
            OptionValue[msuInput][[3]]

        }
];

Options[FSTHDMIIMSSMBCSet] = Options[FSTHDMIIMSSMBCOpenHandle];

FSTHDMIIMSSMBCSet[handle_Integer, a___, (fsSettings | fsSMParameters | fsModelParameters) -> s_List, r___] :=
    FSTHDMIIMSSMBCSet[handle, a, Sequence @@ s, r];

FSTHDMIIMSSMBCSet[handle_Integer, p:OptionsPattern[]] :=
    FSTHDMIIMSSMBCSetLib[
        handle,
        ReleaseHold[Hold[FSTHDMIIMSSMBCCheckIsNumeric /@ {
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

            (* THDMIIMSSMBC input parameters *)
            ,
            OptionValue[TanBeta],
            OptionValue[MSUSY],
            OptionValue[MEWSB],
            OptionValue[MuInput],
            OptionValue[M1Input],
            OptionValue[M2Input],
            OptionValue[MAInput],
            OptionValue[LambdaLoopOrder],
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
            OptionValue[mslInput][[1]],
            OptionValue[mslInput][[2]],
            OptionValue[mslInput][[3]],
            OptionValue[mseInput][[1]],
            OptionValue[mseInput][[2]],
            OptionValue[mseInput][[3]],
            OptionValue[msqInput][[1]],
            OptionValue[msqInput][[2]],
            OptionValue[msqInput][[3]],
            OptionValue[msdInput][[1]],
            OptionValue[msdInput][[2]],
            OptionValue[msdInput][[3]],
            OptionValue[msuInput][[1]],
            OptionValue[msuInput][[2]],
            OptionValue[msuInput][[3]]

        }] /. HoldPattern[OptionValue[param_]] :> param /.
        { p } /.
        FSTHDMIIMSSMBCGetSettings[handle] /.
        FSTHDMIIMSSMBCGetSMInputParameters[handle] /.
        FSTHDMIIMSSMBCGetInputParameters[handle]]];
