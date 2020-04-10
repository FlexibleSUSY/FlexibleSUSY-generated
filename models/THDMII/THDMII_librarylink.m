Print["================================"];
Print["FlexibleSUSY 2.4.2"];
Print["THDMII"];
Print["http://flexiblesusy.hepforge.org"];
Print["================================"];

libTHDMII = FileNameJoin[{Directory[], "models", "THDMII", "THDMII_librarylink.so"}];

FSTHDMIIGetSettings = LibraryFunctionLoad[libTHDMII, "FSTHDMIIGetSettings", LinkObject, LinkObject];
FSTHDMIIGetSMInputParameters = LibraryFunctionLoad[libTHDMII, "FSTHDMIIGetSMInputParameters", LinkObject, LinkObject];
FSTHDMIIGetInputParameters = LibraryFunctionLoad[libTHDMII, "FSTHDMIIGetInputParameters", LinkObject, LinkObject];
FSTHDMIIGetProblems = LibraryFunctionLoad[libTHDMII, "FSTHDMIIGetProblems", LinkObject, LinkObject];
FSTHDMIIGetWarnings = LibraryFunctionLoad[libTHDMII, "FSTHDMIIGetWarnings", LinkObject, LinkObject];
FSTHDMIIToSLHA = LibraryFunctionLoad[libTHDMII, "FSTHDMIIToSLHA", LinkObject, LinkObject];

FSTHDMIIOpenHandleLib = LibraryFunctionLoad[libTHDMII, "FSTHDMIIOpenHandle", {{Real,1}}, Integer];
FSTHDMIICloseHandle = LibraryFunctionLoad[libTHDMII, "FSTHDMIICloseHandle", {Integer}, Void];

FSTHDMIISetLib = LibraryFunctionLoad[libTHDMII, "FSTHDMIISet", {Integer, {Real,1}}, Void];

FSTHDMIICalculateSpectrum = LibraryFunctionLoad[libTHDMII, "FSTHDMIICalculateSpectrum", LinkObject, LinkObject];
FSTHDMIICalculateObservables = LibraryFunctionLoad[libTHDMII, "FSTHDMIICalculateObservables", LinkObject, LinkObject];

FSTHDMIICalculateSpectrum::error = "`1`";
FSTHDMIICalculateSpectrum::warning = "`1`";

FSTHDMIICalculateObservables::error = "`1`";
FSTHDMIICalculateObservables::warning = "`1`";

FSTHDMII::info = "`1`";
FSTHDMII::nonum = "Error: `1` is not a numeric input value!";
FSTHDMIIMessage[s_] := Message[FSTHDMII::info, s];

FSTHDMIICheckIsNumeric[a_?NumericQ] := a;
FSTHDMIICheckIsNumeric[a_] := (Message[FSTHDMII::nonum, a]; Abort[]);

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

fsTHDMIIDefaultInputParameters = {
   Lambda1IN -> 0,
   Lambda2IN -> 0,
   Lambda3IN -> 0,
   Lambda4IN -> 0,
   Lambda5IN -> 0,
   Lambda6IN -> 0,
   Lambda7IN -> 0,
   M122IN -> 0,
   TanBeta -> 0,
   Qin -> 0
};

Options[FSTHDMIIOpenHandle] = {
    Sequence @@ fsDefaultSettings,
    Sequence @@ fsDefaultSMParameters,
    Sequence @@ fsTHDMIIDefaultInputParameters
};

FSTHDMIIOpenHandle[a___, (fsSettings | fsSMParameters | fsModelParameters) -> s_List, r___] :=
    FSTHDMIIOpenHandle[a, Sequence @@ s, r];

FSTHDMIIOpenHandle[OptionsPattern[]] :=
    FSTHDMIIOpenHandleLib[
        FSTHDMIICheckIsNumeric /@ {
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

            (* THDMII input parameters *)
            ,
            OptionValue[Lambda1IN],
            OptionValue[Lambda2IN],
            OptionValue[Lambda3IN],
            OptionValue[Lambda4IN],
            OptionValue[Lambda5IN],
            OptionValue[Lambda6IN],
            OptionValue[Lambda7IN],
            OptionValue[M122IN],
            OptionValue[TanBeta],
            OptionValue[Qin]
        }
];

Options[FSTHDMIISet] = Options[FSTHDMIIOpenHandle];

FSTHDMIISet[handle_Integer, a___, (fsSettings | fsSMParameters | fsModelParameters) -> s_List, r___] :=
    FSTHDMIISet[handle, a, Sequence @@ s, r];

FSTHDMIISet[handle_Integer, p:OptionsPattern[]] :=
    FSTHDMIISetLib[
        handle,
        ReleaseHold[Hold[FSTHDMIICheckIsNumeric /@ {
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

            (* THDMII input parameters *)
            ,
            OptionValue[Lambda1IN],
            OptionValue[Lambda2IN],
            OptionValue[Lambda3IN],
            OptionValue[Lambda4IN],
            OptionValue[Lambda5IN],
            OptionValue[Lambda6IN],
            OptionValue[Lambda7IN],
            OptionValue[M122IN],
            OptionValue[TanBeta],
            OptionValue[Qin]
        }] /. HoldPattern[OptionValue[param_]] :> param /.
        { p } /.
        FSTHDMIIGetSettings[handle] /.
        FSTHDMIIGetSMInputParameters[handle] /.
        FSTHDMIIGetInputParameters[handle]]];
