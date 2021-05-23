Print["================================"];
Print["FlexibleSUSY 2.6.0"];
Print["MSSMNoFVatMGUTHimalaya"];
Print["http://flexiblesusy.hepforge.org"];
Print["================================"];

libMSSMNoFVatMGUTHimalaya = FileNameJoin[{Directory[], "models", "MSSMNoFVatMGUTHimalaya", "MSSMNoFVatMGUTHimalaya_librarylink.so"}];

FSMSSMNoFVatMGUTHimalayaGetSettings = LibraryFunctionLoad[libMSSMNoFVatMGUTHimalaya, "FSMSSMNoFVatMGUTHimalayaGetSettings", LinkObject, LinkObject];
FSMSSMNoFVatMGUTHimalayaGetSMInputParameters = LibraryFunctionLoad[libMSSMNoFVatMGUTHimalaya, "FSMSSMNoFVatMGUTHimalayaGetSMInputParameters", LinkObject, LinkObject];
FSMSSMNoFVatMGUTHimalayaGetInputParameters = LibraryFunctionLoad[libMSSMNoFVatMGUTHimalaya, "FSMSSMNoFVatMGUTHimalayaGetInputParameters", LinkObject, LinkObject];
FSMSSMNoFVatMGUTHimalayaGetProblems = LibraryFunctionLoad[libMSSMNoFVatMGUTHimalaya, "FSMSSMNoFVatMGUTHimalayaGetProblems", LinkObject, LinkObject];
FSMSSMNoFVatMGUTHimalayaGetWarnings = LibraryFunctionLoad[libMSSMNoFVatMGUTHimalaya, "FSMSSMNoFVatMGUTHimalayaGetWarnings", LinkObject, LinkObject];
FSMSSMNoFVatMGUTHimalayaToSLHA = LibraryFunctionLoad[libMSSMNoFVatMGUTHimalaya, "FSMSSMNoFVatMGUTHimalayaToSLHA", LinkObject, LinkObject];

FSMSSMNoFVatMGUTHimalayaOpenHandleLib = LibraryFunctionLoad[libMSSMNoFVatMGUTHimalaya, "FSMSSMNoFVatMGUTHimalayaOpenHandle", {{Real,1}}, Integer];
FSMSSMNoFVatMGUTHimalayaCloseHandle = LibraryFunctionLoad[libMSSMNoFVatMGUTHimalaya, "FSMSSMNoFVatMGUTHimalayaCloseHandle", {Integer}, Void];

FSMSSMNoFVatMGUTHimalayaSetLib = LibraryFunctionLoad[libMSSMNoFVatMGUTHimalaya, "FSMSSMNoFVatMGUTHimalayaSet", {Integer, {Real,1}}, Void];

FSMSSMNoFVatMGUTHimalayaCalculateSpectrum = LibraryFunctionLoad[libMSSMNoFVatMGUTHimalaya, "FSMSSMNoFVatMGUTHimalayaCalculateSpectrum", LinkObject, LinkObject];
FSMSSMNoFVatMGUTHimalayaCalculateObservables = LibraryFunctionLoad[libMSSMNoFVatMGUTHimalaya, "FSMSSMNoFVatMGUTHimalayaCalculateObservables", LinkObject, LinkObject];

FSMSSMNoFVatMGUTHimalayaCalculateSpectrum::error = "`1`";
FSMSSMNoFVatMGUTHimalayaCalculateSpectrum::warning = "`1`";

FSMSSMNoFVatMGUTHimalayaCalculateObservables::error = "`1`";
FSMSSMNoFVatMGUTHimalayaCalculateObservables::warning = "`1`";

FSMSSMNoFVatMGUTHimalaya::info = "`1`";
FSMSSMNoFVatMGUTHimalaya::nonum = "Error: `1` is not a numeric input value!";
FSMSSMNoFVatMGUTHimalayaMessage[s_] := Message[FSMSSMNoFVatMGUTHimalaya::info, s];

FSMSSMNoFVatMGUTHimalayaCheckIsNumeric[a_?NumericQ] := a;
FSMSSMNoFVatMGUTHimalayaCheckIsNumeric[a_] := (Message[FSMSSMNoFVatMGUTHimalaya::nonum, a]; Abort[]);

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

fsMSSMNoFVatMGUTHimalayaDefaultInputParameters = {
   TanBeta -> 0,
   SignMu -> 0,
   M1 -> 0,
   M2 -> 0,
   M3 -> 0,
   AtIN -> 0,
   AbIN -> 0,
   AtauIN -> 0,
   AcIN -> 0,
   AsIN -> 0,
   AmuonIN -> 0,
   AuIN -> 0,
   AdIN -> 0,
   AeIN -> 0,
   mHd2IN -> 0,
   mHu2IN -> 0,
   ml11IN -> 0,
   ml22IN -> 0,
   ml33IN -> 0,
   me11IN -> 0,
   me22IN -> 0,
   me33IN -> 0,
   mq11IN -> 0,
   mq22IN -> 0,
   mq33IN -> 0,
   mu11IN -> 0,
   mu22IN -> 0,
   mu33IN -> 0,
   md11IN -> 0,
   md22IN -> 0,
   md33IN -> 0,
   Mlow -> 0
};

Options[FSMSSMNoFVatMGUTHimalayaOpenHandle] = {
    Sequence @@ fsDefaultSettings,
    Sequence @@ fsDefaultSMParameters,
    Sequence @@ fsMSSMNoFVatMGUTHimalayaDefaultInputParameters

};

FSMSSMNoFVatMGUTHimalayaOpenHandle[a___, (fsSettings | fsSMParameters | fsModelParameters ) -> s_List, r___] :=
    FSMSSMNoFVatMGUTHimalayaOpenHandle[a, Sequence @@ s, r];

FSMSSMNoFVatMGUTHimalayaOpenHandle[OptionsPattern[]] :=
    FSMSSMNoFVatMGUTHimalayaOpenHandleLib[
        FSMSSMNoFVatMGUTHimalayaCheckIsNumeric /@ {
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

            (* MSSMNoFVatMGUTHimalaya input parameters *)
            ,
            OptionValue[TanBeta],
            OptionValue[SignMu],
            OptionValue[M1],
            OptionValue[M2],
            OptionValue[M3],
            OptionValue[AtIN],
            OptionValue[AbIN],
            OptionValue[AtauIN],
            OptionValue[AcIN],
            OptionValue[AsIN],
            OptionValue[AmuonIN],
            OptionValue[AuIN],
            OptionValue[AdIN],
            OptionValue[AeIN],
            OptionValue[mHd2IN],
            OptionValue[mHu2IN],
            OptionValue[ml11IN],
            OptionValue[ml22IN],
            OptionValue[ml33IN],
            OptionValue[me11IN],
            OptionValue[me22IN],
            OptionValue[me33IN],
            OptionValue[mq11IN],
            OptionValue[mq22IN],
            OptionValue[mq33IN],
            OptionValue[mu11IN],
            OptionValue[mu22IN],
            OptionValue[mu33IN],
            OptionValue[md11IN],
            OptionValue[md22IN],
            OptionValue[md33IN],
            OptionValue[Mlow]

        }
];

Options[FSMSSMNoFVatMGUTHimalayaSet] = Options[FSMSSMNoFVatMGUTHimalayaOpenHandle];

FSMSSMNoFVatMGUTHimalayaSet[handle_Integer, a___, (fsSettings | fsSMParameters | fsModelParameters) -> s_List, r___] :=
    FSMSSMNoFVatMGUTHimalayaSet[handle, a, Sequence @@ s, r];

FSMSSMNoFVatMGUTHimalayaSet[handle_Integer, p:OptionsPattern[]] :=
    FSMSSMNoFVatMGUTHimalayaSetLib[
        handle,
        ReleaseHold[Hold[FSMSSMNoFVatMGUTHimalayaCheckIsNumeric /@ {
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

            (* MSSMNoFVatMGUTHimalaya input parameters *)
            ,
            OptionValue[TanBeta],
            OptionValue[SignMu],
            OptionValue[M1],
            OptionValue[M2],
            OptionValue[M3],
            OptionValue[AtIN],
            OptionValue[AbIN],
            OptionValue[AtauIN],
            OptionValue[AcIN],
            OptionValue[AsIN],
            OptionValue[AmuonIN],
            OptionValue[AuIN],
            OptionValue[AdIN],
            OptionValue[AeIN],
            OptionValue[mHd2IN],
            OptionValue[mHu2IN],
            OptionValue[ml11IN],
            OptionValue[ml22IN],
            OptionValue[ml33IN],
            OptionValue[me11IN],
            OptionValue[me22IN],
            OptionValue[me33IN],
            OptionValue[mq11IN],
            OptionValue[mq22IN],
            OptionValue[mq33IN],
            OptionValue[mu11IN],
            OptionValue[mu22IN],
            OptionValue[mu33IN],
            OptionValue[md11IN],
            OptionValue[md22IN],
            OptionValue[md33IN],
            OptionValue[Mlow]
        }] /. HoldPattern[OptionValue[param_]] :> param /.
        { p } /.
        FSMSSMNoFVatMGUTHimalayaGetSettings[handle] /.
        FSMSSMNoFVatMGUTHimalayaGetSMInputParameters[handle] /.
        FSMSSMNoFVatMGUTHimalayaGetInputParameters[handle]]];
