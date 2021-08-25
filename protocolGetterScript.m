dataDir = 'Z:\Jesus\Jittering\201016_jiterring_TRN_ap1500_ml1500_dv1370_and_bad_S1\TRN';
pgObj = ProtocolGetter(dataDir);
pgObj.getConditionSignals;
pgObj.getSignalEdges;
pgObj.getFrequencyEdges;
pgObj.pairStimulus;
% In between, you should verify the Conditions property of the pgObj
% variable if it matches the real stimulation protocol. If it mathces, then
% continue.
pgObj.saveConditions;