dataDir = 'Path\to\your\data';
pgObj = ProtocolGetter(dataDir);
pgObj.getConditionSignals;
pgObj.getSignalEdges;
pgObj.getFrequencySubs;
pgObj.pairStimulus;
% In between, you should verify the Conditions property of the pgObj
% variable if it matches the real stimulation protocol. If it mathces, then
% continue.
pgObj.saveConditions;