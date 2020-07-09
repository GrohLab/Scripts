%calling all scripts
msgbox('select folder')
NewFigureDir=uigetdir;
figureDirJesus = fullfile(NewFigureDir,'New figures Jesus\');
DataDirL200=mkdir(figureDirJesus,'\ControlandL200');
DataDirL100=mkdir(figureDirJesus,'\ControlandL100');
DataDirL50=mkdir(figureDirJesus,'\ControlandL50');
DataDirL10=mkdir(figureDirJesus,'\ControlandL10');
dataDirL1=mkdir(figureDirJesus,'\ControlandL1');
msgbox('select ""relativeSpkTmsStruct" file to load')
JesusFiles=uigetfile;

load(JesusFiles);

run("C:\Users\jesus\Documents\MATLAB\Scripts\Jesus\clusterstomat.m")
%%
run("C:\Users\jesus\Documents\MATLAB\Scripts\Jesus\PSTHandBoxPlotsCONTROLandLASERconditions.m")
%%
run("C:\Users\jesus\Documents\MATLAB\Scripts\Jesus\BeforeandAfterJesusPlots.m")
%%
run("C:\Users\jesus\Documents\MATLAB\Scripts\Jesus\cdfForAllConditions.m")


