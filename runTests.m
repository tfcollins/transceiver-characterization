
%% Run individual tests
%results = run(TransceiverCharacterizationTests,'testEVMR7_SDRtoPXA');
%results = run(TransceiverCharacterizationTests,'testACLRR7_SDRtoPXA');
%results = run(TransceiverCharacterizationTests,'testEVMR7_SigGenToPXA');
results = run(TransceiverCharacterizationTests,'testEVMR7_SigGenToSDR');

%% Run full harness
% import matlab.unittest.TestRunner;
% import matlab.unittest.TestSuite;
% import matlab.unittest.plugins.TestReportPlugin;
% 
% %suite = testsuite('TransceiverCharacterizationTests');
% suite = TestSuite.fromClass(?TransceiverCharacterizationTests);
% 
% runner = TestRunner.withTextOutput;
% pdfFile = 'MyTestReport2.pdf';
% plugin = TestReportPlugin.producingPDF(pdfFile,...
%     'IncludingPassingDiagnostics',true,'IncludingCommandWindowText',true);
% 
% %results1 = run(DocPolynomTest,'testMultiplication');
% 
% runner.addPlugin(plugin);
% result = runner.run(suite,'testEVMR7_SDRtoPXAObjs');


