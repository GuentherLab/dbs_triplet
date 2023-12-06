
PROTOCOL_PATH = 'Z:\Users\busha\Analysis\2020-05-04-TimeFreq-DBS3000';
PROTOCOL_FUNCTION = 'P10_time_frequency_arC_lock_speech_20201019';
PROTOCOL_TABLE = 'P10_Subjets.txt';
exe_daytime = datestr(now,'yyyymmdd_HHMM');
diary([PROTOCOL_PATH filesep 'batch_' PROTOCOL_FUNCTION '_' exe_daytime '.log'])
addpath(PROTOCOL_PATH);

PATH_DATA = 'Z:\DBS';
cd(PROTOCOL_PATH)

subject_table = readtable(PROTOCOL_TABLE);         
fprintf('=== Running protocol %s ===\n',PROTOCOL_FUNCTION)

for i=1:height(subject_table)
  SUBJECT = subject_table.subject{i};
  fprintf('Running protocol.')
  %running protocol
  try
    proto = str2func(PROTOCOL_FUNCTION);
    proto(SUBJECT);
    fprintf('OK\n')
  catch err
    fprintf('FAILED: %s\n',err.message)
  end
end

diary('off')
