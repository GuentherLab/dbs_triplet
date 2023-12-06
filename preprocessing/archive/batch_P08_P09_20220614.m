
PROTOCOL_PATH = 'C:\Users\amsmeier\Documents\MATLAB\P09_artifact_criteria_E';
PROTOCOL_FUNCTION = 'P09_detect_artifact_criteria_E';
PROTOCOL_TABLE = 'P09_Subjects.txt';
exe_daytime = datestr(now,'yyyymmdd_HHMM');
addpath(PROTOCOL_PATH);
diary([PROTOCOL_PATH filesep 'batch_' PROTOCOL_FUNCTION '_' exe_daytime '.log'])

SKIP_OK = false; %Should previous protocols run by this script successfully be skipped
FORCE = true; %Archive all previous versions of the script and run current 
              %overrides any manual modification
%SUBJECTS = {'DBS3001'};

PATH_DATA = 'Z:\DBS';
cd(PROTOCOL_PATH)
% checksum = Simulink.getFileChecksum([PROTOCOL_PATH filesep PROTOCOL_FUNCTION '.m']);

subject_table = readtable(PROTOCOL_TABLE);         
fprintf('=== Running protocol %s ===\n',PROTOCOL_FUNCTION)
if FORCE; fprintf('Forced run, overwritting any manual change.\n'); end
if SKIP_OK; fprintf('Skipping previously successfully executed protocols.\n'); end

sub_inds_to_run = 1:height(subject_table);
% sub_inds_to_run = 10;

for isub = sub_inds_to_run
  SUBJECT = subject_table.subject{isub};
  timetol = subject_table.timetol(isub);
  
  PATH_SUBJECT=[PATH_DATA filesep SUBJECT];
  PATH_PREPROCESSED=[PATH_SUBJECT filesep 'Preprocessed Data'];
  PATH_SYNC=[PATH_PREPROCESSED filesep 'Sync'];
  PATH_ARCHIVE=[PATH_SYNC filesep 'archive'];
  fprintf('\n* Subject %s: ',SUBJECT)
  
  %loading simple file verification
  
  
        % AM changed path because permissions for sync and annot folders are not available
% % % %   sfv_fname = [PATH_SYNC filesep 'protocols.sfv'];
    sfv_fname = [PROTOCOL_PATH filesep 'protocols.sfv'];

  
  
  
  if isfile(sfv_fname)
    opts = detectImportOptions(sfv_fname,'FileType','text');
    opts.Delimiter = {'\t'};
    opts.VariableNamesLine = 0;
    opts.DataLine = 1;
    opts.VariableNames = {'filename','checksum'};
    sfv = readtable(sfv_fname,opts);    
    fprintf('Loaded sfv file. ')
  else
    sfv = cell2table(cell(0,2),'VariableNames',{'filename','checksum'});
    fprintf('No sfv file, ')
  end
  
	new_proto_fname = [SUBJECT '_' PROTOCOL_FUNCTION '_' exe_daytime]; 
  old_proto = dir([PATH_SYNC filesep SUBJECT '_' PROTOCOL_FUNCTION '*.m']);
  
  merge = false; 
  skip = false;
  if ~isempty(old_proto)
    fprintf('Previous script(s) found. ')
    %archiving previous runs of this script
    for j=1:length(old_proto)
      if FORCE
        fprintf('Archiving. ')
        movefile([old_proto(j).folder filesep old_proto(j).name],...
                 [PATH_SYNC filesep 'archive' filesep old_proto(j).name])         
      elseif ismember(old_proto(j).name,sfv.filename)
        %file created by this script
        old_checksum_j = Simulink.getFileChecksum([old_proto(j).folder filesep old_proto(j).name]);       
        if strcmp(sfv.checksum(strcmp(old_proto(j).name,sfv.filename)),old_checksum_j)
          if ~endsWith(old_proto(j).name,["_TO_MERGE.m","_FAILED.m"]) && SKIP_OK
            skip = true;
          else
            fprintf('Archiving. ')
            movefile([old_proto(j).folder filesep old_proto(j).name],...
                     [PATH_SYNC filesep 'archive' filesep old_proto(j).name])           
          end
        else
          merge=true;
        end
      else  
        merge=true;
      end
    end
  end
  
  if skip
    fprintf('SKIPPING\n')
  elseif merge
    fprintf('\nProtocol in subject''s folder modified manually. Copying protocol to folder for manual merging.\n')
    
    copyfile([PROTOCOL_PATH filesep PROTOCOL_FUNCTION '.m'],...
             [PATH_SYNC filesep new_proto_fname '_TO_MERGE.m']);

    %updating sfv
    sfv = [sfv; 
           cell2table({[new_proto_fname '.m'],checksum},'VariableNames',{'filename','checksum'}); 
           cell2table({[new_proto_fname '_TO_MERGE.m'],checksum},'VariableNames',{'filename','checksum'})];
  else
    fprintf('Running protocol.')
       
    %running protocol
    try
        get_highgamma_from_denoised(SUBJECT,timetol); 
        
% % % % % % % % % % % % % %       proto = str2func(PROTOCOL_FUNCTION);
% % % % % % % % % % % % % %       proto(SUBJECT,timetol);
% % % % % % % % % % % % % %       fprintf('OK\n')
      
      
      % AM commented out because permissions for sync and annot folders are not available
      
% % % % % % %       %copying current version of the script to the subjects folder
% % % % % % %       copyfile([PROTOCOL_PATH filesep PROTOCOL_FUNCTION '.m'],...
% % % % % % %                [PATH_SYNC filesep new_proto_fname '.m']);  
           
           
            % AM commented out because simulink checksum generating errors
           
% % % % % % %       %updating sfv
% % % % % % %       sfv = [sfv; cell2table({[new_proto_fname '.m'],checksum},'VariableNames',{'filename','checksum'})];    
      
    catch err
      fprintf('FAILED: %s\n',err.message)
      copyfile([PROTOCOL_PATH filesep PROTOCOL_FUNCTION '.m'],...
               [PATH_SYNC filesep new_proto_fname '_FAILED.m']);  
      %updating sfv
      sfv = [sfv; 
             cell2table({[new_proto_fname '.m'],checksum},'VariableNames',{'filename','checksum'}); 
             cell2table({[new_proto_fname '_FAILED.m'],checksum},'VariableNames',{'filename','checksum'})];
    end

  end
  
	writetable(sfv,sfv_fname,'FileType','text','Delimiter','\t','WriteVariableNames',false);
end
diary('off')
