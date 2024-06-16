% Base directory
baseDir = 'Z:\DBS';

% Get a list of all subdirectories in the base directory
d = dir(baseDir);
isub = [d(:).isdir]; % Returns logical vector
subFolders = {d(isub).name}';
subFolders(ismember(subFolders,{'.','..'})) = []; % Remove . and ..

% Initialize skipped folders
skippedFolders = {};
results = {};

% Loop through each subfolder
for i = 1:length(subFolders)
    if startsWith(subFolders{i}, 'DBS') && ~isempty(regexp(subFolders{i}, '^DBS\d+', 'once'))
        currentFolder = fullfile(baseDir, subFolders{i});
        preprocessedDir = fullfile(currentFolder, 'preprocessed data');
       
        % Check if the Preprocessed Data folder exists
        if exist(preprocessedDir, 'dir')
            annotFile = fullfile(preprocessedDir, 'sync', 'annot', [subFolders{i}, '_produced_phoneme.txt']);
            % Check if the specific file exists
            if exist(annotFile, 'file')
                fid = fopen(annotFile, 'r');
                onTargetCount = 0;
                if fid ~= -1
                    header = fgets(fid); % Read the header line
                    headerCols = strsplit(strtrim(header), '\t'); % Split header into columns
                    phoneticOnTargetCol = find(strcmp(headerCols, 'phonetic_ontarget')); % Find index of the 'phonetic_ontarget' column
                   
                    % Check if 'phonetic_ontarget' column exists
                    if isempty(phoneticOnTargetCol)
                        skippedFolders{end+1} = subFolders{i};
                        fclose(fid);
                        continue;
                    end
                   
                    % Read the file line by line after the header
                    while ~feof(fid)
                        line = fgets(fid);
                        data = strsplit(strtrim(line), '\t'); % Split data into columns
                        if length(data) >= phoneticOnTargetCol % Check if 'phonetic_ontarget' column exists in this line
                            % Increment counts based on the value in the 'phonetic_ontarget' column
                            if strcmp(data{phoneticOnTargetCol}, '1')
                                onTargetCount = onTargetCount + 1;
                            end
                        end
                    end
                    fclose(fid);
                end
                results{end+1} = sprintf('%s, %d', subFolders{i}, onTargetCount);
            else
                skippedFolders{end+1} = subFolders{i};
            end
        else
            skippedFolders{end+1} = subFolders{i};
        end
    end
end

% Output skipped folders
if ~isempty(skippedFolders)
    fprintf('Folders skipped: %s\n', strjoin(skippedFolders, ', '));
end

% Write results to a file
resultsFile = fullfile(baseDir, 'phonetic_ontarget_results.txt');
fid = fopen(resultsFile, 'wt');
if fid ~= -1
    fprintf(fid, 'DBS#, On-Target Count\n');
    for i = 1:length(results)
        fprintf(fid, '%s\n', results{i});
    end
    fclose(fid);
else
    error('Failed to open results file.');
end
