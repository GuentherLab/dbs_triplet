%This mainscript is what will % Base directory
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
                accurateCount = 0;
                inaccurateCount = 0;
                if fid ~= -1
                    header = fgets(fid); % Read the header line
                    headerCols = strsplit(strtrim(header), '\t'); % Split header into columns
                    accuracyCol = find(strcmp(headerCols, 'accuracy')); % Find index of the 'accuracy' column
                   
                    % Read the file line by line after the header
                    while ~feof(fid)
                        line = fgets(fid);
                        data = strsplit(strtrim(line), '\t'); % Split data into columns
                        if length(data) >= accuracyCol % Check if 'accuracy' column exists in this line
                            % Increment counts based on the value in the 'accuracy' column
                            if strcmp(data{accuracyCol}, 'accurate')
                                accurateCount = accurateCount + 1;
                            elseif strcmp(data{accuracyCol}, 'inaccurate')
                                inaccurateCount = inaccurateCount + 1;
                            end
                        end
                    end
                    fclose(fid);
                end
                results{end+1} = sprintf('%s, %d, %d', subFolders{i}, inaccurateCount, accurateCount);
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
resultsFile = fullfile(baseDir, 'accuracy_results.txt');
fid = fopen(resultsFile, 'wt');
if fid ~= -1
    fprintf(fid, 'DBS#, Inaccurate, Accurate\n');
    for i = 1:length(results)
        fprintf(fid, '%s\n', results{i});
    end
    fclose(fid);
else
    error('Failed to open results file.');
end
