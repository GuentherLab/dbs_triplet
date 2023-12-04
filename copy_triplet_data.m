% copy data from the triplet server to new drive location

dbpath = 'Z:\DBS'; 

% topdir_to_copy_to = '\\tsclient\D\triplet'; 
topdir_to_copy_to = '\\tsclient\E\triplet'; 
% topdir_to_copy_to = 'D:\triplet';

 dd = struct2table(dir(dbpath));
 subs = dd.name(contains(dd.name,'DBS3')); %%% subjects from the 3000 (triplet) series
 nsubs = length(subs);

 for isub = 7:nsubs
     thissub = subs{isub}
     dir_to_copy = [dbpath filesep thissub filesep 'Preprocessed Data']; 
     new_sub_dir = [topdir_to_copy_to filesep thissub];
     mkdir(new_sub_dir); % create subject directory in the new location
     copyfile(dir_to_copy,new_sub_dir)
 end
     
     