% copy data from the triplet server to new drive location

% to do:
%%% make option for use in copying files back from external drive to the server; if files already exist on the server, only replace them if I am the owner of the versions on the server and they were created after ~2021
%%% also copy group analysis folder
%%% add option for copying only small files (sort by filesize or filetype), to copy updates that matteo and others are adding to phon coding, speech error coding, artifact marking, etc


clear

setpaths_dbs_triplet()


% topdir_to_copy_to = '\\tsclient\D\triplet'; 
% topdir_to_copy_to = '\\tsclient\E\triplet'; 
topdir_to_copy_to = 'D:\triplet';

%%%%% this version just copies all 3000-series subs rather than using a predefined subject table
dbpath = 'Z:\DBS'; 
 % % % % % % % dd = struct2table(dir(dbpath));
 % % % % % % % subject = dd.name(contains(dd.name,'DBS3')); %%% subjects from the 3000 (triplet) series
 % % % % % % % subs = table(subject); clear subject


subs = readtable([PATH_ARTIFACT filesep 'P08_Subjects_to_analyze.txt']);  % created by generate_triplet_subject_list.m

 nsubs = height(subs);

 for isub = 29:nsubs
     subtic = tic; 
     thissub = subs.subject{isub}
     dir_to_copy = [dbpath filesep thissub filesep 'Preprocessed Data']; 
     new_sub_dir = [topdir_to_copy_to filesep thissub filesep 'Preprocessed Data'];
     mkdir(new_sub_dir); % create subject directory in the new location
     copyfile(dir_to_copy,new_sub_dir)
     subtoc = toc(subtic)
 end
     
     