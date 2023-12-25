% a prior version of copy_triplet_data.m forget to put data in a 'Preprocessed Data' folder;....
% ....this script will create 'Preprocessed Data' folder and move data there

% do not run this script on TURBO - it involves chaning folder names

setpaths_dbs_triplet()

% topdir_to_copy_to = '\\tsclient\D\triplet'; 
% topdir_to_copy_to = '\\tsclient\E\triplet'; 
topdir_to_copy_to = PATH_DATA;

 dd = struct2table(dir(PATH_DATA));
 subs = dd.name(contains(dd.name,'DBS3')); %%% subjects from the 3000 (triplet) series
 nsubs = length(subs);

 for isub = 2:nsubs
     subtic = tic; 

     thissub = subs{isub}
     dir_to_move = [PATH_DATA filesep thissub]; 
     new_sub_dir = [topdir_to_copy_to filesep thissub filesep 'Preprocessed Data'];
     movefile(dir_to_move,new_sub_dir)

     subtoc = toc(subtic)
 end
     
     