%% Loading paths
ft_defaults
bml_defaults
format long
clear

%% Configuration Variables and Paths
PATH_ANALYSIS = '/Volumes/Nexus/Commits/Vibration_artifacts/audio_p-coherence_syllable_triplet';
PATH_DATA='/Volumes/Nexus/DBS';
PATH_AVERAGE_MNI = [PATH_DATA '/DBS_subject_lists/MNI_ICBM_2009b_NLIN_ASYM/cortex/CortexLowRes_15000V.mat'];
PATH_SUBCORT_ATLAS = '/Volumes/Nexus/Resources/STN-Atlas/atlas_index.mat';

cd(PATH_ANALYSIS)
electrode = readtable('data/electrode_session_coherence.tsv','Delimiter', '\t', 'TreatAsEmpty', 'NA','FileType','text');

%loading cortical reconstructions
average_mni = load(PATH_AVERAGE_MNI);

subcort = load(PATH_SUBCORT_ATLAS);


%% 

%% ecog

sel_electrode = electrode(strcmp(electrode.type,'ecog') & (electrode.mni_nonlinear_x < 0) & contains(electrode.subject,'DBS3') ,:);
cmap = colormap('hot');

figure;
hold on
patch('vertices', average_mni.Vertices, 'faces', average_mni.Faces,...
'FaceColor', [.5 .5 .6], 'EdgeColor', 'none', 'FaceAlpha',1, ...
'facelighting', 'gouraud', 'specularstrength', 0, 'ambientstrength', 0.5, 'diffusestrength', 0.5)
hold on
for i = 1:25
    electrode_cmap = sel_electrode((sel_electrode.rNC > i-1) & (sel_electrode.rNC <= i),:);
    scatter3(electrode_cmap.mni_nonlinear_x-15-i, electrode_cmap.mni_nonlinear_y, electrode_cmap.mni_nonlinear_z,'filled',...
      'MarkerFaceAlpha',0.2+0.8*i/25,'MarkerFaceColor',cmap(i*10,:))
end
electrode_cmap = sel_electrode((sel_electrode.rNC > 25),:);
scatter3(electrode_cmap.mni_nonlinear_x-15-26, electrode_cmap.mni_nonlinear_y, electrode_cmap.mni_nonlinear_z,'filled',...
     'MarkerFaceAlpha',1,'MarkerFaceColor',cmap(255,:))
view(-90,0)
axis off; axis equal
camlight('headlight','infinite');
cmap = colormap('hot');
colorbar;
scalebar(0,70,-50, 10, 'mm')

saveas(gcf,[PATH_ANALYSIS '/fig/A16_01_mni_nonlin_electrodes_rNC_overlay-with-colorbar.png'])


%% micro

sel_electrode = electrode(strcmp(electrode.type,'micro') & (electrode.mni_nonlinear_x < 0) & contains(electrode.subject,'DBS3') ,:);

side=2; %left
figure;
hold on
patch('vertices', subcort.atlases.fv{1,side}.vertices, 'faces', subcort.atlases.fv{1,side}.faces,...
'FaceColor', [.5 .5 .6], 'EdgeColor', 'none', 'FaceAlpha',0.5, ...
'facelighting', 'gouraud', 'specularstrength', 0, 'ambientstrength', 0.5, 'diffusestrength', 0.5)
hold on

hold on
for i = 1:25
    electrode_cmap = sel_electrode((sel_electrode.rNC > i-1) & (sel_electrode.rNC <= i),:);
    scatter3(electrode_cmap.mni_nonlinear_x, electrode_cmap.mni_nonlinear_y, electrode_cmap.mni_nonlinear_z,'filled','o',...
      'MarkerFaceColor',cmap(i*10,:))
end
electrode_cmap = sel_electrode((sel_electrode.rNC > 25),:);
scatter3(electrode_cmap.mni_nonlinear_x, electrode_cmap.mni_nonlinear_y, electrode_cmap.mni_nonlinear_z,'filled','o',...
     'MarkerFaceColor',cmap(255,:))

view(-90,0)
axis off; axis equal
%light('Position',[-1.101435411230083e+02 11.998238489140135 51.697756317087205],'Style','local')
camlight('headlight');
% lightangle(0,-90);
cmap = colormap('hot');
colorbar;
scalebar(-10,-7,-12, 2, 'mm')

saveas(gcf,[PATH_ANALYSIS '/fig/A16_02_mni_nonlin_electrodes_STN_micro.png'])
close(gcf)





%% macro

sel_electrode = electrode(strcmp(electrode.type,'macro') & (electrode.mni_nonlinear_x < 0) & contains(electrode.subject,'DBS3') ,:);

side=2; %left
figure;
hold on
patch('vertices', subcort.atlases.fv{1,side}.vertices, 'faces', subcort.atlases.fv{1,side}.faces,...
'FaceColor', [.5 .5 .6], 'EdgeColor', 'none', 'FaceAlpha',0.5, ...
'facelighting', 'gouraud', 'specularstrength', 0, 'ambientstrength', 0.5, 'diffusestrength', 0.5)
hold on

hold on
for i = 1:25
    electrode_cmap = sel_electrode((sel_electrode.rNC > i-1) & (sel_electrode.rNC <= i),:);
    scatter3(electrode_cmap.mni_nonlinear_x, electrode_cmap.mni_nonlinear_y, electrode_cmap.mni_nonlinear_z,'filled','o',...
      'MarkerFaceColor',cmap(i*10,:))
end
electrode_cmap = sel_electrode((sel_electrode.rNC > 25),:);
scatter3(electrode_cmap.mni_nonlinear_x, electrode_cmap.mni_nonlinear_y, electrode_cmap.mni_nonlinear_z,'filled','o',...
     'MarkerFaceColor',cmap(255,:))

view(-90,0)
axis off; axis equal
%light('Position',[-1.101435411230083e+02 11.998238489140135 51.697756317087205],'Style','local')
camlight('headlight');
% lightangle(0,-90);
cmap = colormap('hot');
colorbar;
scalebar(-10,-5,-13, 2, 'mm')

saveas(gcf,[PATH_ANALYSIS '/fig/A16_03_mni_nonlin_electrodes_STN_macro.png'])
close(gcf)



%% dbs

sel_electrode = electrode(strcmp(electrode.type,'dbs') & (electrode.mni_nonlinear_x < 0) & contains(electrode.subject,'DBS3') ,:);

side=2; %left
figure;
hold on
patch('vertices', subcort.atlases.fv{1,side}.vertices, 'faces', subcort.atlases.fv{1,side}.faces,...
'FaceColor', [.5 .5 .6], 'EdgeColor', 'none', 'FaceAlpha',0.5, ...
'facelighting', 'gouraud', 'specularstrength', 0, 'ambientstrength', 0.5, 'diffusestrength', 0.5)
hold on

hold on
for i = 1:25
    electrode_cmap = sel_electrode((sel_electrode.rNC > i-1) & (sel_electrode.rNC <= i),:);
    scatter3(electrode_cmap.mni_nonlinear_x, electrode_cmap.mni_nonlinear_y, electrode_cmap.mni_nonlinear_z,'filled','o',...
      'MarkerFaceColor',cmap(i*10,:))
end
electrode_cmap = sel_electrode((sel_electrode.rNC > 25),:);
scatter3(electrode_cmap.mni_nonlinear_x, electrode_cmap.mni_nonlinear_y, electrode_cmap.mni_nonlinear_z,'filled','o',...
     'MarkerFaceColor',cmap(255,:))

view(-90,0)
axis off; axis equal
%light('Position',[-1.101435411230083e+02 11.998238489140135 51.697756317087205],'Style','local')
camlight('headlight');
% lightangle(0,-90);
cmap = colormap('hot');
colorbar;
scalebar(-10,-5,-13, 2, 'mm')

saveas(gcf,[PATH_ANALYSIS '/fig/A16_04_mni_nonlin_electrodes_STN_dbs.png'])
close(gcf)







sel_electrode = electrode(strcmp(electrode.type,'dbs') & (electrode.mni_nonlinear_x < 0) & contains(electrode.subject,'DBS3') & (electrode.mni_nonlinear_z >3) ,:);




