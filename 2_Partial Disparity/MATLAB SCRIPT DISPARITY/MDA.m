function mda
%*************************************************************************
%*   MDA: Morphospace-Disparity Analysis                                 *
%*   MDA Version 1.2(C)2001 Nicolas Navarro                              *
%* ***********************************************************************
global filename         % Default name of input file (scores on principal components)
global Eigenfile        % Default name of input file (eigenvalues for standardization of PCs)
global Occurrencefile                % Default name of input file (temporal or geographical...occurrences)
global Groupfile                    % Default name of input file (group member)

if not(exist('IntrMDA.mat'))
   [IntrImage, IntrMap] = imread( 'imagMDAvs1.2.bmp' );
   save IntrMDA IntrImage IntrMap;
   delete('imagMDAvs1.2.bmp');
end;
load( 'IntrMDA.mat' );
[NbL, NbC, NbP] = size( IntrImage );

IntroFigure = figure( ...
		'Visible', 'off', ...
		'NumberTitle', 'off', ...
		'Units', 'pixels', ...
	   'Position',[300 500 475 200], ...
		'Color', [0.8 0.8 0.8], ...
		'Colormap', IntrMap, ...
		'KeyPressFcn', 'close(''gcf'')', ...
		'WindowStyle', 'modal' );
set( gca,'Position', [0 0 1 1] ); 
IntroImage = image( IntrImage );  axis off; axis image;
set( IntroImage, 'ButtonDownFcn', 'close(''gcf'')' );
set( IntroFigure, 'Visible', 'on' );

% wait user action
uiwait(IntroFigure);
clear all;
init=1;

while init==1
stp=not(exist('stop.mat'));
if stp~=1
   delete('stop.mat');
   break;
end;

% MDA version and copyright
MDA_version;
pause;

% loading default values of various parameters of analysis
MDA_prefs;

%Text user interface for Morphospace-Disparity Analysis
if not(exist(filename))
   [filename, pathname]=uigetfile('*.txt','Load TXT file: scores on principal components');
   if filename~=0
      filename=[pathname filename];
      PCs=load (filename);
      if isempty(PCs)
         string=['File' filename ' is empty!!!.'];
         title='File Control Failed';
         errordlg(string,title,'modal');
      end; 
   end;
else
   PCs= load (filename);
end;
[lig coll]=size(PCs);

%Option1: Rescaling PCs
disp('   PCs need rescaling PCs to eigenvalues: y/n');
disp('   Help:');
disp('   If you use data from some statistical package, all variances on axes can have be  ');
disp('   scale to unity. In this case the rescaling option (using the root squared of the  ');
disp('   eigenvalue)is used to obtain variance on PCs equal to their eigenvalue and so equal ');
fprintf('   to their initial proportion of variance\n\n');

Analysis_Charact.Standardization = input('   Rescaling on eigenvalue: ','s');                
if Analysis_Charact.Standardization == 'y'
   if not(exist(Eigenfile))
      [Eigenfile, pathname]=uigetfile('*.txt','Load TXT file: eigenvalues for rescaling');
      if filename~=0
            Eigenfile=[pathname Eigenfile];
            eigenvalues=load (Eigenfile);
         if isempty(eigenvalues)
            string=['File' Eigenfile ' is empty!!!.'];
            title='File Control Failed';
            errordlg(string,title,'modal');
         end; 
      end;
   else
      eigenvalues= load (Eigenfile);
   end;
   for yy=1:coll
      PCs(:,yy) = PCs(:,yy).*sqrt(eigenvalues(:,yy));
   end;
elseif Analysis_Charact.Standardization == 'n'
   PCs=PCs;
else
   PCs=PCs;
   warning('no standardization by eigenvalues is used');
end;
[l p] = size(PCs);

%Option2: Reduce number of PCs
Analysis_Charact.PCs_retain = input('\n   number of PCs selected: ');
if Analysis_Charact.PCs_retain>p
   Analysis_Charact.PCs_retain=p;
   fprintf('\n   !!! warning choice of PCs retained upper that number of PCs present\n '); 
   fprintf('  in the initial data matrix; analysis performs on this number of PCs present:\n '); 
   fprintf('  ==>> Number of max available PCs = %i\t',p); 
end;     
PCs = PCs(:,1:Analysis_Charact.PCs_retain); %New data matrix

%Option3: Choice of one type of analysis
fprintf('\n   Four analyses are available:\n');
fprintf('\n   (1) SGA: single group analysis');
fprintf('\n   (2) MGA: multiple group analysis');
fprintf('\n   (3) PDA: Partial Disparity Analysis');
fprintf('\n   (4) BTailTest: Bootstrap Tail Test\n\n');

analysis_perform = input('   choice: ');
if analysis_perform==1
   Analysis_Charact.analysis_perform='SGA';
   SGA(PCs,Analysis_Charact);
elseif analysis_perform==2
   Analysis_Charact.analysis_perform='MGA';
   MGA(PCs,Analysis_Charact);
elseif analysis_perform==3
   Analysis_Charact.analysis_perform='PDA';
   PDA(PCs,Analysis_Charact);
elseif analysis_perform==4
   Analysis_Charact.analysis_perform='BTailTest';
   BTailTest(PCs,Analysis_Charact);
end;
end;
%-----------------------------------------------------
function SGA(PCs,Analysis_Charact);
[sp npc]=size(PCs);
global Occurrencefile                % Default name of input file (temporal or geographical...occurrences)

%loading binary matrix of occurrence (temporal, geographical...occurrence) 
if not(exist(Occurrencefile))
   [filename, pathname]=uigetfile('*.txt','Load TXT file: Temporal, geographical...occurrence');
   if filename~=0
      filename=[pathname filename];
      occurrence=load (filename);
      if isempty(occurrence)
         string=['File' filename ' is empty!!!.'];
         title='File Control Failed';
         errordlg(string,title,'modal');
      end; 
   end;
else
   occurrence=load(Occurrencefile);
end;
[ind interv]=size(occurrence);
if ind~=sp
   error('    number of observation in occurrence file different that this in data_matrix file');
end;
min_obs=min(sum(occurrence));

%Option3: Rarefaction for some quantifiers
[Analysis_Charact.rmenu,samplesize_used]=Rarefaction_Menu(min_obs);

%Option4: Bootstrap menu Nb_Bootstrap_resampling
Analysis_Charact.Nb_Bootstrap_resampling = input('\n   Number of bootstrap resamplings: ');

%Perform Disparity analysis
res=main_analysis(PCs,occurrence,1,samplesize_used,interv,Analysis_Charact);

%Option5: Bootstrap estimator [mean standard deviation Confidence Interval (upper value lower value)]
CL = input('\n   Do you want a confidence interval (e.g. 95%) on the bootstrap y/n\n   case n: error bar = +/- 1 std\n   choice: ','s');
if CL=='n'
   Analysis_Charact.CI='+/- 1 stdev';
elseif CL=='y'
   CL = input('\n   input level of confidence interval (e.g., 95): ');
   CI=sprintf('%i',CL);
   Analysis_Charact.CI=strcat(CI,'%');
end;
for i=1:10
   [res.DispEstim(i).mean,res.DispEstim(i).std,res.DispEstim(i).lowval,res.DispEstim(i).uppval,res.DispEstim(i).lowrang,res.DispEstim(i).upprang]=bootstat(res.DispEstim(i).Bootsampl,CL);
end;
for i=1:2
   [res.Loc_Estim(i).mean,res.Loc_Estim(i).std,res.Loc_Estim(i).lowval,res.Loc_Estim(i).uppval,res.Loc_Estim(i).lowrang,res.Loc_Estim(i).upprang]=bootstat(res.Loc_Estim(i).Bootsampl,CL);
end;

%Plot Menu
%graph: disparity vs time
figDisp=figure( ...
        'Name','Disparity on N first PCs: mean and +/- one std of bootstrap sample', ...
        'Units','centimeters', ...
        'Position',[0.2 1 15 18]);
strg=char('Sum of univariates ranges','Root Product of univ. ranges','Range as max Eucl. Distance','AreaCHull !just 2 first PCs',...
      'PCO Volume','Sum of variances','Root Product of variances','Mean Distance Pairwise','Median Distance Pairwise','Mean Distance Centroid');
Dmin=0;
Tmin=0.5;
Tmax=interv+0.5;
for spt=1:10
   Dmax=max(res.DispEstim(spt).uppval)+ceil(max(res.DispEstim(spt).std)./2);
   Unit=Dmax./12;
   plotting(res.DispEstim(spt).mean,res.DispEstim(spt).lowrang,res.DispEstim(spt).upprang,strg(spt,:),spt,Tmin,Tmax,Dmin,Dmax,Unit);
end;
%graph: Pattern of occupation vs time
if npc>10
   nsubplt=10;
else
   nsubplt=npc;
end;
figMinMax=figure( ...
   'Name','Rarefied Min & Max on N first PCs (max of 10 PCs): mean and +/- one std of bootstrap sample',...
        'Units','centimeters', ...
        'Position',[10 1 15 18]);
for i=1:2
    res.Loc_Estim(i).mean=permute(res.Loc_Estim(i).mean,[3,2,1]);
    res.Loc_Estim(i).std=permute(res.Loc_Estim(i).std,[3,2,1]);
    res.Loc_Estim(i).lowval=permute(res.Loc_Estim(i).lowval,[3,2,1]);
    res.Loc_Estim(i).uppval=permute(res.Loc_Estim(i).uppval,[3,2,1]);
    res.Loc_Estim(i).lowrang=permute(res.Loc_Estim(i).lowrang,[3,2,1]);
    res.Loc_Estim(i).upprang=permute(res.Loc_Estim(i).upprang,[3,2,1]);
end;
for nsplt=1:nsubplt
    strg2=sprintf('Min & Max on PC%i',nsplt);
    MMmax=max(res.Loc_Estim(2).uppval(:,nsplt))+ceil(max(res.Loc_Estim(2).std(:,nsplt))./2);
    MMmin=min(res.Loc_Estim(1).lowval(:,nsplt))-ceil(max(res.Loc_Estim(1).std(:,nsplt))./2);
    Unit=MMmin+(MMmax-MMmin)./12;
    plotting(res.Loc_Estim(1).mean(:,nsplt),res.Loc_Estim(1).lowrang(:,nsplt),res.Loc_Estim(1).upprang(:,nsplt),strg2,nsplt,Tmin,Tmax,MMmin,MMmax,Unit);
    plotting(res.Loc_Estim(2).mean(:,nsplt),res.Loc_Estim(2).lowrang(:,nsplt),res.Loc_Estim(2).upprang(:,nsplt),strg2,nsplt,Tmin,Tmax,MMmin,MMmax,Unit);
end;

%Save Menu
fprintf('\n\n   ===================================\n');
Save_choice = input('   Save numerical results?\n   -> data are appended if file name alredy exists\n   (y/n): ','s');
if Save_choice=='y'
   Save_function(Analysis_Charact,1,interv,res);
elseif Save_choice=='n'
   repeat_choice = input('   Start again analysis(y) or exit (x)?: ','s');
   if repeat_choice=='x'
      stop=1;
      save stop
   end;
end;

%-----------------------------------------------------
function MGA(PCs,Analysis_Charact);
[sp npc]=size(PCs);
global Occurrencefile                % Default name of input file (temporal or geographical...occurrences)
global Groupfile                    % Default name of input file (group member)

%loading binary matrix of occurrence (temporal, geographical...occurrence) 
if not(exist(Occurrencefile))
   [filename, pathname]=uigetfile('*.txt','Load TXT file: Temporal, geographical...occurrence');
   if filename~=0
      filename=[pathname filename];
      occurrence=load (filename);
      if isempty(occurrence)
         string=['File' filename ' is empty!!!.'];
         title='File Control Failed';
         errordlg(string,title,'modal');
      end; 
   end;
else
   occurrence=load(Occurrencefile);
end;
[ind interv]=size(occurrence);
if ind~=sp
   error('    number of observation in occurrence file different that this in data_matrix file');
end;
%loading binary matrix of group member(taxinomic, geographic...groups) 
if not(exist(Groupfile))
   [filename, pathname]=uigetfile('*.txt','Load TXT file: Taxinomic, geographical...groups');
   if filename~=0
      filename=[pathname filename];
      group=load (filename);
      if isempty(group)
         string=['File' filename ' is empty!!!.'];
         title='File Control Failed';
         errordlg(string,title,'modal');
      end; 
   end;
else
   group=load(Groupfile);
end;
[ind ngrps]=size(group);
if ind~=sp
   error('    number of observation in occurrence file different that this in data_matrix file');
end;
%search minimum sample size 
min_obs=ind;
for G=1:ngrps  
   index_group = find(group(:,G));
   Goccur= occurrence(index_group,:);
   for j=1:interv  
      index_goccur = find(Goccur(:,j));
      min_ssize=length(index_goccur);
      if min_ssize>1 & min_ssize<min_obs
         min_obs=min_ssize;
      end;
   end;
end;

%Option3: Rarefaction for some quantifiers
[Analysis_Charact.rmenu,samplesize_used]=Rarefaction_Menu(min_obs);

%Option4: Bootstrap menu Nb_Bootstrap_resampling
Analysis_Charact.Nb_Bootstrap_resampling = input('\n   Number of bootstrap resamplings: ');

%Perform Disparity analysis on several groups
for G=1:ngrps  
   fprintf('\n\n   =====> studied group n°%i',G);
   fprintf(' on %i',ngrps);
   index_group = find(group(:,G));
   Gdata = PCs(index_group,:);
   Goccur= occurrence(index_group,:);
   res(G)=main_analysis(Gdata,Goccur,G,samplesize_used,interv,Analysis_Charact);
end;   

%Option5: Bootstrap estimator [mean standard deviation Confidence Interval (upper value lower value)]
CL = input('\n   Do you want a confidence interval (e.g. 95%) on bootstrap y/n\n   case n: error bar = +/- 1 std\n   choice: ','s');
if CL=='n'
   Analysis_Charact.CI='+/- 1 stdev';
elseif CL=='y'
   CL = input('\n   input level of Confidence interval (e.g., 95): ');
   CI=sprintf('%i',CL);
   Analysis_Charact.CI=strcat(CI,'%');
end;
for G=1:ngrps
   for i=1:10
   	[res(G).DispEstim(i).mean,res(G).DispEstim(i).std,res(G).DispEstim(i).lowval,res(G).DispEstim(i).uppval,res(G).DispEstim(i).lowrang,res(G).DispEstim(i).upprang]=bootstat(res(G).DispEstim(i).Bootsampl,CL);
	end;
	for i=1:2
   	[res(G).Loc_Estim(i).mean,res(G).Loc_Estim(i).std,res(G).Loc_Estim(i).lowval,res(G).Loc_Estim(i).uppval,res(G).Loc_Estim(i).lowrang,res(G).Loc_Estim(i).upprang]=bootstat(res(G).Loc_Estim(i).Bootsampl,'n');
	end;
end;

%Plot Menu
%graph: disparity vs time
for G=1:ngrps
strings=sprintf(' >>> group %i ',G);
figDisp=figure( ...
        'Name','Disparity on N first PCs: mean and +/- one std of bootstrap sample', ...
        'Units','centimeters', ...
        'Position',[0.2 1 15 18]);
   DispGroups= uicontrol('Parent',figDisp, ...
	'Units','centimeters', ...
	'FontName','Arial', ...
	'FontSize',10, ...
	'FontWeight','bold', ...
	'HorizontalAlignment','left', ...
	'ListboxTop',0, ...
	'Position',[0.2 0.2 2.2 0.5], ...
	'String',strings, ...
	'Style','text', ...
	'Tag','N°group');
strg=char('Sum of univariates ranges','Root Product of univ. ranges','Range as max Eucl. Distance','AreaCHull !just 2 first PCs',...
      'PCO Volume','Sum of variances','Root Product of variances','Mean Distance Pairwise','Median Distance Pairwise','Mean Distance Centroid');
Dmin=0;
Tmin=0.5;
Tmax=interv+0.5;
for spt=1:10
   Dmax=max(res(G).DispEstim(spt).uppval)+ceil(max(res(G).DispEstim(spt).std)./2);
   Unit=Dmax./12;
   plotting(res(G).DispEstim(spt).mean,res(G).DispEstim(spt).lowrang,res(G).DispEstim(spt).upprang,strg(spt,:),spt,Tmin,Tmax,Dmin,Dmax,Unit);
end;
%graph: Pattern of occupation vs time
if npc>10
   nsubplt=10;
else
   nsubplt=npc;
end;
figMinMax=figure( ...
   'Name','Rarefied Min & Max on N first PCs (max of 10 PCs): mean and +/- one std of bootstrap sample',...
        'Units','centimeters', ...
        'Position',[10 1 15 18]);
    MinMaxGroups= uicontrol('Parent',figMinMax, ...
	'Units','centimeters', ...
	'FontName','Arial', ...
	'FontSize',10, ...
	'FontWeight','bold', ...
	'HorizontalAlignment','left', ...
	'ListboxTop',0, ...
	'Position',[0.2 0.2 2.2 0.5], ...
	'String',strings, ...
	'Style','text', ...
	'Tag','N°group');
for i=1:2
    res(G).Loc_Estim(i).mean=permute(res(G).Loc_Estim(i).mean,[3,2,1]);
    res(G).Loc_Estim(i).std=permute(res(G).Loc_Estim(i).std,[3,2,1]);
    res(G).Loc_Estim(i).uppval=permute(res(G).Loc_Estim(i).uppval,[3,2,1]);
    res(G).Loc_Estim(i).lowval=permute(res(G).Loc_Estim(i).lowval,[3,2,1]);
    res(G).Loc_Estim(i).upprang=permute(res(G).Loc_Estim(i).upprang,[3,2,1]);
    res(G).Loc_Estim(i).lowrang=permute(res(G).Loc_Estim(i).lowrang,[3,2,1]);
end;
for nsplt=1:nsubplt
    strg2=sprintf('Min & Max on PC%i',nsplt);
    MMmax=max(res(G).Loc_Estim(2).uppval(:,nsplt))+ceil(max(res(G).Loc_Estim(2).std(:,nsplt))./2);
    MMmin=min(res(G).Loc_Estim(1).lowval(:,nsplt))-ceil(max(res(G).Loc_Estim(1).std(:,nsplt))./2);
    Unit=MMmin+(MMmax-MMmin)./12;
    plotting(res(G).Loc_Estim(1).mean(:,nsplt),res(G).Loc_Estim(1).lowrang(:,nsplt),res(G).Loc_Estim(1).upprang(:,nsplt),strg2,nsplt,Tmin,Tmax,MMmin,MMmax,Unit);
    plotting(res(G).Loc_Estim(2).mean(:,nsplt),res(G).Loc_Estim(1).lowrang(:,nsplt),res(G).Loc_Estim(1).upprang(:,nsplt),strg2,nsplt,Tmin,Tmax,MMmin,MMmax,Unit);
end;
end;
%Save Menu
fprintf('\n\n   ===================================\n');
Save_choice = input('   Save numerical results?\n   -> data are appended if file name already exists\n   (y/n): ','s');
if Save_choice=='y'
   Save_function(Analysis_Charact,ngrps,interv,res);
elseif Save_choice=='n'
   repeat_choice = input('   Start again analysis(y) or exit (x)?: ','s');
   if repeat_choice=='x'
      stop=1;
      save stop
   end;
end;

%-----------------------------------------------------
function PDA(PCs,Analysis_Charact);
[sp npc]=size(PCs);
global Occurrencefile                % Default name of input file (temporal or geographical...occurrences)
global Groupfile                    % Default name of input file (group member)

%loading binary matrix of occurrence (temporal, geographical...occurrence) 
if not(exist(Occurrencefile))
   [filename, pathname]=uigetfile('*.txt','Load TXT file: Temporal, geographical...occurrence');
   if filename~=0
      filename=[pathname filename];
      occurrence=load (filename);
      if isempty(occurrence)
         string=['File' filename ' is empty!!!.'];
         title='File Control Failed';
         errordlg(string,title,'modal');
      end; 
   end;
else
   occurrence=load(Occurrencefile);
end;
[ind interv]=size(occurrence);
if ind~=sp
   error('    number of observation in occurrence file different that this in data_matrix file');
end;
%loading binary matrix of group member(taxinomic, geographic...groups) 
if not(exist(Groupfile))
   [filename, pathname]=uigetfile('*.txt','Load TXT file: Taxinomic, geographical...groups');
   if filename~=0
      filename=[pathname filename];
      group=load (filename);
      if isempty(group)
         string=['File' filename ' is empty!!!.'];
         title='File Control Failed';
         errordlg(string,title,'modal');
      end; 
   end;
else
   group=load(Groupfile);
end;
[ind ngrps]=size(group);
if ind~=sp
   error('    number of observation in occurrence file different that this in data_matrix file');
end;

%initialize a global group (totality of data) and add to initial group data
global_group=ones(ind,1);
group=cat(2,global_group,group); 
[ind ngrps]=size(group);

%Option3: No rarefaction is needed for this analysis
Analysis_Charact.rmenu='none';

%Option4: Bootstrap menu Nb_Bootstrap_resampling
Analysis_Charact.Nb_Bootstrap_resampling = input('\n   Number of bootstrap resamplings: ');

res.Bootsampl=zeros(Analysis_Charact.Nb_Bootstrap_resampling,interv);
for j=1:interv 
   fprintf('\n   >>> studied sample n°%i',j);
   fprintf(' on %i',interv);
   index_occur = find(occurrence(:,j));
   Subdata = PCs(index_occur,:);
   Subgroups = group(index_occur,:);
   [subind col]=size(Subdata);
   for r=1:Analysis_Charact.Nb_Bootstrap_resampling
      [bootsampl bootGroups]=Gresampl(Subdata,Subgroups,subind);   %Bootstrap of subtemporal data and subtemporal group presence  
      M=mean(bootsampl);
      for G=1:ngrps
          index_grp = find(bootGroups(:,G));
          yy=length(index_grp);
          %lacking data
          if yy<1
             res(G).Bootsampl(r,j)=nan;
          elseif yy==1
             res(G).Bootsampl(r,j)=0;
          elseif yy>1
             G_bootsampl = bootsampl(index_grp,:);
             %Main analysis
             %The first group correspond to global group => Total disparity; others correspond to user specified => Partial Disparity other 
             res(G).Bootsampl(r,j)=sum(SEdist_centr(G_bootsampl,M))./(subind-1); %Mean square Euclidean Distance between specimens and centroid (Foote 1993)
          end;
          %Analysis parameters
          res(G).analysis_parameters(j,1)=subind;          %Sample size 
          res(G).analysis_parameters(j,2)=j;               %e.g. temporal subdivision (depends on nature of data contain in occurrence file)
          res(G).analysis_parameters(j,3)=G;               %e.g. group subdivision (depends on nature of data contain in group file; the first is a global group)
          %re-initialize for the next group 
          G_bootsampl=[];
          index_grp=[];
      end;
      %re-initialize for the next bootstrap resampling 
      bootsampl=[];
      bootGroups=[];
      M=[];
   end;
   %re-initialize for the next temporal interval 
   Subdata=[];
   Subgroups=[];
end;

%Option5: Bootstrap estimator [mean standard deviation Confidence Interval (upper value lower value)]
CL = input('\n   Do you want a confidence interval (e.g. 95%) on the bootstrap y/n\n   case n: error bar = +/- 1 std\n   choice: ','s');
if CL=='n'
   Analysis_Charact.CI='+/- 1 stdev';
elseif CL=='y'
   CL = input('\n   input level of confidence interval (e.g., 95): ');
   CI=sprintf('%i',CL);
   Analysis_Charact.CI=strcat(CI,'%');
end;
for G=1:ngrps
   for jj=1:interv
      Bsamp=res(G).Bootsampl(:,jj);
      indNan=isnan(Bsamp);
      BsampnoNan=Bsamp(find(indNan==0),:);
      if isempty(BsampnoNan)
         res(G).mean(jj)=nan;res(G).std(jj)=nan;res(G).lowval(jj)=nan;res(G).uppval(jj)=nan;res(G).lowrang(jj)=nan;res(G).upprang(jj)=nan;
      else
         [res(G).mean(jj),res(G).std(jj),res(G).lowval(jj),res(G).uppval(jj),res(G).lowrang(jj),res(G).upprang(jj)]=bootstat(BsampnoNan,CL);
      end;
      Bsampl=[];
      indNan=[];
      BsampnoNan=[];
   end;
end;
%Plot Menu
%graph: disparity vs time
strg=sprintf('Total Disparity');
for G=2:ngrps
   strg=char(strg,sprintf('Partial Disparity of group %i',G-1));
end;
nfig=ceil((ngrps)./10);
for nfg=1:nfig
    figDisp=figure( ...
            'Name','Disparity on N first PCs: mean and +/- one std of bootstrap sample', ...
            'Units','centimeters', ...
            'Position',[0.2 1 15 18]);
    Dmin=0;
    Tmin=0.5;
    Tmax=interv+0.5;
    grii=((nfg-1)*10)+1;
    if nfg==nfig
       grss=ngrps-((nfg-1)*10);
    else
       grss=nfg*10;
    end;
    for G=grii:grss
        Dmax=max(res(G).uppval)+ceil(max(res(G).std)./2);
        Unit=Dmax./12;
        plotting(res(G).mean,res(G).lowrang,res(G).upprang,strg(G,:),G,Tmin,Tmax,Dmin,Dmax,Unit);
    end;
end;

%Save Menu
fprintf('\n\n   ===================================\n');
Save_choice = input('   Save numeric results?\n   -> data are appended if file name already exists\n   (y/n): ','s');
if Save_choice=='y'
   Save_function(Analysis_Charact,ngrps,interv,res);
elseif Save_choice=='n'
   repeat_choice = input('   Start again analysis(y) or exit (x)?: ','s');
   if repeat_choice=='x'
      stop=1;
      save stop
   end;
end;

%-----------------------------------------------------
function BTailTest(PCs,Analysis_Charact);
[sp npc]=size(PCs);
global Groupfile                    % Default name of input file (group member)

fprintf('\n   Bootstrap Tail Test Help:\n');
fprintf('\n   >> The first group in group file is the one use for comparison');
fprintf('\n   This group is resampled with replacement (bootstrap)[except for one metric see below]');
fprintf('\n   and observed value of other groups are compared with this obtained distribution.');
fprintf('\n   Because disparity estimates can be increase or decrease, three cases arise:');
fprintf('\n     - one-tailed test with either upper tail or lower tail');
fprintf('\n     - two-tailed test.');
fprintf('\n\n   Range-based disparity estimates are sensitive to sample size and unavailable here. ');
fprintf('\n\n   1000 bootstrap resamplings is the minimum value for accuracy.\n');

%loading binary matrix of group member(taxinomic, geographic...groups) 
if not(exist(Groupfile))
   [filename, pathname]=uigetfile('*.txt','Load TXT file: group file with first column is the reference group');
   if filename~=0
      filename=[pathname filename];
      group=load (filename);
      if isempty(group)
         string=['File' filename ' is empty!!!.'];
         title='File Control Failed';
         errordlg(string,title,'modal');
      end; 
   end;
else
   group=load(Groupfile);
end;
[ind ngrps]=size(group);
if ind~=sp
   error('    number of observation in occurrence file different that this in data_matrix file');
end;

%Option: Choice of one disparity estimate => reduce time of analysis, permits more bootstrap replications
fprintf('\n   Choosing one disparity estimate saves time , allowing analysis with a higher number of');
fprintf('\n   bootstrap resamplings.');
fprintf('\n   Disparity estimate available:');
fprintf('\n     (1) Sum of variances');
fprintf('\n     (2) N-root product of variances');
fprintf('\n     (3) Mean pairwise Euclidean distance');
fprintf('\n     (4) Median pairwise Euclidean distance');
fprintf('\n     (5) Mean Euclidean distance to group centroid');
fprintf('\n     (6) Partial Disparity (resampled group is the complete data set)');
DispEstimate= input('\n\n   choice: ','s');
Analysis_Charact.DispEstimate=char(sprintf('Sum of variances'),sprintf('N-root product of variances'),...
   sprintf('Mean pairwise Euclidean distance'),sprintf('Median pairwise Euclidean distance'),...
   sprintf('Mean Euclidean distance to group centroid'),sprintf('Partial Disparity'));
switch DispEstimate
   case {'1'}
      Analysis_Charact.DispEstimate=Analysis_Charact.DispEstimate(1,:);
   case {'2'}
      Analysis_Charact.DispEstimate=Analysis_Charact.DispEstimate(2,:);
   case {'3'}
      Analysis_Charact.DispEstimate=Analysis_Charact.DispEstimate(3,:);
   case {'4'}
      Analysis_Charact.DispEstimate=Analysis_Charact.DispEstimate(4,:);
   case {'5'}
      Analysis_Charact.DispEstimate=Analysis_Charact.DispEstimate(5,:);
   case {'6'}
      Analysis_Charact.DispEstimate=Analysis_Charact.DispEstimate(6,:);
end;

%Option3: Rarefaction for some quantifiers
Analysis_Charact.rmenu='none';

%Option4: Bootstrap menu Nb_Bootstrap_resampling
Analysis_Charact.Nb_Bootstrap_resampling = input('\n   Number of bootstrap resamplings: ');
if Analysis_Charact.Nb_Bootstrap_resampling<1000
   disp('   number of resamplings is low: results may be inaccurate');
end;

%Option5: No statististic on bootstrap sample
Analysis_Charact.CI='none';

switch DispEstimate
   case {'1','2','3','4','5'}
%Perform Disparity Analysis on first group
index_group = find(group(:,1));
Gdata = PCs(index_group,:);
[lig col]=size(Gdata);
step=1;
mstep=Analysis_Charact.Nb_Bootstrap_resampling./10;
bound=0:mstep:Analysis_Charact.Nb_Bootstrap_resampling;
bound(1)=1;
for r=1:Analysis_Charact.Nb_Bootstrap_resampling
   if r==bound(step)
      fprintf('\n   >>> bootstrap iteration n°%i',bound(step));
      fprintf(' on %i',Analysis_Charact.Nb_Bootstrap_resampling);
      step=step+1;
   end;
      bootsampl=resamplk(Gdata,lig);             %bootstrap sample 
      %Dispersion estimate
   switch DispEstimate
      case {'1'}
      res.Bootsampl(r,1)=sum(var(bootsampl));                     %Sum of univariate variance (or Total Variance) (e.g. Wills et 1994, Foote XXX, Eble 2000 Paleobiology)
      case {'2'}
      res.Bootsampl(r,1)=(prod(var(bootsampl),2)).^(1/Analysis_Charact.PCs_retain);       %root N of product of univariate variance (e.g. Wills et 1994, Foote XXX Paleobiology)
      case {'3'}
      pairwdist=(sqrt(sqrtpairwdist(bootsampl))); 
      res.Bootsampl(r,1)=mean(pairwdist);                                %Mean pairwise Euclidean distance (e.g. Foote 1999 Paleobiology)
      case {'4'}
      pairwdist=(sqrt(sqrtpairwdist(bootsampl))); 
      res.Bootsampl(r,1)=median(pairwdist);                              %Median pairwise Euclidean distance (Wagner 1997 Paleobiology)
      case {'5'}
      res.Bootsampl(r,1)=mean(Edistcentroid(bootsampl));          %Mean Euclidean distance between specimens and centroid
   end;
      %re-initialize
      clear pairwdist;
end;
case {'6'}
   %Partial Disparity (group resampling in the complet data set)
[lig col]=size(PCs);
step=1;
mstep=Analysis_Charact.Nb_Bootstrap_resampling./10;
bound=0:mstep:Analysis_Charact.Nb_Bootstrap_resampling;
bound(1)=1;
for r=1:Analysis_Charact.Nb_Bootstrap_resampling
   if r==bound(step)
      fprintf('\n   >>> bootstrap iteration n°%i',bound(step));
      fprintf(' on %i',Analysis_Charact.Nb_Bootstrap_resampling);
      step=step+1;
   end;
      bootsampl=resamplk(PCs,lig);             %bootstrap sample 
      M=mean(bootsampl);
      res.Bootsampl(r,1)=sum(SEdist_centr(bootsampl,M))./(lig-1); %Mean square Euclidean Distance between specimens and centroid (Foote 1993)
      M=[];
   end;
  M_obs=mean(PCs); %mean morphology observed use for estimate observed value   
end;

%Observed values for disparity estimates for for group 1 to N 
for ng=1:ngrps
      index_group = find(group(:,ng));
      ii=length(index_group);
      Gdata = PCs(index_group,:);
      res.analysis_parameters(ng,1)=ng;                                              %index of groups with 1 = group used for comparison 
      %Dispersion estimate
         fprintf('\n   >>> studied group n°%i',ng);
         fprintf(' on %i',ngrps);
   switch DispEstimate
      case {'1'}
      res.obs_value(ng,1)=sum(var(Gdata));                                             %Sum of univariate variance (or Total Variance) (e.g. Wills et 1994, Foote XXX, Eble 2000 Paleobiology)
      case {'2'}
      res.obs_value(ng,1)=(prod(var(Gdata),2)).^(1/Analysis_Charact.PCs_retain);       %root N of product of univariate variance (e.g. Wills et 1994, Foote XXX Paleobiology)
      case {'3'}
      pairwdist=(sqrt(sqrtpairwdist(Gdata))); 
      res.obs_value(ng,1)=mean(pairwdist);                                %Mean pairwise Euclidean distance (e.g. Foote 1999 Paleobiology)
      pairwdist=[];
      case {'4'}
      pairwdist=(sqrt(sqrtpairwdist(Gdata))); 
      res.obs_value(ng,1)=median(pairwdist);                              %Median pairwise Euclidean distance (Wagner 1995 Paleobiology)
      pairwdist=[];
      case {'5'}
      res.obs_value(ng,1)=mean(Edistcentroid(Gdata));                     %Mean Euclidean distance between specimens and centroid
      case {'6'}
      res.obs_value(ng,1)=sum(SEdist_centr(Gdata,M_obs))./(lig-1); %Mean square Euclidean Distance between specimens and centroid (Foote 1993)
   end;
      %re-initialize
      clear pairwdist;
end;

%find probabilities for each descriptor for group 2 to N-1 (suppress group 1 => comparison)
for ng=1:ngrps
      M=mean(res.Bootsampl(:,1));
      index_valueBT = find(abs(res.Bootsampl(:,1)-M) >= abs(res.obs_value(ng,1)-M));
      numbersBT = length(index_valueBT);
      res.p_twotail(ng,1) = numbersBT./Analysis_Charact.Nb_Bootstrap_resampling;
      index_valueUTS = find(res.Bootsampl(:,1) >= res.obs_value(ng,1));
      numbersUTS = length(index_valueUTS);
      res.p_upperonetail(ng,1) = numbersUTS./Analysis_Charact.Nb_Bootstrap_resampling;
      index_valueUTI = find(res.Bootsampl(:,1) <= res.obs_value(ng,1));
      numbersUTI = length(index_valueUTI);
      res.p_loweronetail(ng,1) = numbersUTI./Analysis_Charact.Nb_Bootstrap_resampling;
      
      %re-initialize
      M=[];
      index_valueBT=[]; 
      numbersBT=[];
      index_valueUTS=[];
      numbersUTS=[];
      index_valueUTI=[];
      numbersUTI=[];
end;

%Plot Menu
%graph: Histogram of bootsample of disparity estimate
min_val=min(res.Bootsampl);
min_obs=min(res.obs_value);
if min_obs<min_val
   min_val=min_obs;
end;
max_val=max(res.Bootsampl);
max_obs=max(res.obs_value);
if max_obs>max_val
   max_val=max_obs;
end;
strg=strcat(Analysis_Charact.DispEstimate,sprintf(' (Bootstrap resampling = %i',Analysis_Charact.Nb_Bootstrap_resampling),...
        sprintf(')'));
figDisp=figure( ...
        'Name','Disparity on N first PCs: mean and +/- one std of bootstrap sample', ...
        'Units','pixels', ...
        'Position',[28 302 342 422]);
	  hist(res.Bootsampl,30);
     xlabel(strg,'FontName','Geneva','FontSize',9,'FontWeight','bold');
     hold on;
%Save Menu
fprintf('\n\n   ===================================\n');
Save_choice = input('   Save numerical results?\n   -> data are appended if file name exists\n   (y/n): ','s');
if Save_choice=='y'
   Save_function(Analysis_Charact,ngrps,1,res);
elseif Save_choice=='n'
   repeat_choice = input('   Start again analysis(y) or exit (x)?: ','s');
   if repeat_choice=='x'
      stop=1;
      save stop
   end;
end;

%--------------------------------------------------------
function res=main_analysis(data,Occur,groupindex,samplesize_used,interv,Analysis_Charact);
res.analysis_parameters=cat(2,zeros(interv,3),groupindex*ones(interv,1));
for j=1:interv  
   fprintf('\n   >>> studied sample n°%i',j);
   fprintf(' on %i',interv);
   indoccur = find(Occur(:,j));
   subdata = data(indoccur,:);
   [lig col]=size(subdata);
   if lig>1
   if samplesize_used=='obs'
      Rsamplesize=lig;
   elseif samplesize_used>lig
      Rsamplesize=lig;
      fprintf('\n   for one sample, rarefaction size is upper to sample size:\n   analysis is performed with sample size !!\n ');   
   else
      Rsamplesize=samplesize_used;
   end;
   for r=1:Analysis_Charact.Nb_Bootstrap_resampling
         bootsampl=resamplk(subdata,lig);             %bootstrap sample 
      if samplesize_used~=lig
         bootsampk=resamplk(subdata,Rsamplesize);     %bootstrap sample rarefied to k specimens
      else
         bootsampk=bootsampl;
      end;
      %Dispersion estimate
      res.DispEstim(1).Bootsampl(r,j)=sum(unirange(bootsampk));                %Sum of univariate range (or Total Variance) (e.g. Wills et 1994, Foote 1992, Eble 2000 Paleobiology)
      res.DispEstim(2).Bootsampl(r,j)=prod(unirange(bootsampk),2).^(1/Analysis_Charact.PCs_retain);    %root N of product of univariate range (e.g. Wills et 1994, Foote 1992, Paleobiology)
      res.DispEstim(3).Bootsampl(r,j)=sqrt(max(sqrtpairwdist(bootsampk)));     %Range estimate as the maximum distance Euclidean pairwise (Ciampaglio et al. 2001 Paleobiology)
      res.DispEstim(4).Bootsampl(r,j)=convexhull(bootsampk);                   %Area of the convex hull (e.g. Foote 1999 Paleobiology, Navarro et al in prep  !! warning just in 2 dimensions)
      res.DispEstim(5).Bootsampl(r,j)=PCOvolume(bootsampl);                    %PCO volume (Ciampaglio et al. 2001)
      res.DispEstim(6).Bootsampl(r,j)=sum(var(bootsampl));                     %Sum of univariate variance (or Total Variance) (e.g. Wills et 1994, Foote XXX, Eble 2000 Paleobiology)
      res.DispEstim(7).Bootsampl(r,j)=(prod(var(bootsampl),2)).^(1/Analysis_Charact.PCs_retain);       %root N of product of univariate variance (e.g. Wills et 1994, Foote XXX Paleobiology)
      pairwdist=(sqrt(sqrtpairwdist(bootsampl))); 
      res.DispEstim(8).Bootsampl(r,j)=mean(pairwdist);                                %Mean pairwise Euclidean distance (e.g. Foote 1999 Paleobiology)
      res.DispEstim(9).Bootsampl(r,j)=median(pairwdist);                              %Median pairwise Euclidean distance (Wagner 1997 Paleobiology)
      pairwdist=[];
      res.DispEstim(10).Bootsampl(r,j)=mean(Edistcentroid(bootsampl));         %Mean Euclidean distance between specimens and centroid
      %Location estimated
      res.Loc_Estim(1).Bootsampl(r,:,j)=min(bootsampk);                        %Rarefied minimum value of PC (McShea 1994 Evolution; Navarro et al. in prep) 
      res.Loc_Estim(2).Bootsampl(r,:,j)=max(bootsampk);                        %Rarefied maximum value of PC (McShea 1994 Evolution; Navarro et al. in prep)
      
      %re-initialize
      bootsampl=[];
      bootsampk=[];
   end;
   res.analysis_parameters(j,1)=lig;             %Sample size
   res.analysis_parameters(j,2)=Rsamplesize;     %Rarefied sample size
   res.analysis_parameters(j,3)=j;               %e.g. temporal subdivision (depends on nature of data contain in occurrence file)
   
   %re-initialize for the next temporal subdivision 
   Rsamplesize=samplesize_used;
   %lacking data
   elseif lig<1
      for i=1:10
         res.DispEstim(i).Bootsampl(1,j)=nan; 
      end;
      for i=1:2
         res.Loc_Estim(i).Bootsampl(1,:,j)=nan*ones(1,col);
      end;
   elseif lig==1
      for i=1:10
         res.DispEstim(i).Bootsampl(1,j)=nan; 
      end;
      for i=1:2
         res.Loc_Estim(i).Bootsampl(1,:,j)=subdata;
      end;
   end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-----------------------------------------------------
% subfunctions for bootstrap resampling 
%-----------------------------------------------------
function Dboot=resamplk(data,k)
[l p]=size(data);
Index_random=Intrndunif(l,k,1);
Dboot=data(Index_random,:); %bootstrap resampling of k specimens from l

%------------------------------------------------------
function [out1,out2]=Gresampl(data,grp,k)
%used for Partial disparity analysis
[l p]=size(data);
Index_random=Intrndunif(l,k,1);
out1=data(Index_random,:);
out2=grp(Index_random,:);

%------------------------------------------------------
function indrandom=Intrndunif(a,b,c)
indunif=rand(b,c);
rescalind=a.*indunif;
indrandom=ceil(rescalind);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------------------------------------------------------
% subfunctions for disparity estimates
%------------------------------------------------------
function variance=var(X)
StandDev=std(X);
variance=StandDev.^2;

%------------------------------------------------------
function unrang=unirange(X)
Minima=min(X);
Maxima=max(X);
unrang=Maxima-Minima;

%------------------------------------------------------
function spdist=sqrtpairwdist(data)
%Use vectorization method using code of Roland Bunschoten
%available to http://www.mathsworks.com
%This fast method permit to increase number of
%bootstrap resampling
%This method works really quickly with large sample size 
%as used with BTailTest comparing to the previous used code.

spdist=distance(data',data');

%Extraction of lower triangular portion of the squared distance matrix
%without the diagonal based on trilow.m function of R. E. Strauss
%available to http://www.biol.ttu.edu/Strauss/Matlab/matlab.htm
spdist=trilow(spdist);

%------------------------------------------------------
function d = distance(a,b)
% DISTANCE - computes Euclidean distance matrix
% E = distance(A,B)
%    A - (DxM) matrix 
%    B - (DxN) matrix
% Returns:
%    E - (MxN) Euclidean distances between vectors in A and B
% Description : 
%    This fully vectorized (VERY FAST!) m-file computes the 
%    Euclidean distance between two vectors by:
%                 ||A-B|| = sqrt ( ||A||^2 + ||B||^2 - 2*A.B )
% Example : 
%    A = rand(400,100); B = rand(400,200);
%    d = distance(A,B);
% Author   : Roland Bunschoten
%            University of Amsterdam
%            Intelligent Autonomous Systems (IAS) group
%            Kruislaan 403  1098 SJ Amsterdam
%            tel.(+31)20-5257524
%            bunschot@wins.uva.nl
% Last Rev : Oct 29 16:35:48 MET DST 1999
% Tested   : PC Matlab v5.2 and Solaris Matlab v5.3
% Thanx    : Nikos Vlassis
% Copyright notice: You are free to modify, extend and distribute 
%    this code granted that the author of the original code is 
%    mentioned as the original author of the code.

if (nargin ~= 2)
   error('Not enough input arguments');
end

if (size(a,1) ~= size(b,1))
   error('A and B should be of same dimensionality');
end

aa=sum(a.*a,1); bb=sum(b.*b,1); ab=a'*b; 
%d = sqrt(abs(repmat(aa',[1 size(bb,2)]) + repmat(bb,[size(aa,2) 1]) - 2*ab));
%Modifications for MDA: computing squared distance:
d = abs(repmat(aa',[1 size(bb,2)]) + repmat(bb,[size(aa,2) 1]) - 2*ab);

%----------------------------------
function [c,i,j] = trilow(x)
% TRILOW: Extracts the lower triangular portion (without diagonal) of a
%         square symmetric matrix into a columnwise column vector.
%         Optionally returns corresponding subscripts.
%         If given a scalar, returns it.
%         Use trisqmat() to reverse the extraction.
%     Usage: [c,i,j] = trilow(x)
% RE Strauss, 1/13/96
%   8/20/99 - subscripting changes due to new Matlab v5 conventions.
  [n,p] = size(x);
  if (n~=p)
    error('  Matrix must be square');
  end;
  if (n==1)
    c = x;
  else
    t = tril(ones(n)) - eye(n);
    t = t(:);
    x = x(:);
    c = x(t==1);
  end;
  if (nargout>1)
    i = ones(length(c),1);
    j = i;
    k = 0;
    for jj = 1:(n-1)
      for ii = (jj+1):n
        k = k+1;
        i(k) = ii;
        j(k) = jj;
      end;
    end;
  end;
  %return;

%-------------------------------------------------------
function Edistcentroid=Edistcentroid(X)
[l p]=size(X);
M=mean(X);
Differ=X-M(ones(l,1),:);
differ2=Differ.^2;              
Edistcentroid=sqrt(sum((differ2),2));    

%-------------------------------------------------------
function SEdist_centr=SEdist_centr(X,M)
%use in Partial disparity analysis
[l p]=size(X);
Differ=X-M(ones(l,1),:);
differ2=Differ.^2;              
SEdist_centr=sum((differ2),2); 

%-------------------------------------------------------
function PCOvolume=PCOvolume(X)
%see Ciampaglio et al. 2001 Paleobiology
%Standardization by root square of number of specimens
%gives same results that utilization of rarefaction procedure
%without standardization by square of number of specimens
%=> Initial formulae of Ciampaglio are conserved here.
[l p]=size(X);
for i=1:l
   for j=i:l
      Cross_dist_Matrix(i,j)=sqrt(sum((X(i,:)-X(j,:)).^2));
      Cross_dist_Matrix(j,i)=Cross_dist_Matrix(i,j);
   end;
end;
[u, s, v]=svd(Cross_dist_Matrix);
eigen=sort(diag(s,0));
lambda1=eigen(l);
lambda2=eigen(l-1);
PCOvolume=lambda1*lambda2/(l.^2);

%-------------------------------------------------------
function areaconvexhull=convexhull(data)
%search convex hull in 2D and compute area
data_2D=data(:,1:2);
[l p]=size(data);
Nmin=searchdif(data);
%because convex hull in 2D, 2+1 points are necessary
if Nmin>2
   CHindex = convhull(data(:,1),data(:,2));
   convexhull=data_2D(CHindex,:);
   areaconvexhull=polyarea(convexhull(:,1),convexhull(:,2));
else
   areaconvexhull=0;
end;

%-------------------------------------------------------
function indexdif=searchdif(data)
indexdif=0;
aa=0;
while aa<3
   dif=[];
   [uu pp]=size(data);
   if uu==0
      break;
   end;
   kk=0;
   for jj=2:uu
      if data(1,:)~=data(jj,:)
         kk=kk+1;
         dif(kk,:)=data(jj,:);
      end;
   end;
   if kk~=0
      indexdif=indexdif+1;
   else
      break;
   end;
   aa=indexdif;
   data=dif;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------------------------------------------------------
% subfunctions for bootstrap statistic (mean std CL)
%------------------------------------------------------
function [Mean,stdev,lowval,uppval,LowRang,UppRang]=bootstat(data,CL)
[l p]=size(data);
Mean=mean(data);
stdev=std(data);
if CL=='n'
   uppval=Mean+stdev;
   lowval=Mean-stdev;
   UppRang=stdev;
   LowRang=stdev;
else
   [lowval, uppval]=ConfInterv(data, CL);
   UppRang=uppval-Mean;
   LowRang=Mean-lowval;
end;
%------------------------------------------------------
function [Lower, Upper]=ConfInterv(data, CL)
%Confidence interval on bootstrap estimator based on 
%percentile method
[re ax int]=size(data);
val=(((100-CL)./100)./2);
q1=ceil(re*val);
q2=ceil(re-re*val);
for ii=1:ax
   for j=1:int
       ranked=sort(data(:,ii,j));
       Upper(1,ii,j)=ranked(q2);
       Lower(1,ii,j)=ranked(q1);
   end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------------------------------------   
%Plotting functions
%-------------------------------------------------------
function plotting(Yvalue,lower,upper,leg,nsplot,Xmin,Xmax,Ymin,Ymax,Unit)
X=1:length(Yvalue);
subplot(5,2,nsplot);
sbp=errorbar(X,Yvalue,lower,upper,'k'); 
set(sbp,'LineWidth',1);
axis([Xmin Xmax Ymin Ymax])
hold on;
text(Xmin+0.05,Unit,leg);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------------------------------------------------------
%Rarefaction Menu
%------------------------------------------------------
function [rmenu,samplesize_used]=Rarefaction_Menu(min_obs)
fprintf('\n   Rarefaction Menu:\n');
fprintf('\n     Help:');
fprintf('\n     Some estimate of morphospace occupation are sensitive to sample size:');
fprintf('\n\t-Total range,\n\t-N Root of Product of ranges,\n\t-Range,\n\t-Area of convex hull,\n\t-Min & Max on PC  ');
fprintf('\n     Variance-based estimates are not sensitive to sample size and their analyses');
fprintf('\n     are not based on rarefied sample');
fprintf('\n     The rarefaction procedure is used to eliminated this link (Foote 1992 Paleobiology)');
fprintf('\n     The rarefaction is performed using bootstrapping (Eble 2000 Paleobiology).\n');
fprintf('\n   (1) no rarefaction');
fprintf('\n   (2) rarefaction with minimum observed sample size = %i\t',min_obs);
fprintf('\n   (3) rarefaction with a user-choice sample size\n');
rarefaction = input('\n   choice:  ');
if rarefaction == 1
   samplesize_used = 'obs';
   rmenu=sprintf('no rarefaction: size used: obs');
elseif rarefaction == 2
   samplesize_used = min_obs;  
   rmenu=sprintf('y: rarefaction size used: %i',min_obs);
elseif rarefaction == 3
   samplesize_used = input('\n   sample size of rarefaction: ');  
   rmenu=sprintf('y: rarefaction size used: %i',samplesize_used);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------------------------------------
%Save as:
%-------------------------------------------------------
function Save_function(Analysis_Charact,ngroup,interv,res);
%erase files if necessary
[file, path] = uiputfile( '.txt', 'Save results as:');
outfile=fullfile(path, file);
%Analysis Characteristic
Summary_string=strcat(sprintf('\n\n   Analysis Summary:'),...
   sprintf('\n   ==================================='),...
   sprintf('\n   Analysis performed: %s',Analysis_Charact.analysis_perform),...
   sprintf('\n   Number of PCs selected: %i',Analysis_Charact.PCs_retain),...
   sprintf('\n   PCs rescaling to eigenvalues: %s',Analysis_Charact.Standardization),...
   sprintf('\n   Rarefaction: %s',Analysis_Charact.rmenu),...
   sprintf('\n   Number of bootstrap resamplings: %i',Analysis_Charact.Nb_Bootstrap_resampling),...
   sprintf('\n   Upper-lower values: %s',Analysis_Charact.CI));
switch Analysis_Charact.analysis_perform
   case {'BTailTest'}
      Summary_string=strcat(Summary_string,sprintf('\n   Disparity metric used: %s',Analysis_Charact.DispEstimate));
end;
%Print Analysis Characteristics
disp(Summary_string);
Summary_outfile=strcat(Summary_string,sprintf('\n\n')); 
fid = fopen(outfile,'a');
fprintf(fid,Summary_outfile,'uint8');
fclose(fid);
switch Analysis_Charact.analysis_perform
   case {'SGA','MGA'}
   %format data
   string_outfile=strcat(sprintf('\n\nSample Size'),...
      sprintf('\tRarefaction Size'),...
      sprintf('\tInterval'),...
      sprintf('\tGroup'),...
      sprintf('\tSum of Ranges'),sprintf('\tstdev SR'),sprintf('\tupper value SR'),sprintf('\tlower value SR'),...
      sprintf('\tR Product of Ranges'),sprintf('\tstdev PR'),sprintf('\tupper value PR'),sprintf('\tlower value PR'),...
      sprintf('\tRange (Max Dist.)'),sprintf('\tstdev R'),sprintf('\tupper value R'),sprintf('\tlower value R'),...
      sprintf('\tArea Conv. Hull'),sprintf('\tstdev ACH'),sprintf('\tupper value ACH'),sprintf('\tlower value ACH'),...
      sprintf('\tPCO volume'),sprintf('\tstdev PCOvol'),sprintf('\tupper value PCOvol'),sprintf('\tlower value PCOvol'),...
      sprintf('\tSum of Variances'),sprintf('\tstdev SV'),sprintf('\tupper value SV'),sprintf('\tlower value SV'),...
      sprintf('\tR Product of Variances'),sprintf('\tstdev PV'),sprintf('\tupper value PV'),sprintf('\tlower value PV'),...
      sprintf('\tMean Pairw. Dist'),sprintf('\tstdev MPD'),sprintf('\tupper value MPD'),sprintf('\tlower value MPD'),...
      sprintf('\tMedian Pairw. Dist'),sprintf('\tstdev MdPD'),sprintf('\tupper value MdPD'),sprintf('\tlower value MdPD'),...
      sprintf('\tMean Dist. Centroid'),sprintf('\tstdev MDc'),sprintf('\tupper value MDc'),sprintf('\tlower value MDc'));

   col=Analysis_Charact.PCs_retain;
   if col>10
      col=10;%limit Min Max values to 10 first PCs
   end;
   for Naxis=1:col
      Min_string=sprintf('Min%i',Naxis);
      Min_string=strcat(sprintf('\t%s',Min_string),sprintf('\tstdev%s',Min_string),sprintf('\tupper value%s',Min_string),sprintf('\tlower value%s',Min_string));
      Max_string=sprintf('Max%i',Naxis);
      Max_string=strcat(sprintf('\t%s',Max_string),sprintf('\tstdev%s',Max_string),sprintf('\tupper value%s',Max_string),sprintf('\tlower value%s',Max_string));
	   string_outfile=strcat(string_outfile,Min_string,Max_string);
   end

   for G=1:ngroup
      for Naxis=1:col
         for i=1:2
             temp(i).mm=cat(2,res(G).Loc_Estim(i).mean(:,Naxis),res(G).Loc_Estim(i).std(:,Naxis),res(G).Loc_Estim(i).uppval(:,Naxis),res(G).Loc_Estim(i).lowval(:,Naxis));
         end;
      temp(Naxis).Location=cat(2,temp.mm);
      end;
      tmp(G).Loc=cat(2,temp.Location);
   end;
   Location=cat(1,tmp.Loc);
   [uu vv]=size(Location);
   for G=1:ngroup
      for i=1:10
          for nint=1:interv
              temp_result(i).dd(nint,:)=cat(2,res(G).DispEstim(i).mean(nint),res(G).DispEstim(i).std(nint),res(G).DispEstim(i).uppval(nint),res(G).DispEstim(i).lowval(nint));
          end;
          tmp(G).result=cat(2,temp_result.dd);
      end;
   end;
   all_result=cat(1,tmp.result);
   res_analysis_parameters=cat(1,res.analysis_parameters);
   [uuu vvv]=size(all_result);
   for nint=1:interv*ngroup
      string_res=sprintf('\n%i\t%i\t%i\t%i',res_analysis_parameters(nint,:));

      string_outfile=strcat(string_outfile,string_res);
      all_result_temp=all_result(nint,:);
      for i=1:vvv
         string_outfile=strcat(string_outfile,sprintf('\t%f',all_result_temp(i)));
      end;  
      temp_MinMax=Location(nint,:);
      for i=1:vv
         string_outfile=strcat(string_outfile,sprintf('\t%f',temp_MinMax(i)));
      end;
   end;
   string_outfile=strcat(string_outfile,sprintf('\n\n')); 
        
   fid = fopen(outfile,'a');
   fprintf(fid,string_outfile,'uint8');
   fclose(fid);
   
   repeat_choice = input('\n   Start again analysis(y) or exit (x)?: ','s');
   if repeat_choice=='x'
      stop=1;
      save stop
   end;
   
   case {'PDA'}
   %format data
   string_outfile=strcat(sprintf('\n\nSample Size'),...
      sprintf('\tInterval'),...
      sprintf('\tGroup'),...
      sprintf('\tMean Square Eucl. D. centroid'),...
      sprintf('\tstdev'),...
      sprintf('\tupper value'),...
      sprintf('\tlower value'));
   for G=1:ngroup
      for nint=1:interv
          tmp(G).result(nint,:)=cat(2,res(G).mean(nint),res(G).std(nint),res(G).uppval(nint),res(G).lowval(nint));
      end;
   end;
   all_result=cat(1,tmp.result);
   res_analysis_parameters=cat(1,res.analysis_parameters);
   [uuu vvv]=size(all_result);
   for nint=1:interv*ngroup
      string_res=sprintf('\n%i\t%i\t%i',res_analysis_parameters(nint,:));
      string_outfile=strcat(string_outfile,string_res);
      all_result_temp=all_result(nint,:);
      for i=1:vvv
         string_outfile=strcat(string_outfile,sprintf('\t%f',all_result_temp(i)));
      end;  
   end;
   string_outfile=strcat(string_outfile,sprintf('\n\n')); 
        
   fid = fopen(outfile,'a');
   fprintf(fid,string_outfile,'uint8');
   fclose(fid);
   
   repeat_choice = input('\n   Start again analysis(y) or exit (x)?: ','s');
   if repeat_choice=='x'
      stop=1;
      save stop
   end;
   
   case {'BTailTest'}
   %format data
   string_outfile=strcat(sprintf('\n\nGroup'),...
      sprintf('\tObserved value'),...
      sprintf('\tp two-tail'),sprintf('\tp upper one-tail'),sprintf('\tp lower one-tail'));
   string_outfile=strcat(string_outfile,sprintf('\n'));
   all_result=cat(2,res.obs_value(:,1),res.p_twotail(:,1),res.p_upperonetail(:,1),res.p_loweronetail(:,1));
   [uuu vvv]=size(all_result);
   for ng=1:ngroup
      string_res=sprintf('\n%i',res.analysis_parameters(ng,1));
      string_outfile=strcat(string_outfile,string_res);
      all_result_temp=all_result(ng,:);
      for i=1:vvv
         string_outfile=strcat(string_outfile,sprintf('\t%f',all_result_temp(i)));
      end;  
   end;
   string_outfile=strcat(string_outfile,sprintf('\n\n')); 
   fprintf('\n   -----------------------------------\n');   
   disp(string_outfile);
   fid = fopen(outfile,'a');
   fprintf(fid,string_outfile,'uint8');
   fclose(fid);
   
   repeat_choice = input('\n   Start again analysis(y) or exit (x)?: ','s');
   if repeat_choice=='x'
      stop=1;
      save stop
   end;
end;

%*******************************************************
