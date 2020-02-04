function [flag,olap] = landpoints(lon,lat,file,part,j,plotall,crit)
% LANDPOINTS	Test for on-land lon/lat positions.
% Uses both a provided 6 minute mask and fine scale vector coastline
% data to flag positions on land.
% 
% [flag,olap] = landpoints(lon,lat,file,part,j,plotall,crit)
% 
% lon,lat = Numeric input for positions to test. The output flag will
%           have the same size as these. For large data sets it is
%           recommended to reduce precision to single.
% file    = String with full path and name of file the positions come
%           from, in order to put and name the output mat-files
%           accordingly.
%           (default = '~/Downloads/file'; can be edited below).
% part    = Three element logical, for which of the three parts of
%           the test (se below) to run (default = logical([1 1 0]);
%           full testing). Plotting test maps cannot be done in the
%           same call as the test. Hence, [0 0 1] is the correct
%           input for making plots for inspection of results. 
% j       = Integer used as numbering of output messages and figures.
%           Useful when processing several files of data (default=0).
% plotall = Logical, set as true to force maps to be plotted also when
%           everything is fine (default is to only plot maps when
%           there are points on land). 
% crit    = Setting for optimization of clustering of Part 1. These
%           are empirically chosen using large and global sets of
%           positions, so do not input if you are not experienced
%           in the use of this function.
%
% flag    = Logical the same size as input lon/lat, true for on-land
%           positions. These are also stored in an output mat-file.
% olap    = Numerical of value one for clusters with coastline
%           polygions that overlap eachother, which are normally
%           lakes, but worthwhile checking out. 
%
%
% THE THREE PARTS:
%
% There are three parts that can be run independently, but in sequence
% and with the exact same dataset. Default is to run the test and
% provide the flags (Part 1 and 2).
%
% 1) First, coastal and inland masks at 6 arcminutes are used to
% eliminate open ocean positions, and flag clearly inland
% positions. This mask is provided with this function in
% coastal_and_inland_masks6.mat. Then, a coastline file for the
% remaining coastal-zone positions is made, using the full resolution
% GSHHS coastline. If positions are widespread or many, clustering of
% positions will be done to reduce load. This making of test data is a
% time consuming process, hence separated from the actual testing of
% positions.
%
% 2) Test the postitions against the corresponding coastline
% (clusters). INPOLYGON is used to find data inside the polygons. If a
% position falls inside two polygons, it's likely in a lake. The clearly
% inland positions identified in Part 1 are also assigned to the flag
% variable. On the other hand, on-land positions closer to the coastline
% than half the resolution of the nearest coastline points, are not
% flagged.
%
% 3) Maps for visual inspection. (This part must be run separately.)
% Prints maps for each cluster file and shows red dots for flagged
% positions and green dots for not flagged. Uses input positions,
% coastlines from (1), and flags from (2). By default plots maps only
% for clusters with flagged positions, but can be set to plot for all
% positions (see above).
%
%
% OUTPUT FILES:
%
% (from part 1) Mat-files with map polygons and indices to which data
% belongs in each cluster, for each data
% file and cluster. Written alongside data files and named 
% '<file>.landtest.coast.cluster*.mat'. 
%
% (from part 1) Mat-files with lon and lat ranges, indices to positions
% in the costal zone and inland, and the size of both the original and
% the reduced dataset, named
% '<file>.landtest.parameters.mat'
%
% (from part 2) Mat-files with logical matrix called 'flag' matching the
% lon/lat matrices of the input file. File names are
% '<file>.landtest_flagged.mat'
%
% (from part 3) All maps in same figure, but saved to file
% '<file>.landtest_flagged.cluster*.png'.
%
% (Always) Indices to files with warnings (currently only about
% overlapping polygons; see below), in the file
% 'landtest_analytics.mat' in the present working directory. 
%
%
% MESSAGES DISPLAYED: 
%
% - Time stamped info for each file written (as well as number of
% clusters, number of positions, length indication of coastline,
% longitude range, latitude range, and an internal code for reason for
% or type of clustering done (c = long coastline; w = widespread
% positions; CK = by the thousands of positions; C2 = two clusters;
% MP/LP= medium/large dataset primitively split by longitude; 1L = one
% location only).  
% - Time stamp and name of every flag file written. 
% - The warning 'two polygons have overlapped' when a point is inside
% two polygons, which means it is likely in a lake or inland sea. More
% severely 'SEVERAL POLYGONS HAVE OVERLAPPED!' when there is a need to
% investigate.
%
%
% EXAMPLE OF USE:
%
% for j=1:FN							% Loop each file.
%   lon=single(ncread(files{j},'LONGITUDE'));			% Load the data positions from file.
%   lat=single(ncread(files{j},'LATITUDE'));			% Load the data positions from file.
%   [ans,olap(j,:)]=landpoints(lon,lat,files{j},[1 1 0],j);	% Test the positions.
%   landpoints(lon,lat,files{j},[0 0 1],j);			% Plot the flagged positions. 
% end
%
%
% REQUIREMENTS:
%
% M_MAP with M_GSHHS with coarse and full coast data installed, and
% the file coastal_and_inland_masks6.mat put alongside this function. 
%
%
% See also CLUSTERDATA INPOLYGON M_MAP M_GSHHS EVENMAT

%% Set your own parameters here: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
temp='~/Downloads/';			% Temporary directory
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

error(nargchk(2,7,nargin));
if nargin<7 | isempty(crit),	crit=[30 50 400];	end	% [maxclust, maxlat, maxcoast] empirically chosen internal criteria.
if nargin<6 | isempty(plotall),	plotall=logical(0);	end	% Logical whether to plot all or only on-land cases.
if nargin<5 | isempty(j),	j=0;			end	% External (file) counter.
if nargin<4 | isempty(part),	part=logical([1 1 0]),	end	% Turn on and off which parts of script to run.
if nargin<3 | isempty(file),	file=[temp,'file'];	end	% Dummy path and filename to name and place the files by.
mes='';								% Short message to add to output messages.
olap=zeros(1,crit(1));						% Notification flag in case of overlapping land patches.
flag=[];							% Empty output for Part 1 and 2.

if (part(1)|part(2))&part(3), 
  part=logical([1 1 0]); 
  warning('LANDPOINTS can either test or plot! Defaulting to testing.');
end

if part(1)|part(2)		% Either these two or one of them or plotting (Part 3)
  
  if part(1) % -------------1) ELIMINATE OPEN OCEAN POINTS AND FLAG INLAND POINTS ---------------
    %if exist('zlon'),inpolygon(lon,lat,zlon,zlat); lon=lon(ans); lat=lat(ans);end	% JUST A TEST FOR SMALLER AREAS OF A FILE
    lom=[min(lon) max(lon)]; lam=[min(lat) max(lat)];					% Ranges of positions in file.
    [oM,oN]=size(lon);									% Store original size.
    LO=lon;LA=lat;
    if oM>1&oN>1, disp('      [ positions are in a matrix!!! ]'); end			% Not supported!
    [ans,IA,IC]=unique([lon(:),lat(:)],'rows','stable'); lon=ans(:,1); lat=ans(:,2);	% Unique positions; ans = A(IA) and A = ans(IC)
    load coastal_and_inland_masks6 LON LAT COASTAL INLAND				% Load the masks
    coi=find(ismember(int16(floor([lon(:),lat(:)]*10)),[LON(COASTAL(:)),LAT(COASTAL(:))],'rows'));	% Indices to unique positions in the coastal zone 
    lai=find(ismember(int16(floor([lon(:),lat(:)]*10)),[LON(INLAND(:)),LAT(INLAND(:))],'rows'));	% Indices to unique positions clearly inland.
    LIC=IA(lai);									% Indices to the original placement of the inland positions
    CIC=IA(coi);									% Indices to the original placement of the positions to be tested
    lon=lon(coi); lat=lat(coi);								% Unique positions in coastal zone.
    % Note: now exactly the same as lon=lon(CIC) on the original loaded lon.
    [M,N]=size(lon);									% Size of the reduced dataset
    system(['rm -f ',file,'.landtest*']);						% Remove old landtest files belonging to this datafile.
    save([file,'.landtest.parameters.mat'],'lom','lam','oM','oN','M','N','CIC','LIC');

    % -------------- EAST-WEST CHECK (possibly missing sign): --------------------- 
    %
    % This could be activated when LIC is a large fraction of oM*oN
    % and in the east. When many points are way inland, there is high
    % possibility that the position cloud fits better in the east.
    % Could make testplots for lon = -lon in special folders for each
    % file to evaluate all positions visually. Possibly the pattern of
    % all position may fit better to a coastline in the west. We have
    % the advantage when testing whole files, that we may see
    % patterns. The landtest could run regardless, and figures checked
    % and reversal of sign done in a separate new run of the landtest.
    % LO=single(ncread(file,'LONGITUDE')); LA=single(ncread(file,'LATITUDE'));	% Load the data positions from file
    % if LIC/(oM*oN)>.01 & LO(LIC)>0
    %   % figure(12); clf
    %   % m_proj('mercator','lon',flipdim(-double(lom),2),'lat',double(lam));
    %   % m_grid; m_gshhs_h('patch',[.8 .8 .8]);
    %   % hli=m_line(-LO,LA,'linestyle','none','marker','.','color','g','markersize',8); % Plot good points
    %   % ht=title([int2str(j),' - ',file,' eastwest_check']); 
    %   % set(ht,'interpreter','none','fontsize',10);set(1,'renderer','painters');
    %   % print('-dpng','-r100',[file,'.landtest.eastwest_check.png']);
    %   disp([datestr(now),' - ',int2str(j),' - ',file,'.landtest.eastwest_check.png plotted!']);
    % end
    
    if ~isempty(coi) % ----- Positions in the coastal zone => CLUSTERING: ---------------
      coslat=max(cosd(lat(:)),1e-5);							% Avoid division by zero at poles
      lom=[min(lon) max(lon)]; lam=[min(lat) max(lat)];					% Ranges of the positions to be tested (coastal zone)
      I=ones(size(lon));		% Cluster number for each observation (init.).
      C=1;				% Set of cluster numbers (init.).
      CN=1;				% Number of clusters (init.).
      if diff(lom)==0 & diff(lam)==0						% It is a time series at same location.
	mes=[mes,'1L']; % one location
	I=1; C=1; CN=1;
      elseif diff(lom) > 15/cosd(min(abs(lam))) | diff(lam) > 15 |  M*N > 1e5	% Locations are wide spread or many.
        % ----- INITIAL CLUSTERING: -----
	if M*N > 0.8*1e5		% There is too much data for clusterdata (1) and inpolygon (2) to handle.
	  mes=[mes,'LP'];		% Split large dataset primitively.
	  C=[1:0.8*1e5:M*N-1 M*N]; CN=length(C);
	  if lom./cosd((lam(2)-lam(1))/2) > lam	% Sort before splitting to give different regions, hence less
	    [ans,is]=sort(lon);			% coastlines, hence faster inpolygon in (2).  
	  else					% Sort by the longest direction.
	    [ans,is]=sort(lat);	
	  end
	  for c=1:CN-1	        
	    I(is(C(c):C(c+1)))=c;		% Assign the first group as first cluster, etc.
	  end
	elseif M*N > 1e4			% There is still too much data for inpolygon (2) to handle. (cluster of >1e4 in 2675 was slow in (2)).
	  mes=[mes,'MP'];			% Split medium dataset primitively.
	  C=[1:1e4:M*N-1 M*N]; CN=length(C);
	  if lom./cosd((lam(2)-lam(1))/2) > lam	% Sort before splitting to give different regions, hence less
	    [ans,is]=sort(lon);			% coastlines, hence faster inpolygon in (2).  
	  else					% Sort by the longest direction.
	    [ans,is]=sort(lat);	
	  end
	  for c=1:CN-1 
	    I(is(C(c):C(c+1)))=c;		% Assign the first group as first cluster, etc.
	  end
	elseif M*N < 1000			% There is less data than normal max for cluster.
	  mes=[mes,'C2'];			% Make 2 clusters.
	  I=clusterdata([lon(:)./coslat,lat(:)],'distance','euclidean','linkage','ward','criterion','distance','maxclust',2);
	else
	  mes=[mes,'CK'];			% Or just cluster by the 1000s.
	  I=clusterdata([lon(:)./coslat,lat(:)],'distance','euclidean','linkage','ward','criterion','distance','cutoff',1000);
	end
	C=unique(I)'; CN=length(C);
	if any(isnan(I)), return; end
	% ----- ADDITIONAL CLUSTERING: -----
	c=1;					% Counter for cluster number.
	while CN<crit(1)			% Maybe split clusters (iterative process).
	  ci=find(I==C(c))';			% The indices to positions belonging in cluster.
	  lom=[min(lon(ci))-.06 max(lon(ci))+.05]; lam=[min(lat(ci))-.03 max(lat(ci))+.04]; % Ranges with carefully added space (do not change).
	  if (diff(lom) > crit(2)/cosd(min(abs(lam))) | diff(lam) > crit(2) ) & length(ci) > 1
	    mes=[mes,'w'];			% Locations are widespread. Make two new clusters:
	    I(ci)=max(C)+clusterdata([lon(ci)./coslat(ci),lat(ci)],'distance','euclidean','linkage','ward','criterion','distance','maxclust',2);
	    C=unique(I)'; CN=length(C); c=1;
	  elseif length(ci) > 1			% Or maybe coastline is too long. 
	    m_proj('albers','lon',lom,'lat',lam);  m_gshhs('c','save',[temp,'temp.mat']); load([temp,'temp'],'ncst');
	    if length(ncst)>crit(3), clear ncst	% Coastline is too long. 
	      mes=[mes,'c']; 
						% Add new cluster number to existing I at the ci points:
	      I(ci)=max(C)+clusterdata([lon(ci)./coslat(ci),lat(ci)],'distance','euclidean','linkage','ward','criterion','distance','maxclust',2);
	      C=unique(I)'; CN=length(C); c=1;
	    else
	      c=c+1;
	    end
	  else
	    c=c+1; 
	  end
	  if c>CN, break; end
	end
	C=unique(I)'; CN=length(C);
      end
      % ----- WRITE COASTFILES AND DISPLAY MESSAGE: -------
      for c=1:CN						% Loop clusters and find their fine coastlines.
	ci=find(I==C(c));					% The indices to positions in cluster (also saved and used later).
	lom=[min(lon(ci))-.06 max(lon(ci))+.05]; lam=[min(lat(ci))-.03 max(lat(ci))+.04]; 
	m_proj('albers','lon',double(lom),'lat',double(lam));
	m_gshhs('f','save',[file,'.landtest.coast.cluster',num2str(c,'%2.2d'),'.mat']); 
                      save([file,'.landtest.coast.cluster',num2str(c,'%2.2d'),'.mat'],'ci','lom','lam','mes','-append');
		      load([file,'.landtest.coast.cluster',num2str(c,'%2.2d'),'.mat'],'ncst'); ncst=single(ncst);
		      save([file,'.landtest.coast.cluster',num2str(c,'%2.2d'),'.mat'],'ncst','-append'); clear ncst
	% These files are needed both for the testing (Part 2) and the plotting (Part 3).
	m_gshhs('c','save',[temp,'temp.mat']); load([temp,'temp'],'ncst'); LC=length(ncst); clear ncst
	[CN,length(ci),LC,ceil(diff(lom)),ceil(diff(lam))];
	disp([datestr(now),' - ',int2str(j),' - ',file, ...
	      '.landtest.coast.cluster',num2str(c,'%2.2d'),'.mat (',num2str(ans,'%2u\t%6u\t%6u\t%3u\t%2u'),') ',mes]);
      end
    end
  end
  
  
  
  
  if part(2) % -------------- 2) TEST AND FLAG POSITIONS -------------------------------
    load([file,'.landtest.parameters.mat']);
    flag=logical(zeros(oM,oN));						% Flag matrix for the original positions (init.).
    flagg=zeros(M,N);							% Flag matrix for the unique coastal zone positions (init.).
    if ~isempty(CIC)
      if ~part(1)
	lon=lon(CIC); lat=lat(CIC);		% Reduce to the unique coastal zone positions if not already reduced in Part 1.
      end
      clusterfiles=dir([file,'.landtest.coast.cluster*']);		% Read all the corresponding coastfiles.
      for c=1:length(clusterfiles)	% Test all positions cluster by cluster, flagging the same flag-matrix.
	flg=zeros([M,N]);							% Full size flag matrix, but internal to this loop.
	load([clusterfiles(c).folder,filesep,clusterfiles(c).name]);	% Load coastlines and indices.
	for i=1:length(k)-1						% Flag the on-land data positions in this cluster.
	  flg(ci)=flg(ci)+inpolygon(lon(ci),lat(ci),ncst(k(i)+1:k(i+1)-1,1),ncst(k(i)+1:k(i+1)-1,2));
	end
	if any(flg>1,'all')		% Test if flg>1 anywhere, because then polygons have overlapped.
	  disp([datestr(now),' - ',int2str(j),' - ',file,'-cluster',int2str(c),' – Two polygons have overlapped.']);
	  flg(flg==2)=0;							% Invert the logical (it's a lake)
	  olap(1,c)=1;
	  if any(flg>2,'all')		% Safety message in case of several overlaps. Likely island in lake.
	    disp([datestr(now),' - ',int2str(j),' - ',file,'-cluster',int2str(c),' – SEVERAL POLYGONS HAVE OVERLAPPED!']);
	    olap(1,c)=2;
	  end 
	else	
	  olap(1,c)=0;
	end
	if strcmp(mes,'1L') & length(ci)==1
	  flg=repmat(flg(ci),[M,N]);		% expand to all points if tested for one position
	elseif strcmp(mes,'1L') & length(ci)~=1				
	  warning('1L but not single ci!');	% Warning (and stop) if something wrong with indexing single locations.
	end
	~flagg; flagg(ans)=flg(ans);		% Add the cluster's flg to the file's flagg, but avoid double
						% flagging between the rectangular cluster areas.
      end
      flagg=logical(flagg);			% Make logical and keep these flags to the unique coastal zone positions (for plotting).
      
      % -------------- CLOSENESS TO COASTLINE: ----------------------------- 
      % Unflag the flagged points close to the coastline. Criterion will be closer than half the resolution of the nearest coastline points.
      ii=find(flagg);
      for i=1:length(ii)						% Loop the flagged points in the coastal zone of this file.
	A=[lon(ii(i)),lat(ii(i))];					% The position in question
        % 1) Quickly find the nearest coastal points: There are ncst distances up to 5 km ≈ 0.05 deg lat, but we are using the
	% criteria of half a coastpoint distance, so we use .025 ≈ 2.8 km here. Furthermore longitudes are normalised on the
	% latitudes, to search in as equidistant a space as possible.
	jB=find(A(1)-.025/cosd(A(2))<ncst(:,1) & ncst(:,1)<A(1)+.025/cosd(A(2)) & A(2)-.025<ncst(:,2) & ncst(:,2)<A(2)+.025); 
	if length(jB)>=3
	  B=ncst(jB,:);	nB=size(B,1);						% The coastpoints to use initially
          % 2) Now test each found coastpoint using real distances:
	  reshape([repmat(shiftdim(A,-1),1,nB,1);shiftdim(B,-1)],nB*2,2);	% Order into A,B,A,B,...
	  sw_dist(ans(:,2),ans(:,1));						% Distances between A,B,A,B,... 
	  [closest,j1]=sort(ans(1:2:end));					% Sort the A to B differences only
      	  % 3) Find the closest coast point's closed polygon: 
	  find(k>=jB(j1(1)));
	  sta=k(ans(1)-1)+1; sto=k(ans(1))-2;					% The closed polygon (not closed)
	  if jB(j1(1))==k(ans(1)-1)+1,		cp=[sto sta+[0 1]];		% Facilitate for endpoints: first,
	  elseif jB(j1(1))==k(ans(1))-2,	cp=[sto+[-1 0] 1];		% last,
	  else					cp=jB(j1(1))+[-1:1];		% between.	     
	  end									
	  B=ncst(cp,:);								% The coastpoints to test against
	  coastres=sw_dist(B(:,2),B(:,1));					% Distances between the 3 nearest coast points 
	  if closest(1)<max(coastres)/2						% Evaluate closeness
	    flagg(ii(i))=0;							% If too close, unflag the position.   
      	    % % Make control figures for this routine:
	    % figure(11);clf;
	    % m_proj('Azimuthal Equidistant','lon',double(A(1)),'lat',double(A(2)),'rad',double(A)+[.025 -.025]); 
	    % m_grid; m_gshhs_f('patch',[.7 .7 .7]);
	    % hgA=m_line(A(1),A(2),'linestyle','none','marker','.','color','g','markersize',8);
	    % hgC=m_line(B(:,1),B(:,2),'linestyle','-','marker','.','color','r','markersize',8);
	    % print('-dpng','-r100',[file,'.reversed_landpoint',num2str(i,'%3.3d'),'.png']);  
	  end
	end
      end 
    end
    flag(CIC)=flagg;
    flag(LIC)=logical(1);			% Now, flag the clearly inland positions.	    
    % ----- WRITE FLAG FILES AND DISPLAY MESSAGE: -------
    system(['rm -f ',file,'.landtest_flagged*']);	% Remove old flag files belonging to this datafile.
    if any(flag,'all')					% Save both these logical flag matrices, etc.:
      save([file,'.landtest_flagged.mat'],'flagg','flag');
      disp([datestr(now),' - ',int2str(j),' - ',file,'.landtest_flagged.mat']);
    end
  end
  
  
  
  
elseif part(3) % -------------- 3) MAKE PLOTS FOR CONTROL -----------------------------------
  screensize=get(0,'Screensize');
  load([file,'.landtest.parameters.mat'],'oM','oN');
  noflags=logical(0);						% (used for plotall)
  try
    load([file,'.landtest_flagged.mat']);			% Load flags ...
  catch	
    noflags=logical(1);
    if ~plotall, return; end					% ... but return if no flag/file.
  end
  system(['rm -f ',file,'.landtest_flagged*.png']);		% Remove old figure files belonging to this datafile.
  if noflags, flag=logical(zeros(oM,oN)); end
  %if noflags, flag=logical(zeros(size(lon))); end
  clusterfiles=dir([file,'.landtest.coast.cluster*']);		% Read all the corresponding coastfiles.
  for c=1:length(clusterfiles)					% Loop cluster by cluster.
    load([clusterfiles(c).folder,filesep,clusterfiles(c).name],'lom','lam','ci');	% Load the cluster's region and indices.
    i=find(lom(1)<=lon & lon<lom(2) & lam(1)<lat & lat<lam(2));				% All positions in that cluster's region.
    if any(flag(i),'all') | plotall | olap(1,c)			% Only plot if flagged data in region, overlap, or plotall.  
      figure(1); set(gcf,'OuterPosition',screensize); clf;						% Make map 
      m_proj('mercator','lon',double(lom),'lat',double(lam));
      m_grid;
      m_usercoast([clusterfiles(c).folder,filesep,clusterfiles(c).name],'patch',[.7 .7 .7]);
      hg=m_line(lon(i(~flag(i))),lat(i(~flag(i))),'linestyle','none','marker','.','color','g','markersize',8); % Plot good points
      hb=m_line(lon(i( flag(i))),lat(i( flag(i))),'linestyle','none','marker','.','color','r','markersize',8); % Plot bad points
      ht=title([int2str(j),' - ',file,' - Cluster ',int2str(c)]); 
      set(ht,'interpreter','none','fontsize',10);
      % ----- PRINT FIGURES AND DISPLAY MESSAGE: -------
      set(1,'renderer','painters');
      if olap(1,c), nadd='_overlap'; else nadd=''; end		% Mark figures with ovelapping polygons
      print('-dpng','-r100',[file,'.landtest_flagged.cluster',num2str(c,'%2.2d'),nadd,'.png']);  
      [length(ci),ceil(diff(lom)),ceil(diff(lam))];
      disp([datestr(now),' - ',int2str(j),' - ',file,'.landtest_flagged.cluster',num2str(c,'%2.2d'),nadd,'.png (',num2str(ans,'%6u\t%3u\t%2u'),')']);
    end
  end
end
  