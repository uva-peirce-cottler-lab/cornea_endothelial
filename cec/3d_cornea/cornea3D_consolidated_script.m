%% Getting matrix data from .lsm file
if true
	% Cleaning the workspace
	clearvars

	% Getting the file and path to the .lsm file
	[filename,pathname] = uigetfile('*.lsm') ;

	% Reading the lsm data into a matrix
	imageData = lsmread([pathname,filename]) ;

	% Rearranging the dimensions
	imageData = squeeze(imageData) ;
	imageData = permute(imageData,[2,3,1]) ; % 2 and 3 can be swapped here but it's also symmetric

	% Saving the data
	save([filename(1:end-4),'-imageData.mat'],'imageData')
end

%% Create binary data
if true
	delay = 0.01 ;
	threshold = 20 ;
	SE = strel('disk',3,0) ;
	
	binaryData = logical(imageData) ;

%%%%% Start: Constants
	if true
		averagingFilterSize = 64 ;
		grey2BWthreshold = 0.01 ;
		minCCArea = 50 ;

		ConwaysIterations	 = 1 ;
		bridgeFillIterations = 1 ;

		bridgeIterations   = 1 ;
		fillIterations     = 1 ;
		majorityIterations = 1;
		thinIterations     = inf ;

		survive = [4,5,6,7,8] ;
		born    = [4,5,6,7,8] ;

		row = 1024 ;
		col = 1024 ;

		pR = [1 1:row-1];
		qR = [2:row row];
		pC = [1 1:col-1];
		qC = [2:col col];
	end
%%%%% End: Constants

%%%%% Start: Analysis
	grey = double(imageData(:,:, 1 )) ;

	greynbrs = ...
		grey(:,pC)  + grey(:,qC) + ...
		grey(pR,:)  + grey(qR,:) + ...
		grey(pR,pC) + grey(qR,qC) + ...
		grey(pR,qC) + grey(qR,pC);

	greynbrs = greynbrs / (255*8) ;

	meanValue = imfilter( greynbrs , ...
		fspecial('average', [averagingFilterSize , ...
		averagingFilterSize]) , 'replicate' ) ;

	BW_1 = ( greynbrs - meanValue ) > grey2BWthreshold ;
	BW_1 = bwmorph( BW_1 , 'majority' , 4 ) ;

	BW_1 = imerode(BW_1,SE) ;

	binaryData(:,:,1) = BW_1 ;

	for i = (2:(size(imageData,3)))

		grey = double(imageData(:,:, i )) ;

		greynbrs = ...
			grey(:,pC)  + grey(:,qC) + ...
			grey(pR,:)  + grey(qR,:) + ...
			grey(pR,pC) + grey(qR,qC) + ...
			grey(pR,qC) + grey(qR,pC);

		greynbrs = greynbrs / (255*8) ;

		meanValue = imfilter( greynbrs , ...
			fspecial('average', [averagingFilterSize , ...
			averagingFilterSize]) , 'replicate' ) ;

		BW_1 = ( greynbrs - meanValue ) > grey2BWthreshold ;
		BW_1 = bwmorph( BW_1 , 'majority' , 4 ) ;
		
		BW_1 = imerode(BW_1,SE) ;
		
		binaryData(:,:,i) = BW_1 ;

		disp(i/size(imageData,3))

	end
	
%%%%% End: Analysis

	% Saving the binary data
	save([filename(1:end-4),'-binaryData.mat'],'binaryData')

	disp('Done creating binary data!')
end

%% Connected Component filtering
if true

%%%%% Start: Filter out CC with Volume < 100 and identify by layer

	%%%%%%%%%% THIS IS SET BY THE USER 
	LayerCutoffs = [21,299] ;
	%%%%%%%%%% THIS IS SET BY THE USER
		
	z=zeros(1024,'logical');
	x=zeros(1024,1,404,'logical');
	y=zeros(1,1026,404,'logical');
	binaryData = cat( 3 , z,binaryData,z ) ;
	binaryData = cat( 2 , x,binaryData,x ) ;
	binaryData = cat( 1 , y,binaryData,y ) ;
		
	CC = bwconncomp(binaryData,6) ;
	stats = regionprops3(CC,'Centroid','Volume','VoxelIdxList') ;
	stats( stats.Volume<100 , : ) = [] ;

	Layer = ones(length(stats.Volume),1) ;
	Layer = Layer + (stats.Centroid(:,3)>LayerCutoffs(1)) ;
	Layer = Layer + (stats.Centroid(:,3)>LayerCutoffs(2)) ;

	stats = [stats,table(Layer)] ;
%%%%% End: Filter out CC with Volume < 100 and identify by layer	

	
%%%%% Start: Create Seperate Layers
		botLayer = logical(binaryData.*false) ;
		midLayer = botLayer ;
		topLayer = botLayer ;

		for i = 1:length(stats.Layer)
			if stats.Layer(i) == 1
				botLayer( stats.VoxelIdxList{i} ) = true ;
			end
		end

		for i = 1:length(stats.Layer)
			if stats.Layer(i) == 2
				midLayer( stats.VoxelIdxList{i} ) = true ;
			end
		end

		for i = 1:length(stats.Layer)
			if stats.Layer(i) == 3
				topLayer( stats.VoxelIdxList{i} ) = true ;
			end
		end
%%%%% End: Create Seperate Layers
	
	% Saving the layer data
	save([filename(1:end-4),'-layerData.mat'],'botLayer','midLayer','topLayer')
	
	disp('Done creating layer data!')
end

%% Create isosurface New
if true
	% Include boundary slices in the lower layers
	botLayer(:,:,1:27) = botLayer(:,:,1:27) + midLayer(:,:,1:27) ;
	midLayer(:,:,1:25) = 0 ;
	topLayer(:,:,299:end) = topLayer(:,:,299:end) + midLayer(:,:,299:end) ;
	midLayer(:,:,301:end) = 0 ;
	
	startBounds = [1,129, 1,129, 1,52] ;
	b = [ 1,1 , 1,1 , 1,1 ] ;
	kBounds = [1,31;
		29,52;
		50,152;
		150,202;
		200,252;
		250,301;
		299,352;
		350,404] ;
	
	jBounds = [1,128;127,257;255,385;383,512;511,641;639,769;767,897;895,1026];
	iteration = 1 ;

	for k = 1:8 % Z sets in groups of 50 layers
		b(5:6) = kBounds(k,:) ;

		for j = 1:8 % X sets in groups of 128 pixels
			b(3:4) = jBounds(j,:) ;

			for i = 1:8 % Y sets in groups of 128 pixels

				% b stands for bounds
% 				b(1:4) = [ [((i-1)*128)-1,(i*128)+1] , [((j-1)*128)-1,(j*128)+1] ] ;
% 				b( b<1 ) = 1 ;
% 				b( b>1024 ) = 1024 ;
				b(1:2) = jBounds(i,:) ;

				smoothedB = smooth3(smooth3(smooth3(smooth3(botLayer(b(1):b(2),b(3):b(4),b(5):b(6)))))) ;
				smoothedM = smooth3(smooth3(smooth3(smooth3(midLayer(b(1):b(2),b(3):b(4),b(5):b(6)))))) ;
				smoothedT = smooth3(smooth3(smooth3(smooth3(topLayer(b(1):b(2),b(3):b(4),b(5):b(6)))))) ;
				
				
				fvB( iteration ) = isosurface(smoothedB,0.5) ; %#ok<*SAGROW>
				fvM( iteration ) = isosurface(smoothedM,0.5) ; 
				fvT( iteration ) = isosurface(smoothedT,0.5) ; 
				
				if ~isempty(fvB(iteration).vertices)
					fvB(iteration).vertices(:,1) = fvB(iteration).vertices(:,1) + (b(3)-startBounds(3)) ;
					fvB(iteration).vertices(:,2) = fvB(iteration).vertices(:,2) + (b(1)-startBounds(1)) ;
					fvB(iteration).vertices(:,3) = fvB(iteration).vertices(:,3) + (b(5)-startBounds(5)) ;
				end
				
				if ~isempty(fvM(iteration).vertices)
					fvM(iteration).vertices(:,1) = fvM(iteration).vertices(:,1) + (b(3)-startBounds(3)) ;
					fvM(iteration).vertices(:,2) = fvM(iteration).vertices(:,2) + (b(1)-startBounds(1)) ;
					fvM(iteration).vertices(:,3) = fvM(iteration).vertices(:,3) + (b(5)-startBounds(5)) ;
				end
				
				if ~isempty(fvT(iteration).vertices)
					fvT(iteration).vertices(:,1) = fvT(iteration).vertices(:,1) + (b(3)-startBounds(3)) ;
					fvT(iteration).vertices(:,2) = fvT(iteration).vertices(:,2) + (b(1)-startBounds(1)) ;
					fvT(iteration).vertices(:,3) = fvT(iteration).vertices(:,3) + (b(5)-startBounds(5)) ;
				end

				% Shift the first column of vertices for a shift in y
				% Shift the second column of vertices for a shift in x

				iteration = iteration + 1 ;

			end

		end
		disp(k/8)
	end
	
	% Saving the display data
	save([filename(1:end-4),'-isosurfaceData.mat'],'fvB','fvM','fvT')
	
	disp('Done creating isosurface data!')

end






