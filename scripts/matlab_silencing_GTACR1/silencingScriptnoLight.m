%% Script to make struct array for ACR1 tracked data
% This only works with speicific folder,subfolder arrangement as in
% selected_60-30_powdered_ACR1
subFolders = dir(pwd);
realSubFolders = subFolders([subFolders.isdir]');
realSubFolders = realSubFolders(3:end);

% Initialise video temporal parameters
fps = 30; %frames per second

maxLastFrameNumber = 160*fps;

%Initialize smoothing parameters for velocity attributes
smoothWindow = 15; %frames
angSmoothWindow= 15; %frames
walkThreshold = 1.5; %mm/s
%%
data={};

%Collect per genotype data in the following for loop
for i = 1:length(realSubFolders)
    [currGenotype, ACR] = strtok(realSubFolders(i).name,'_');
    currPath = [pwd '\' realSubFolders(i).name];
    
    data(i).Genotype = currGenotype;
    data(i).Path = currPath;
    
    genotypeSubFolders = dir(currPath);
    genotypeSubFolders = genotypeSubFolders([genotypeSubFolders.isdir]');
    genotypeSubFolders = genotypeSubFolders (3:end);
    data(i).numFlies = length(genotypeSubFolders);
    
    featureFiles = getAllExtFiles(currPath,'feat.mat');
    trackFiles = getAllExtFiles(currPath,'track.mat');
    CalFiles = getAllExtFiles(currPath,'calibration.mat');
    
    
    %Collect per fly data in the following for loop
    for j=1:data(i).numFlies
        currFeatureFile = char(featureFiles(j));
        tempMatFileFeat = matfile(currFeatureFile);
        tempMatFileTrk = matfile(char(trackFiles(j)));
        tempMatFileCal = matfile(char(CalFiles(j)));
        tempCal = tempMatFileCal.calib;
        tempFeature = tempMatFileFeat.feat;
        tempTrk = tempMatFileTrk.trk;
        tempVel = tempFeature.data(:,:,strcmp(tempFeature.names,'vel'));
        if(isempty(tempVel))
            break
        end
        tempAngVel = tempFeature.data(:,:,strcmp(tempFeature.names,'ang_vel'));
        tempOri = tempTrk.data(:,:,strcmp(tempTrk.names,'ori'));
        tempPosX = tempTrk.data(:,:,strcmp(tempTrk.names,'pos x'));
        tempPosY = tempTrk.data(:,:,strcmp(tempTrk.names,'pos y'));
        
        %at times velocity arrays are 1 or 2 frames longer/shorter than
        %expected. In such cases the following condition makes vel and
        %andVel arrays into consistent size as manually provided in variable
        %maxLastFrameNumber. 
        
        if length(tempVel)<maxLastFrameNumber
                   fff=2222 
           currFeatLastFrame = length(tempVel);
           tempVel(end+1:maxLastFrameNumber)=mean(tempVel);%(currFeatLastFrame);
           tempAngVel(end+1:maxLastFrameNumber)=mean(tempAngVel);%(currFeatLastFrame);
        else
            fff=111
           tempVel = tempVel(1:maxLastFrameNumber);
           tempAngVel = tempAngVel(1:maxLastFrameNumber);
        end
        
    % Since Trk and Feat files created seperately do if statement once more
        if length(tempOri)<maxLastFrameNumber
           currTrkLastFrame = length(tempOri);
           tempOri(end+1:maxLastFrameNumber)=tempOri(currTrkLastFrame);
           tempPosX(end+1:maxLastFrameNumber)=tempPosX(currTrkLastFrame);
           tempPosY(end+1:maxLastFrameNumber)=tempPosY(currTrkLastFrame);
        else
           tempOri = tempOri(1:maxLastFrameNumber);
           tempPosX = tempPosX(1:maxLastFrameNumber);
           tempPosY = tempPosY(1:maxLastFrameNumber);
        end
        
        timetags = ([1:length(tempVel)])/fps;
        
        
        data(i).perFlyData(j).feat = tempFeature;
        data(i).perFlyData(j).trk = tempTrk;
        data(i).perFlyData(j).ppm = tempCal.PPM;
        
        %The following if statement converts manually annotated backward walking time periods as
        %negative velocity time points
        if exist(strcat(currFeatureFile(1:length(currFeatureFile)-9),'-actions.mat'),'file');
                
           tempMatFileActions = matfile(strcat(currFeatureFile(1:length(currFeatureFile)-9),'-actions.mat'));
           tempbehs = tempMatFileActions.behs;
           tempbouts = tempMatFileActions.bouts;
           backBouts = cell2mat(tempbouts(1,1));
           backBouts = backBouts(:,[1 2],:);
           backBouts = backBouts/fps;
           tempCondArr = zeros(length(backBouts),length(timetags));
           for ctrtemp=(1:length(backBouts))
               tempCondArr(ctrtemp,:)=transpose(timetags>backBouts(ctrtemp,1) &timetags<backBouts(ctrtemp,2));
           end
           condBack = sum(tempCondArr,1);
           tempVel(condBack>0) = tempVel(condBack>0)*-1;
        end
        
        
        theta2 = tempOri(2:end);
        theta1 = tempOri(1:end-1);
        smooth_kernel = [1 2 1]/4;
        ori_diffRAW = mod(theta1+pi/2 - theta2,pi)-pi/2;
        tempAngVelSigned = [ori_diffRAW(1) ori_diffRAW]*fps;
        tempAngVelSigned(2:end-1) = conv(tempAngVelSigned,smooth_kernel,'valid');
        
            
        data(i).perFlyData(j).vel = tempVel;
        data(i).perFlyData(j).angVel = tempAngVel;
        data(i).perFlyData(j).angVelSigned = tempAngVelSigned;
        data(i).perFlyData(j).ori = tempOri;
        data(i).perFlyData(j).PosX = tempPosX;
        data(i).perFlyData(j).PosY = tempPosY;
        
        data(i).perFlyData(j).smoothVel = (smooth(tempVel,smoothWindow))';
        
        walked = data(i).perFlyData(j).smoothVel>walkThreshold;
        
        
        data(i).perFlyData(j).meanWalkedVel = mean(data(i).perFlyData(j).smoothVel(walked));
        
        
        data(i).perFlyData(j).smoothAngVel = (smooth(tempAngVel,angSmoothWindow))';
       
        data(i).perFlyData(j).smoothAngVelSigned = (smooth(tempAngVelSigned,angSmoothWindow))';
        
       
        
    end
    data(i).catVel = cat(1,data(i).perFlyData.vel);
    data(i).catSmoothVel = cat(1,data(i).perFlyData.smoothVel);
    
    data(i).catAngVel = cat(1,data(i).perFlyData.angVel);
    data(i).catSmoothAngVel = cat(1,data(i).perFlyData.smoothAngVel);
    
    data(i).catAngVelSigned = cat(1,data(i).perFlyData.angVelSigned);
    data(i).catSmoothAngVelSigned = cat(1,data(i).perFlyData.smoothAngVelSigned);
   
   
    
end

%% Create full data heat map
fullSmoothVelData = cat(1,data.catSmoothVel);
labels = {data.Genotype};
minVel = 0;
maxVel = 20;
n1 = cumsum([data.numFlies])';
n2 = ([data.numFlies]/2)';
n=n1-n2;
n1StringArr = cell(1,length(n1));
n1StringArr(:) = {'----'};
newN=vertcat(n,n1);
[newN,index]=sort(newN);
newLabels = horzcat(labels,n1StringArr);
newLabels = newLabels(index);

clims=[minVel,maxVel];
figure;
imagesc(fullSmoothVelData,clims);title('CsChrimson Activation');
set(gca,'YTick',newN','YTickLabel',newLabels,'XTick',1:15*fps:length(fullSmoothVelData),'XTickLabel',0:15:length(fullSmoothVelData)/fps);
%set(gca,'YTick',n1','YTickLabel',n1StringArr);
ylabel('Genotypes');xlabel('time');colorbar;
%savefig(['modeError - ' num2str(ZPLANE) '.fig']);
%% walk Bout analysis
walkThreshold = 1.5;%mm/s
framesArray = 1:maxLastFrameNumber;
timeArray = 0:1/fps:(maxLastFrameNumber-1)/fps;
%timeLimitMax = 5;%min
%timeLimitMin = 0; %min
frameLimitMax = maxLastFrameNumber;
frameLimitMin = 1;
frameLimitWindow = framesArray>frameLimitMin & framesArray<frameLimitMax;
boutDurationLimit = 15;%frames
jumpLimit = 30; %mm/s


for i = 1:length(data)
    data(i).realFlyNumber = size(data(i).catSmoothVel,1);
    data(i).walked = data(i).catSmoothVel>walkThreshold;% & data(i).catSmoothVel<jumpLimit;
    data(i).walkBegin = data(i).walked & ~cat(2,false(data(i).realFlyNumber,1),data(i).walked(:,1:end-1));
    data(i).walkEnd = data(i).walked & ~cat(2,data(i).walked(:,2:end),false(data(i).realFlyNumber,1));
    
    for j = 1:data(i).realFlyNumber
             
    
        data(i).perRealFlyData(j).beginFrame = find(data(i).walkBegin(j,:));
        data(i).perRealFlyData(j).endFrame = find(data(i).walkEnd(j,:));
        %data(i).perRealFlyData(j).beginFrame = data(i).perRealFlyData(j).beginFrameAll(frameLimitWindow);
        %data(i).perRealFlyData(j).beginFrame = data(i).perRealFlyData(j).beginFrameAll(frameLimitWindow);
        
        if length(data(i).perRealFlyData(j).beginFrame) > length(data(i).perRealFlyData(j).endFrame)
            data(i).perRealFlyData(j).endFrame = cat(2,data(i).perRealFlyData(j).endFrame,frameLimitMax);
        else
            if length(data(i).perRealFlyData(j).beginFrame) < length(data(i).perRealFlyData(j).endFrame)
                data(i).perRealFlyData(j).beginFrame = cat(2,frameLimitMin,data(i).perRealFlyData(j).beginFrame);
            end
        end
        data(i).perRealFlyData(j).numBouts = length(data(i).perRealFlyData(j).beginFrame);        
        data(i).perRealFlyData(j).boutDuration=data(i).perRealFlyData(j).endFrame-data(i).perRealFlyData(j).beginFrame;
        boutFilter = data(i).perRealFlyData(j).boutDuration>boutDurationLimit;
        data(i).perRealFlyData(j).beginFrame = data(i).perRealFlyData(j).beginFrame(boutFilter);
        data(i).perRealFlyData(j).endFrame = data(i).perRealFlyData(j).endFrame(boutFilter);
        data(i).perRealFlyData(j).numBouts = length(data(i).perRealFlyData(j).beginFrame);
        data(i).perRealFlyData(j).boutDuration = data(i).perRealFlyData(j).boutDuration(boutFilter);
        
        data(i).perRealFlyData(j).boutMedianVel = zeros(1,data(i).perRealFlyData(j).numBouts,1);
        data(i).perRealFlyData(j).boutMedianAngVel = zeros(1,data(i).perRealFlyData(j).numBouts,1);
        for k = 1:data(i).perRealFlyData(j).numBouts
            data(i).perRealFlyData(j).boutMedianVel(k) = mean(data(i).catSmoothVel(j,data(i).perRealFlyData(j).beginFrame(k):data(i).perRealFlyData(j).endFrame(k)));
            data(i).perRealFlyData(j).boutMedianAngVel(k) = mean(data(i).catSmoothAngVel(j,data(i).perRealFlyData(j).beginFrame(k):data(i).perRealFlyData(j).endFrame(k)));
        end
        data(i).perRealFlyData(j).meanBoutMedianVel = mean(data(i).perRealFlyData(j).boutMedianVel,'omitnan');
        data(i).perRealFlyData(j).meanBoutMedianAngVel = mean(data(i).perRealFlyData(j).boutMedianAngVel,'omitnan');
        data(i).perRealFlyData(j).meanBoutDuration = mean(data(i).perRealFlyData(j).boutDuration,'omitnan');
        
    end
    data(i).catBoutVel = cat(2,data(i).perRealFlyData.boutMedianVel);
    data(i).catBoutAngVel = cat(2,data(i).perRealFlyData.boutMedianAngVel);
    data(i).catBoutDuration = cat(2,data(i).perRealFlyData.boutDuration);
    data(i).meanCatBoutVel = cat(2,data(i).perRealFlyData.meanBoutMedianVel);
    data(i).meanCatBoutAngVel = cat(2,data(i).perRealFlyData.meanBoutMedianAngVel);
    data(i).meanCatBoutDuration = cat(2,data(i).perRealFlyData.meanBoutDuration);
    data(i).catNumBouts = cat(2,data(i).perRealFlyData.numBouts);

end
%%
figure;
hold on;
scatter3(data(1).catBoutVel,data(1).catBoutAngVel,data(1).catBoutDuration,'MarkerEdgeColor','b');
scatter3(data(2).catBoutVel,data(2).catBoutAngVel,data(2).catBoutDuration,'MarkerEdgeColor','r');
hold off;
axis([0,20,0,8,0,2500]);
%%
figure;
hold on;
scatter3(data(1).meanCatBoutVel,data(1).meanCatBoutAngVel,data(1).meanCatBoutDuration,'MarkerEdgeColor','b');
scatter3(data(2).meanCatBoutVel,data(2).meanCatBoutAngVel,data(2).meanCatBoutDuration,'MarkerEdgeColor','r');
hold off;
%axis([-2,5,-2,5,-2,5]);

%% merge and cluster
superCatBoutVel = zscore(cat(2,data.catBoutVel));
superCatBoutAngVel = zscore(cat(2,data.catBoutAngVel));
superCatBoutDuration = zscore(cat(2,data.catBoutDuration));
f1Bouts = length(data(1).catBoutVel);
f2Bouts = length(data(2).catBoutVel);
figure;
hold on;
scatter3(superCatBoutVel(1:f1Bouts),superCatBoutAngVel(1:f1Bouts),superCatBoutDuration(1:f1Bouts),'MarkerEdgeColor','b');
scatter3(superCatBoutVel(f1Bouts+1:end),superCatBoutAngVel(f1Bouts+1:end),superCatBoutDuration(f1Bouts+1:end),'MarkerEdgeColor','r');
hold off;
axis([-2,5,-2,5,-2,12]);
figure;
hold on;
scatter(superCatBoutVel(f1Bouts+1:end),superCatBoutDuration(f1Bouts+1:end),'MarkerEdgeColor','r');
scatter(superCatBoutVel(1:f1Bouts),superCatBoutDuration(1:f1Bouts),'MarkerEdgeColor','b');

line([2,2],ylim,'LineStyle','--','Color','k');
line(xlim,[2,2],'LineStyle','--','Color','k');
hold off;
axis([-2,5,-2,12]);

z_threshold = 2;
zNegativeThreshold = 2;
boutClass1 = superCatBoutVel >= z_threshold & superCatBoutDuration >= z_threshold;
boutClass2 = superCatBoutVel >= z_threshold & superCatBoutDuration < zNegativeThreshold;
boutClass3 = superCatBoutVel < zNegativeThreshold & superCatBoutDuration >= z_threshold;
boutClass4 = superCatBoutVel < z_threshold & superCatBoutDuration < z_threshold;

flyIdentity = [repmat(1,1,f1Bouts),repmat(2,1,f2Bouts)]; 
% fly1Class1 = (nnz(flyIdentity(boutClass1) == 1)/f1Bouts)*100;
% fly1Class2 = (nnz(flyIdentity(boutClass2) == 1)/f1Bouts)*100;
% fly1Class3 = (nnz(flyIdentity(boutClass3) == 1)/f1Bouts)*100;
% fly1Class4 = (nnz(flyIdentity(boutClass4) == 1)/f1Bouts)*100;
% fly2Class1 = (nnz(flyIdentity(boutClass1) == 2)/f2Bouts)*100;
% fly2Class2 = (nnz(flyIdentity(boutClass2) == 2)/f2Bouts)*100;
% fly2Class3 = (nnz(flyIdentity(boutClass3) == 2)/f2Bouts)*100;
% fly2Class4 = (nnz(flyIdentity(boutClass4) == 2)/f2Bouts)*100;
fly1Class1 = (nnz(flyIdentity(boutClass1) == 1));
fly1Class2 = (nnz(flyIdentity(boutClass2) == 1));
fly1Class3 = (nnz(flyIdentity(boutClass3) == 1));
fly1Class4 = (nnz(flyIdentity(boutClass4) == 1));
fly2Class1 = (nnz(flyIdentity(boutClass1) == 2));
fly2Class2 = (nnz(flyIdentity(boutClass2) == 2));
fly2Class3 = (nnz(flyIdentity(boutClass3) == 2));
fly2Class4 = (nnz(flyIdentity(boutClass4) == 2));
figure;
bar([[fly1Class1;fly2Class1],[fly1Class2;fly2Class2],[fly1Class3;fly2Class3],[fly1Class4;fly2Class4]]');

%% k means
X = [superCatBoutVel',superCatBoutAngVel',superCatBoutDuration'];
opts = statset('Display','final');
[idx,C] = kmeans(X,4,'Distance','sqeuclidean',...
    'Replicates',1000,'Options',opts);


figure;
scatter(X(idx==1,1),X(idx==1,2),'b')
hold on
scatter(X(idx==2,1),X(idx==2,2),'g')
hold on;
scatter(X(idx==3,1),X(idx==3,2),'r')
hold on;
scatter(X(idx==4,1),X(idx==4,2),'m')
xlabel('vel');ylabel('duration');
hold off;
%% Create total duration metrics for no light condition, i.e. no trial averaging
fps=30;
walkThreshold = 1.5;
totalduration = length(data(1).catSmoothVel)/fps;
for i = 1:length(data)
    data(i).Fulldist = mean(data(i).catSmoothVel,2)*totalduration/fps;
    data(i).meanAngVel = mean(data(i).catSmoothAngVel,2);
    data(i).Walked = data(i).catSmoothVel>walkThreshold;
    data(i).WalkDurationFraction = sum(data(i).Walked,2)/length(data(i).Walked);
    for j = 1:data(i).numFlies
        data(i).perFlyData(j).meanWalkedVel = mean(data(i).perFlyData(j).smoothVel(data(i).perFlyData(j).smoothVel>walkThreshold));
        data(i).perFlyData(j).meanWalkedAngVel = mean(data(i).perFlyData(j).smoothAngVel(data(i).perFlyData(j).smoothVel>walkThreshold));
    end
    data(i).catMeanWalkedVel = cat(1,data(i).perFlyData.meanWalkedVel);
    data(i).catMeanWalkedAngVel = cat(1,data(i).perFlyData.meanWalkedAngVel);
end
%% Create Distance file for plotting in Graphpad
totalDistAvgONData=cat(2,{data.catTotalDistAvgON});
maxLength = max(cellfun(@length,totalDistAvgONData));
totalDistAvgONMat= cell2mat(cellfun(@(x)cat(1,x,nan(maxLength-length(x),1)),totalDistAvgONData,'UniformOutput',false));


totalDistAvgOFFData=cat(2,{data.catTotalDistAvgOFF});
maxLength = max(cellfun(@length,totalDistAvgOFFData));
totalDistAvgOFFMat= cell2mat(cellfun(@(x)cat(1,x,nan(maxLength-length(x),1)),totalDistAvgOFFData,'UniformOutput',false));

%% Create mean Walked Vel file for plotting in Graphpad
meanWalkedVelAvgONData=cat(2,{data.catMeanWalkedVelAvgON});
maxLength = max(cellfun(@length,meanWalkedVelAvgONData));
meanWalkedVelAvgONMat= cell2mat(cellfun(@(x)cat(1,x,nan(maxLength-length(x),1)),meanWalkedVelAvgONData,'UniformOutput',false));


meanWalkedVelAvgOFFData=cat(2,{data.catMeanWalkedVelAvgOFF});
maxLength = max(cellfun(@length,meanWalkedVelAvgOFFData));
meanWalkedVelAvgOFFMat= cell2mat(cellfun(@(x)cat(1,x,nan(maxLength-length(x),1)),meanWalkedVelAvgOFFData,'UniformOutput',false));
%% Create mean Straightness file for plotting in Graphpad
meanWalkedSpatialStrAvgONData=cat(2,{data.catMeanWalkedSpatialStrAvgON});
maxLength = max(cellfun(@length,meanWalkedSpatialStrAvgONData));
meanWalkedSpatialStrAvgONMat= cell2mat(cellfun(@(x)cat(1,x,nan(maxLength-length(x),1)),meanWalkedSpatialStrAvgONData,'UniformOutput',false));


meanWalkedSpatialStrAvgOFFData=cat(2,{data.catMeanWalkedSpatialStrAvgOFF});
maxLength = max(cellfun(@length,meanWalkedSpatialStrAvgOFFData));
meanWalkedSpatialStrAvgOFFMat= cell2mat(cellfun(@(x)cat(1,x,nan(maxLength-length(x),1)),meanWalkedSpatialStrAvgOFFData,'UniformOutput',false));
%% Create time Walked file for plotting in Graphpad
timeWalkedAvgONData=cat(2,{data.catTimeWalkedAvgON});
maxLength = max(cellfun(@length,timeWalkedAvgONData));
timeWalkedAvgONMat= cell2mat(cellfun(@(x)cat(1,x,nan(maxLength-length(x),1)),timeWalkedAvgONData,'UniformOutput',false));


timeWalkedAvgOFFData=cat(2,{data.catTimeWalkedAvgOFF});
maxLength = max(cellfun(@length,timeWalkedAvgOFFData));
timeWalkedAvgOFFMat= cell2mat(cellfun(@(x)cat(1,x,nan(maxLength-length(x),1)),timeWalkedAvgOFFData,'UniformOutput',false));
%%
%% plot bounded line plots for straightness
gen={data.Genotype};
figure;hold on;%plot average of all trials per genotype and std error bounds


p1=boundedline((1:trialDuration*fps)',mean(data(strcmp({data.Genotype},char(gen(1)))).catSmoothSpatialStrAvgTrial)',std(data(strcmp({data.Genotype},char(gen(1)))).catSmoothSpatialStrAvgTrial)/sqrt(data(strcmp({data.Genotype},char(gen(1)))).numFlies)','m','alpha');
p2=boundedline((1:trialDuration*fps)',mean(data(strcmp({data.Genotype},char(gen(2)))).catSmoothSpatialStrAvgTrial,'omitnan')',std(data(strcmp({data.Genotype},char(gen(2)))).catSmoothSpatialStrAvgTrial,'omitnan')/sqrt(data(strcmp({data.Genotype},char(gen(2)))).numFlies)','b','alpha');
p3=boundedline((1:trialDuration*fps)',mean(data(strcmp({data.Genotype},char(gen(3)))).catSmoothSpatialStrAvgTrial)',std(data(strcmp({data.Genotype},char(gen(3)))).catSmoothSpatialStrAvgTrial)/sqrt(data(strcmp({data.Genotype},char(gen(3)))).numFlies)','r','alpha');
p4=boundedline((1:trialDuration*fps)',mean(data(strcmp({data.Genotype},char(gen(4)))).catSmoothSpatialStrAvgTrial)',std(data(strcmp({data.Genotype},char(gen(4)))).catSmoothSpatialStrAvgTrial)/sqrt(data(strcmp({data.Genotype},char(gen(4)))).numFlies)','g','alpha');
%p5=boundedline((1:trialDuration*fps)',mean(data(strcmp({data.Genotype},char(gen(5)))).catSmoothSpatialStrAvgTrial)',std(data(strcmp({data.Genotype},char(gen(5)))).catSmoothSpatialStrAvgTrial)/sqrt(data(strcmp({data.Genotype},char(gen(5)))).numFlies)','k','alpha');
%p6=boundedline((1:trialDuration*fps)',mean(data(strcmp({data.Genotype},char(gen(6)))).catSmoothSpatialStrAvgTrial,'omitnan')',std(data(strcmp({data.Genotype},char(gen(6)))).catSmoothSpatialStrAvgTrial,'omitnan')/sqrt(data(strcmp({data.Genotype},char(gen(6)))).numFlies)','c','alpha');
%p7=boundedline((1:trialDuration*fps)',mean(data(strcmp({data.Genotype},char(gen(7)))).catSmoothSpatialStrAvgTrial)',std(data(strcmp({data.Genotype},char(gen(7)))).catSmoothSpatialStrAvgTrial)/sqrt(data(strcmp({data.Genotype},char(gen(7)))).numFlies)','r','alpha');
%p8=boundedline((1:trialDuration*fps)',mean(data(strcmp({data.Genotype},char(gen(8)))).catSmoothSpatialStrAvgTrial)',std(data(strcmp({data.Genotype},char(gen(8)))).catSmoothSpatialStrAvgTrial)/sqrt(data(strcmp({data.Genotype},char(gen(8)))).numFlies)','g','alpha');
%p9=boundedline((1:trialDuration*fps)',mean(data(strcmp({data.Genotype},char(gen(9)))).catSmoothSpatialStrAvgTrial)',std(data(strcmp({data.Genotype},char(gen(9)))).catSmoothSpatialStrAvgTrial)/sqrt(data(strcmp({data.Genotype},char(gen(9)))).numFlies)','g','alpha');
%p10=boundedline((1:trialDuration*fps)',mean(data(strcmp({data.Genotype},char(gen(10)))).catSmoothSpatialStrAvgTrial)',std(data(strcmp({data.Genotype},char(gen(10)))).catSmoothSpatialStrAvgTrial)/sqrt(data(strcmp({data.Genotype},char(gen(10)))).numFlies)','k','alpha');

%ylim([0,10]);
line([30*fps,30*fps],ylim,'LineStyle','--','Color','k');
line(xlim,[0,0],'LineStyle','-','Color','k');
set(gca,'YTick', min(ylim):5:max(ylim),'XTick',1:10*fps:max(xlim),'XTickLabel',0:10:max(xlim)/fps);
ylabel('Avg Trial SpatialStrocity +-SEM (mm/s)');xlabel('time (s)');legend([p1,p2,p3,p4],[gen(1),gen(2),gen(3),gen(4)]);
%% plot full length vel

gen={data.Genotype};
figure;hold on;%plot average of all trials per genotype and std error bounds


p1=boundedline((1:maxLastFrameNumber)',mean(data(2).catSmoothVel)',std(data(2).catSmoothVel)/sqrt(data(2).numFlies)','b','alpha');
p1=boundedline((1:maxLastFrameNumber)',mean(data(4).catSmoothVel)',std(data(4).catSmoothVel)/sqrt(data(4).numFlies)','r','alpha');
hold off;
%% plot bounded line plots
gen={data.Genotype};
figure;hold on;%plot average of all trials per genotype and std error bounds


p1=boundedline((1:trialDuration*fps)',mean(data(strcmp({data.Genotype},char(gen(1)))).catSmoothVelAvgTrial)',std(data(strcmp({data.Genotype},char(gen(1)))).catSmoothVelAvgTrial)/sqrt(data(strcmp({data.Genotype},char(gen(1)))).numFlies)','b','alpha');
%p2=boundedline((1:trialDuration*fps)',mean(data(strcmp({data.Genotype},char(gen(2)))).catSmoothVelAvgTrial)',std(data(strcmp({data.Genotype},char(gen(2)))).catSmoothVelAvgTrial)/sqrt(data(strcmp({data.Genotype},char(gen(2)))).numFlies)','b','alpha');
%p3=boundedline((1:trialDuration*fps)',mean(data(strcmp({data.Genotype},char(gen(3)))).catSmoothVelAvgTrial)',std(data(strcmp({data.Genotype},char(gen(3)))).catSmoothVelAvgTrial)/sqrt(data(strcmp({data.Genotype},char(gen(3)))).numFlies)','k','alpha');
%p4=boundedline((1:trialDuration*fps)',mean(data(strcmp({data.Genotype},char(gen(4)))).catSmoothVelAvgTrial)',std(data(strcmp({data.Genotype},char(gen(4)))).catSmoothVelAvgTrial)/sqrt(data(strcmp({data.Genotype},char(gen(4)))).numFlies)','g','alpha');
p5=boundedline((1:trialDuration*fps)',mean(data(strcmp({data.Genotype},char(gen(5)))).catSmoothVelAvgTrial)',std(data(strcmp({data.Genotype},char(gen(5)))).catSmoothVelAvgTrial)/sqrt(data(strcmp({data.Genotype},char(gen(5)))).numFlies)','r','alpha');
%p6=boundedline((1:trialDuration*fps)',mean(data(strcmp({data.Genotype},char(gen(6)))).catSmoothVelAvgTrial)',std(data(strcmp({data.Genotype},char(gen(6)))).catSmoothVelAvgTrial)/sqrt(data(strcmp({data.Genotype},char(gen(6)))).numFlies)','c','alpha');
%p7=boundedline((1:trialDuration*fps)',mean(data(strcmp({data.Genotype},char(gen(7)))).catSmoothVelAvgTrial)',std(data(strcmp({data.Genotype},char(gen(7)))).catSmoothVelAvgTrial)/sqrt(data(strcmp({data.Genotype},char(gen(7)))).numFlies)','r','alpha');
%p8=boundedline((1:trialDuration*fps)',mean(data(strcmp({data.Genotype},char(gen(8)))).catSmoothVelAvgTrial)',std(data(strcmp({data.Genotype},char(gen(8)))).catSmoothVelAvgTrial)/sqrt(data(strcmp({data.Genotype},char(gen(8)))).numFlies)','g','alpha');
%p9=boundedline((1:trialDuration*fps)',mean(data(strcmp({data.Genotype},char(gen(9)))).catSmoothVelAvgTrial)',std(data(strcmp({data.Genotype},char(gen(9)))).catSmoothVelAvgTrial)/sqrt(data(strcmp({data.Genotype},char(gen(9)))).numFlies)','g','alpha');
%p10=boundedline((1:trialDuration*fps)',mean(data(strcmp({data.Genotype},char(gen(10)))).catSmoothVelAvgTrial)',std(data(strcmp({data.Genotype},char(gen(10)))).catSmoothVelAvgTrial)/sqrt(data(strcmp({data.Genotype},char(gen(10)))).numFlies)','k','alpha');

ylim = [0,15];
line([30*fps,30*fps],ylim,'LineStyle','--','Color','k');
line(xlim,[0,0],'LineStyle','-','Color','k');
set(gca,'YTick', min(ylim):5:max(ylim),'XTick',1:10*fps:max(xlim),'XTickLabel',0:10:max(xlim)/fps);
ylabel('Avg Trial Velocity +-SEM (mm/s)');xlabel('time (s)');legend([p1,p5],[gen(1),gen(5)]);

%%
%% plot bounded line plots
gen={data.Genotype};
figure;hold on;%plot average of all trials per genotype and std error bounds


p1=boundedline((1:trialDuration*fps)',mean(data(strcmp({data.Genotype},char(gen(1)))).catSmoothAngVelAvgTrial)',std(data(strcmp({data.Genotype},char(gen(1)))).catSmoothAngVelAvgTrial)/sqrt(data(strcmp({data.Genotype},char(gen(1)))).numFlies)','m','alpha');
%p2=boundedline((1:trialDuration*fps)',mean(data(strcmp({data.Genotype},char(gen(2)))).catSmoothAngVelAvgTrial)',std(data(strcmp({data.Genotype},char(gen(2)))).catSmoothAngVelAvgTrial)/sqrt(data(strcmp({data.Genotype},char(gen(2)))).numFlies)','b','alpha');
p3=boundedline((1:trialDuration*fps)',mean(data(strcmp({data.Genotype},char(gen(3)))).catSmoothAngVelAvgTrial)',std(data(strcmp({data.Genotype},char(gen(3)))).catSmoothAngVelAvgTrial)/sqrt(data(strcmp({data.Genotype},char(gen(3)))).numFlies)','r','alpha');
%p4=boundedline((1:trialDuration*fps)',mean(data(strcmp({data.Genotype},char(gen(4)))).catSmoothAngVelAvgTrial)',std(data(strcmp({data.Genotype},char(gen(4)))).catSmoothAngVelAvgTrial)/sqrt(data(strcmp({data.Genotype},char(gen(4)))).numFlies)','g','alpha');
%p5=boundedline((1:trialDuration*fps)',mean(data(strcmp({data.Genotype},char(gen(5)))).catSmoothAngVelAvgTrial)',std(data(strcmp({data.Genotype},char(gen(5)))).catSmoothAngVelAvgTrial)/sqrt(data(strcmp({data.Genotype},char(gen(5)))).numFlies)','k','alpha');
%p6=boundedline((1:trialDuration*fps)',mean(data(strcmp({data.Genotype},char(gen(6)))).catSmoothAngVelAvgTrial)',std(data(strcmp({data.Genotype},char(gen(6)))).catSmoothAngVelAvgTrial)/sqrt(data(strcmp({data.Genotype},char(gen(6)))).numFlies)','c','alpha');
%p7=boundedline((1:trialDuration*fps)',mean(data(strcmp({data.Genotype},char(gen(7)))).catSmoothAngVelAvgTrial)',std(data(strcmp({data.Genotype},char(gen(7)))).catSmoothAngVelAvgTrial)/sqrt(data(strcmp({data.Genotype},char(gen(7)))).numFlies)','r','alpha');
%p8=boundedline((1:trialDuration*fps)',mean(data(strcmp({data.Genotype},char(gen(8)))).catSmoothAngVelAvgTrial)',std(data(strcmp({data.Genotype},char(gen(8)))).catSmoothAngVelAvgTrial)/sqrt(data(strcmp({data.Genotype},char(gen(8)))).numFlies)','g','alpha');
%p9=boundedline((1:trialDuration*fps)',mean(data(strcmp({data.Genotype},char(gen(9)))).catSmoothAngVelAvgTrial)',std(data(strcmp({data.Genotype},char(gen(9)))).catSmoothAngVelAvgTrial)/sqrt(data(strcmp({data.Genotype},char(gen(9)))).numFlies)','g','alpha');
%p10=boundedline((1:trialDuration*fps)',mean(data(strcmp({data.Genotype},char(gen(10)))).catSmoothAngVelAvgTrial)',std(data(strcmp({data.Genotype},char(gen(10)))).catSmoothAngVelAvgTrial)/sqrt(data(strcmp({data.Genotype},char(gen(10)))).numFlies)','k','alpha');

ylim = [0,5];
line([30*fps,30*fps],ylim,'LineStyle','--','Color','k');
line(xlim,[0,0],'LineStyle','-','Color','k');
set(gca,'YTick', min(ylim):5:max(ylim),'XTick',1:10*fps:max(xlim),'XTickLabel',0:10:max(xlim)/fps);
ylabel('Avg Trial Ang Vel +-SEM (r/s)');xlabel('time (s)');legend([p1,p3],[gen(1),gen(3)]);

%%
subplot(2,1,2);
hold on;
p1=boundedline((1:trialDuration*fps)',mean(data(strcmp({data.Genotype},char(gen(1)))).catSmoothAngVelSignedAvgTrial)',std(data(strcmp({data.Genotype},char(gen(1)))).catSmoothAngVelSignedAvgTrial)/sqrt(data(strcmp({data.Genotype},char(gen(1)))).numFlies)','m');
p2=boundedline((1:trialDuration*fps)',mean(data(strcmp({data.Genotype},char(gen(2)))).catSmoothAngVelSignedAvgTrial)',std(data(strcmp({data.Genotype},char(gen(2)))).catSmoothAngVelSignedAvgTrial)/sqrt(data(strcmp({data.Genotype},char(gen(2)))).numFlies)','b','alpha');
p3=boundedline((1:trialDuration*fps)',mean(data(strcmp({data.Genotype},char(gen(3)))).catSmoothAngVelSignedAvgTrial)',std(data(strcmp({data.Genotype},char(gen(3)))).catSmoothAngVelSignedAvgTrial)/sqrt(data(strcmp({data.Genotype},char(gen(3)))).numFlies)','r','alpha');
p4=boundedline((1:trialDuration*fps)',mean(data(strcmp({data.Genotype},char(gen(4)))).catSmoothAngVelSignedAvgTrial)',std(data(strcmp({data.Genotype},char(gen(4)))).catSmoothAngVelSignedAvgTrial)/sqrt(data(strcmp({data.Genotype},char(gen(4)))).numFlies)','g','alpha');
%p5=boundedline((1:trialDuration*fps)',mean(data(strcmp({data.Genotype},char(gen(5)))).catSmoothAngVelSignedAvgTrial)',std(data(strcmp({data.Genotype},char(gen(5)))).catSmoothAngVelSignedAvgTrial)/sqrt(data(strcmp({data.Genotype},char(gen(5)))).numFlies)','m');
%p6=boundedline((1:trialDuration*fps)',mean(data(strcmp({data.Genotype},char(gen(6)))).catSmoothAngVelSignedAvgTrial)',std(data(strcmp({data.Genotype},char(gen(6)))).catSmoothAngVelSignedAvgTrial)/sqrt(data(strcmp({data.Genotype},char(gen(6)))).numFlies)','b','alpha');
%p7=boundedline((1:trialDuration*fps)',mean(data(strcmp({data.Genotype},char(gen(7)))).catSmoothAngVelSignedAvgTrial)',std(data(strcmp({data.Genotype},char(gen(7)))).catSmoothAngVelSignedAvgTrial)/sqrt(data(strcmp({data.Genotype},char(gen(7)))).numFlies)','r','alpha');
%p9=boundedline((1:trialDuration*fps)',mean(data(strcmp({data.Genotype},char(gen(8)))).catSmoothAngVelSignedAvgTrial)',std(data(strcmp({data.Genotype},char(gen(8)))).catSmoothAngVelSignedAvgTrial)/sqrt(data(strcmp({data.Genotype},char(gen(8)))).numFlies)','g','alpha');

ylim = [-10,10];
line([30*fps,30*fps],ylim,'LineStyle','--','Color','k');
line(xlim,[0,0],'LineStyle','-','Color','k');
set(gca,'YTick', min(ylim):5:max(ylim),'XTick',1:10*fps:max(xlim),'XTickLabel',0:10:max(xlim)/fps);
%ylabel('Avg Trial Angular Velocity SEM (rad/s)');xlabel('time (s)');%legend([p1,p2,p3],{'left','right','both'});
hold off;
%% plot single fly velocity plots
figure;hold on;
h2=subplot(2,1,1);
hold on;
p1=plot((1:trialDuration*fps)',data(strcmp({data.Genotype},'left')).catSmoothVelAvgTrial,'b','LineWidth',0.5)';
p2=plot((1:trialDuration*fps)',data(strcmp({data.Genotype},'right')).catSmoothVelAvgTrial,'r','LineWidth',0.5)';
p3=plot((1:trialDuration*fps)',data(strcmp({data.Genotype},'both')).catSmoothVelAvgTrial,'g','LineWidth',0.5)';
ylim=[0,25];
line([30*fps,30*fps],ylim,'LineStyle','--','Color','k');
line(xlim,[0,0],'LineStyle','-','Color','k');
set(gca,'YTick', min(ylim):5:max(ylim),'XTick',1:10*fps:max(xlim),'XTickLabel',0:10:max(xlim)/fps);
ylabel('Avg Trial Velocity (mm/s)');xlabel('time (s)');%legend([p1(1),p2(1),p3(1)],{'left','right','both'});

h2=subplot(2,1,2);
hold on;
p1=plot((1:trialDuration*fps)',data(strcmp({data.Genotype},'left')).catSmoothAngVelSignedAvgTrial,'b','LineWidth',0.5)';
p2=plot((1:trialDuration*fps)',data(strcmp({data.Genotype},'right')).catSmoothAngVelSignedAvgTrial,'r','LineWidth',0.5)';
p3=plot((1:trialDuration*fps)',data(strcmp({data.Genotype},'both')).catSmoothAngVelSignedAvgTrial,'g','LineWidth',0.5)';
ylim=[-15,15];
line([30*fps,30*fps],ylim,'LineStyle','--','Color','k');
line(xlim,[0,0],'LineStyle','-','Color','k');
set(gca,'YTick', min(ylim):5:max(ylim),'XTick',1:10*fps:max(xlim),'XTickLabel',0:10:max(xlim)/fps);
ylabel('Avg Trial Velocity (rad/s)');xlabel('time (s)');%legend([p1(1),p2(1),p3(1)],{'left','right','both'});
hold off;

%%

% figure;hold on;%plot average of all trials per genotype and std dev bounds
% p1=boundedline((1:trialDuration*fps)',mean(data(strcmp({data.Genotype},[GenotypeName '-male'])).catSmoothVelAllTrials)',std(data(strcmp({data.Genotype},[GenotypeName '-male'])).catSmoothVelAllTrials)','b','alpha');
% p2=boundedline((1:trialDuration*fps)',mean(data(strcmp({data.Genotype},[GenotypeName '-female'])).catSmoothVelAllTrials)',std(data(strcmp({data.Genotype},[GenotypeName '-female'])).catSmoothVelAllTrials)','r','alpha');
% p3=boundedline((1:trialDuration*fps)',mean(data(strcmp({data.Genotype},'w1118-male')).catSmoothVelAllTrials)',std(data(strcmp({data.Genotype},'w1118-male')).catSmoothVelAllTrials)','g','alpha');
% p4=boundedline((1:trialDuration*fps)',mean(data(strcmp({data.Genotype},'w1118-female')).catSmoothVelAllTrials)',std(data(strcmp({data.Genotype},'w1118-female')).catSmoothVelAllTrials)','m','alpha');
% line([30*fps,30*fps],ylim,'LineStyle','--','Color','k');
% line(xlim,[0,0],'LineStyle','-','Color','k');
% set(gca,'YTick', min(ylim):5:max(ylim),'XTick',1:10*fps:max(xlim),'XTickLabel',0:10:max(xlim)/fps);
% ylabel('Velocity +-SD (mm/s)');xlabel('time (s)');legend([p1,p2,p3,p4],{[GenotypeName '-male'],[GenotypeName '-female'],'w1118-male','w1118-female'});
% hold off;

%% Angular Velocity
GenotypeName = 'Mesa';
figure;hold on;%plot average of all trials per genotype and std error bounds
p1=boundedline((1:trialDuration*fps)',mean(data(strcmp({data.Genotype},'left')).catSmoothAngVelAvgTrial)',std(data(strcmp({data.Genotype},'left')).catSmoothAngVelAvgTrial)/sqrt(data(strcmp({data.Genotype},'left')).numFlies)','b');
p2=boundedline((1:trialDuration*fps)',mean(data(strcmp({data.Genotype},'right')).catSmoothAngVelAvgTrial)',std(data(strcmp({data.Genotype},'right')).catSmoothAngVelAvgTrial)/sqrt(data(strcmp({data.Genotype},'right')).numFlies)','r','alpha');
p3=boundedline((1:trialDuration*fps)',mean(data(strcmp({data.Genotype},'both')).catSmoothAngVelAvgTrial)',std(data(strcmp({data.Genotype},'both')).catSmoothAngVelAvgTrial)/sqrt(data(strcmp({data.Genotype},'both')).numFlies)','g','alpha');

line([30*fps,30*fps],ylim,'LineStyle','--','Color','k');
line(xlim,[0,0],'LineStyle','-','Color','k');
set(gca,'YTick', min(ylim):1:max(ylim),'XTick',1:10*fps:max(xlim),'XTickLabel',0:10:max(xlim)/fps);
ylabel('Avg Trial Angular Velocity SEM (mm/s)');xlabel('time (s)');legend([p1,p2,p3],{'left','right','both'});
hold off;

%% Angular Velocity Signed
%% Angular Velocity
GenotypeName = 'Mesa';
figure;hold on;%plot average of all trials per genotype and std error bounds
p1=boundedline((1:trialDuration*fps)',mean(data(strcmp({data.Genotype},'left')).catSmoothAngVelSignedAvgTrial)',std(data(strcmp({data.Genotype},'left')).catSmoothAngVelSignedAvgTrial)/sqrt(data(strcmp({data.Genotype},'left')).numFlies)','b');
p2=boundedline((1:trialDuration*fps)',mean(data(strcmp({data.Genotype},'right')).catSmoothAngVelSignedAvgTrial)',std(data(strcmp({data.Genotype},'right')).catSmoothAngVelSignedAvgTrial)/sqrt(data(strcmp({data.Genotype},'right')).numFlies)','r','alpha');
p3=boundedline((1:trialDuration*fps)',mean(data(strcmp({data.Genotype},'both')).catSmoothAngVelSignedAvgTrial)',std(data(strcmp({data.Genotype},'both')).catSmoothAngVelSignedAvgTrial)/sqrt(data(strcmp({data.Genotype},'both')).numFlies)','g','alpha');

line([30*fps,30*fps],ylim,'LineStyle','--','Color','k');
line(xlim,[0,0],'LineStyle','-','Color','k');
set(gca,'YTick', min(ylim):1:max(ylim),'XTick',1:10*fps:max(xlim),'XTickLabel',0:10:max(xlim)/fps);
ylabel('Avg Trial Angular Velocity SEM (rad/s)');xlabel('time (s)');legend([p1,p2,p3],{'left','right','both'});
hold off;

