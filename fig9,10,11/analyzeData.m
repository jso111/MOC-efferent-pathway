% Program to load in Ari's data and create time-stamped variables 
% of the awake mice data
% Pupil diameter, gain, Q10dB, mag, phase, etc. 

% This version is modified to bin pupil diameter into different categories

close all;clear;clc;
warning('off','all')

% conditions={'Awake','Anesth','Dead'};
conditions={'Awake','Anesth'};
% conditions={'Awake'};
protocolLabels={'*multitone*','*single-tone*','*TM*'};
% protocolLabels={'*multitone*','*TM*'};
% protocolLabels={'*TM*'};
protocols={'Multitone','Singletone','TM'};
% protocols={'Multitone','TM'};
% protocols={'TM'};
genotypes={'WT','Alpha9KO'};
pupilSmooth=2;
pSize={'Small','Medium','Large'};
percent1=20; %the top and bottom percent of pupil diameters for thresholding
percent2=50; %the top and bottom percent of pupil diameter derivatives for thresholding

useDerivative=0;
%%
%dataAll=loadData(conditions,protocolLabels,protocols,genotypes,pupilSmooth);
load dataAll;
%%
pupilSizeAll=calc_pupilSizes(dataAll,percent1, percent2);

%%
binWidth=10; % the percent of max pupil diameter used for bin width

if length(dataAll)>5
    p1=2;
    p2=6;
elseif length(dataAll)>1
    p1=5;
    p2=1;
else
    p1=1;
    p2=1;
end

delayTime=0; %time delay between vibrometry measurement and pupil measurement (positive is pupil after vibrometry)
minNum=1; %minimum number of points to include in averaging
close all;
v1=table('Size',[1 14], ...
    'VariableTypes',{'int32','string','string','string','string','double', ...
    'double','double','double','double', 'double', ...
    'double','double','double'},'VariableNames', ...
    {'experiment','condition','protocol','genotype','pupilSize', ...
    'freq','level','mag','phi','gain1','gain2','q','bf','diameter'});
v2=table('Size',[1 14], ...
    'VariableTypes',{'int32','string','string','string','string','double', ...
    'double','double','double','double', 'double', ...
    'double','double','double'},'VariableNames', ...
    {'experiment','condition','protocol','genotype','pupilSize', ...
    'freq','level','mag','phi','gain1','gain2','q','bf','diameter'});
v3=table('Size',[1 14], ...
    'VariableTypes',{'int32','string','string','string','string','double', ...
    'double','double','double','double', 'double', ...
    'double','double','double'},'VariableNames', ...
    {'experiment','condition','protocol','genotype','pupilSize', ...
    'freq','level','mag','phi','gain1','gain2','q','bf','diameter'});
v4=table('Size',[1 14], ...
    'VariableTypes',{'int32','string','string','string','string','double', ...
    'double','double','double','double', 'double', ...
    'double','double','double'},'VariableNames', ...
    {'experiment','condition','protocol','genotype','pupilSize', ...
    'freq','level','mag','phi','gain1','gain2','q','bf','diameter'});

v1AA=table('Size',[1 14], ...
    'VariableTypes',{'int32','string','string','string','string','double', ...
    'double','double','double','double', 'double', ...
    'double','double','double'},'VariableNames', ...
    {'experiment','condition','protocol','genotype','pupilSize', ...
    'freq','level','mag','phi','gain1','gain2','q','bf','diameter'});
v2AA=table('Size',[1 14], ...
    'VariableTypes',{'int32','string','string','string','string','double', ...
    'double','double','double','double', 'double', ...
    'double','double','double'},'VariableNames', ...
    {'experiment','condition','protocol','genotype','pupilSize', ...
    'freq','level','mag','phi','gain1','gain2','q','bf','diameter'});
v3AA=table('Size',[1 14], ...
    'VariableTypes',{'int32','string','string','string','string','double', ...
    'double','double','double','double', 'double', ...
    'double','double','double'},'VariableNames', ...
    {'experiment','condition','protocol','genotype','pupilSize', ...
    'freq','level','mag','phi','gain1','gain2','q','bf','diameter'});

iTable=1;

% now calculate averaged data for binned pupil diameters
% for each experiment
threshold=10; % must have at least this many measurements in a bin to use the averaged data

nBins=floor(100/binWidth);
xBin=zeros(1,nBins);
yMagAve=zeros(1,nBins);
yMagLowAve=zeros(1,nBins);
yPhaseAve=zeros(1,nBins);
yGainAve=zeros(1,nBins);
yBFAve=zeros(1,nBins);
yQAve=zeros(1,nBins);

yMagAveData=zeros(length(dataAll),nBins);
yMagLowData=zeros(length(dataAll),nBins);
yPhaseData=zeros(length(dataAll),nBins);
yGainAveData=zeros(length(dataAll),nBins);
yBFAveData=zeros(length(dataAll),nBins);
yQAveData=zeros(length(dataAll),nBins);

for c=1:(length(conditions))  
    for p=1:(length(protocols))   
        for e=1:length(dataAll)
            pupilSize=pupilSizeAll{e};
            clear mag phi gain1 gain 2 Q BF diameter yMag yMagLow yPhase yGain yQ yBF;
            logicalp=strcmp({dataAll(e).vib.protocol},protocols{p});
            logicalc=strcmp({dataAll(e).vib.condition}, conditions{c});
            logicalcp=logicalp & logicalc;
            index=find(logicalcp);
            nAve=sum(logicalcp);
            F=dataAll(e).vib(1).F;
            L=dataAll(e).vib(1).L;
            if nAve>0
                diameter=zeros(1,nAve);
                pupilBin=zeros(1,nAve);
                derivative=zeros(1,nAve);
                timeMeasure=zeros(1,nAve);
                gain1=zeros(length(F),nAve);
                gain2=zeros(length(F),nAve);
                Q=zeros(length(L),nAve);
                BF=zeros(length(L),nAve);
                mag=zeros(length(F),length(L),nAve);
                phi=zeros(length(F),length(L),nAve);
                for n=1:nAve
                    mag(:,:,n)=dataAll(e).vib(index(n)).mag;
                    phi(:,:,n)=dataAll(e).vib(index(n)).phi;
                    gain1(:,n)=dataAll(e).vib(index(n)).gain1;
                    gain2(:,n)=dataAll(e).vib(index(n)).gain2;
                    Q(:,n)=dataAll(e).vib(index(n)).q;
                    BF(:,n)=dataAll(e).vib(index(n)).bf;

                    % get data at the BF for the stim level we are looking at
                    levelIndex=2; %20 dB
                    levelIndexLow=5; %50 dB
                    BFIndex=find(F==max(BF(:,n)));
                    if isempty(BFIndex)
                        yMag(n)=nan;
                        yMagLow(n)=nan;
                        yPhase(n)=nan;
                        yGain(n)=nan;
                        yBF(n)=nan;
                        yQ(n)=nan;
                    else
                        yMag(n)=mag(BFIndex,levelIndex,n);
                        yMagLow(n)=mag(floor(BFIndex/2),levelIndexLow,n);
                        yPhase(n)=phi(BFIndex,levelIndex,n);
                        yGain(n)=gain1(BFIndex,n);
                        yBF(n)=BF(levelIndex,n);
                        yQ(n)=Q(levelIndex,n);
                    end
         
                    %find pupil diameter and derivative
                    timeMeasure(n)=seconds(dataAll(e).vib(index(n)).time-dataAll(e).vib(1).time);
                    tDiff=dataAll(e).pupil.datetime - (dataAll(e).vib(index(n)).time + duration(0,0,delayTime));
                    [x,tDiffIndex]=min(abs(tDiff));               
                    if seconds(x)>2 
                        % if pupil diameter was not measured at the time of the
                        % vibratory measurement, set diameter and derivative to nan
                        diameter(n)=nan;
                        derivative(n)=nan;
                    else
                        diameter(n)=dataAll(e).pupil.diameter(tDiffIndex);
                        derivative(n)=dataAll(e).pupil.derivative(tDiffIndex);
                    end
                end   

                % bin the pupil dimeter
                pupilBin=floor((100*diameter/max(diameter))/binWidth);
                dataScatter(c,p,e).nAve=nAve;
                dataScatter(c,p,e).timeMeasure=timeMeasure;
                dataScatter(c,p,e).diameter=diameter;
                dataScatter(c,p,e).pupilBin=pupilBin;
                dataScatter(c,p,e).derivative=derivative;
                dataScatter(c,p,e).gain1=gain1;
                dataScatter(c,p,e).gain2=gain2;
                dataScatter(c,p,e).Q=Q;
                dataScatter(c,p,e).BF=BF;
                dataScatter(c,p,e).yMag=yMag;
                dataScatter(c,p,e).goodDiam=~isnan(diameter);
                dataScatter(c,p,e).goodMag=~isnan(yMag);
                dataScatter(c,p,e).goodboth=dataScatter(c,p,e).goodDiam & dataScatter(c,p,e).goodMag;
                dataScatter(c,p,e).nGood=sum(dataScatter(c,p,e).goodboth);

                for i=1:nBins
                    xBin(i)=(i-1)*binWidth;
                    if (sum(~isnan(yMag(pupilBin==i)))>threshold)
                        yMagAve(i)=mean(yMag(pupilBin==i),'omitnan');
                    else
                        yMagAve(i)=nan;
                    end
                    if (sum(~isnan(yMagLow(pupilBin==i)))>threshold)
                        yMagLowAve(i)=mean(yMagLow(pupilBin==i),'omitnan');
                    else
                        yMagLowAve(i)=nan;
                    end
                    if (sum(~isnan(yPhase(pupilBin==i)))>threshold)
                        yPhaseAve(i)=mean(yPhase(pupilBin==i),'omitnan');
                    else
                        yPhaseAve(i)=nan;
                    end
                    if (sum(~isnan(yGain(pupilBin==i)))>threshold)
                        yGainAve(i)=mean(yGain(pupilBin==i),'omitnan');
                    else
                        yGainAve(i)=nan;
                    end
                    if (sum(~isnan(yBF(pupilBin==i)))>threshold)
                        yBFAve(i)=mean(yBF(pupilBin==i),'omitnan');
                    else
                        yBFAve(i)=nan;
                    end
                    if (sum(~isnan(yQ(pupilBin==i)))>threshold)
                        yQAve(i)=mean(yQ(pupilBin==i),'omitnan');
                    else
                        yQAve(i)=nan;
                    end                   
                end
                dataScatter(c,p,e).xBin=xBin;
                dataScatter(c,p,e).yMagAve=yMagAve;
                dataScatter(c,p,e).yGainAve=yGainAve;
                dataScatter(c,p,e).BFAve=mean(BF','omitnan');
                dataScatter(c,p,e).QAve=mean(Q','omitnan');

                if (c==1&p==1)
                    % Load data into matrix for later averaging for binned
                    % histograms
                    yMagAveData(e,:)=yMagAve;
                    yMagLowAveData(e,:)=yMagLowAve;
                    yPhaseAveData(e,:)=yPhaseAve;
                    yGainAveData(e,:)=yGainAve;
                    yBFAveData(e,:)=yBFAve;
                    yQAveData(e,:)=yQAve;
                end

                if (c==1|c==2)&(p==1)
                    % Load data into arrays for awake vs anesth comparisons
                    dataAA(c,e).F=F;
                    dataAA(c,e).Fmat=repmat(F,[length(L) 1])';
                    dataAA(c,e).L=L;

                    dataAA(c,e).mag=trimmean(mag,5,3);
                    dataAA(c,e).phi=mean(phi,3,'omitnan');
                    dataAA(c,e).gain1=mean(gain1,2,'omitnan');
                    dataAA(c,e).gain2=mean(gain2,2,'omitnan');
                    dataAA(c,e).Q=mean(Q,2,'omitnan');
                    dataAA(c,e).BF=mean(BF,2,'omitnan');
    
                    dataAA(c,e).magSD=std(mag,0,3,'omitnan');
                    dataAA(c,e).phiSD=std(phi,0,3,'omitnan');
                    dataAA(c,e).gain1SD=std(gain1,0,2,'omitnan');
                    dataAA(c,e).gain2SD=std(gain2,0,2,'omitnan');
                    dataAA(c,e).QSD=std(Q,0,2,'omitnan');
                    dataAA(c,e).BFSD=std(BF,0,2,'omitnan');
    
                    dataAA(c,e).magN=length(mag) - sum(isnan(mag),3);
                    dataAA(c,e).phiN=length(phi) - sum(isnan(phi),3);
                    dataAA(c,e).gain1N=length(gain1) - sum(isnan(gain1),2);
                    dataAA(c,e).gain2N=length(gain2) - sum(isnan(gain2),2);
                    dataAA(c,e).QN=length(Q) - sum(isnan(Q),2);
                    dataAA(c,e).BFN=length(BF) - sum(isnan(BF),2);
                end
                
                % for averaging based on pupil size, remove datapoints where the pupil was not
                % measured
                indNAN=isnan(diameter);
                mag(:,:,indNAN) =[];
                phi(:,:,indNAN) =[];
                gain1(:,indNAN) =[];
                gain2(:,indNAN) =[];
                Q(:,indNAN) =[];
                BF(:,indNAN) =[];
                timeMeasure(indNAN) =[];
                diameter(indNAN) =[];
                derivative(indNAN) =[];   

                if useDerivative==0
                    % average the data into 3 groups by pupil size
                    idx{1}=find(diameter<=pupilSize(1));
                    idx{2}=find(diameter>pupilSize(1)&diameter<pupilSize(2));
                    idx{3}=find(diameter>=pupilSize(2));
                elseif useDerivative==1
                    % average the data into 3 groups by pupil size
                    idx{1}=find(diameter>pupilSize(1)&diameter<pupilSize(2)&derivative<=pupilSize(3));
                    idx{2}=find(diameter>pupilSize(1)&diameter<pupilSize(2)&derivative>pupilSize(3)&derivative<pupilSize(4));
                    idx{3}=find(diameter>pupilSize(1)&diameter<pupilSize(2)&derivative>=pupilSize(4));
                end
                
                for s=1:3
                    data(c,p,e,s).size=pSize(s);
                    data(c,p,e,s).F=F;
                    data(c,p,e,s).Fmat=repmat(F,[length(L) 1])';
                    data(c,p,e,s).L=L;

                    data(c,p,e,s).mag=trimmean(mag(:,:,idx{s}),5,3);
                    data(c,p,e,s).phi=mean(phi(:,:,idx{s}),3,'omitnan');
                    data(c,p,e,s).gain1=mean(gain1(:,idx{s}),2,'omitnan');
                    data(c,p,e,s).gain2=mean(gain2(:,idx{s}),2,'omitnan');
                    data(c,p,e,s).Q=mean(Q(:,idx{s}),2,'omitnan');
                    data(c,p,e,s).BF=mean(BF(:,idx{s}),2,'omitnan');
                    data(c,p,e,s).diameter=mean(diameter(idx{s}),'omitnan');
    
                    data(c,p,e,s).magSD=std(mag(:,:,idx{s}),0,3,'omitnan');
                    data(c,p,e,s).phiSD=std(phi(:,:,idx{s}),0,3,'omitnan');
                    data(c,p,e,s).gain1SD=std(gain1(:,idx{s}),0,2,'omitnan');
                    data(c,p,e,s).gain2SD=std(gain2(:,idx{s}),0,2,'omitnan');
                    data(c,p,e,s).QSD=std(Q(:,idx{s}),0,2,'omitnan');
                    data(c,p,e,s).BFSD=std(BF(:,idx{s}),0,2,'omitnan');
                    data(c,p,e,s).diameterSD=std(diameter(idx{s}),'omitnan');
    
                    data(c,p,e,s).magN=length(mag(:,:,idx{s})) - sum(isnan(mag(:,:,idx{s})),3);
                    data(c,p,e,s).phiN=length(phi(:,:,idx{s})) - sum(isnan(phi(:,:,idx{s})),3);
                    data(c,p,e,s).gain1N=length(gain1(:,idx{s})) - sum(isnan(gain1(:,idx{s})),2);
                    data(c,p,e,s).gain2N=length(gain2(:,idx{s})) - sum(isnan(gain2(:,idx{s})),2);
                    data(c,p,e,s).QN=length(Q(:,idx{s})) - sum(isnan(Q(:,idx{s})),2);
                    data(c,p,e,s).BFN=length(BF(:,idx{s})) - sum(isnan(BF(:,idx{s})),2);
                    data(c,p,e,s).diameterN=length(diameter(:,idx{s})) - sum(isnan(diameter(:,idx{s})),2);
                end
            end
        end

        % ave data for pupil size comparison
        for g=1:2
            logicalg=strcmp({dataAll.genotype}, genotypes{g});
            indexg=find(logicalg);
            nAve=sum(logicalg);
            tempmag=cell(1,3);
            tempphi=cell(1,3);
            tempgain1=cell(1,3);
            tempgain2=cell(1,3);
            tempQ=cell(1,3);
            tempBF=cell(1,3);
            tempdiameter=cell(1,3);            
            for s=1:3
                tempmag{s}=nan(20,8,nAve);
                tempphi{s}=nan(20,8,nAve);
                tempgain1{s}=nan(20,nAve);
                tempgain2{s}=nan(20,nAve);
                tempQ{s}=nan(8,nAve);
                tempBF{s}=nan(8,nAve);
                tempdiameter{s}=nan(1,nAve);
                for n=1:nAve
                    z=data(c,p,indexg(n),s).mag;
                    if length(z)>0
                        tempmag{s}(:,:,n)=z;
                    end
                    z=data(c,p,indexg(n),s).phi;
                    if length(z)>0
                        tempphi{s}(:,:,n)=z;
                    end
                    z=data(c,p,indexg(n),s).gain1;
                    if length(z)>0
                        tempgain1{s}(:,n)=z;
                    end
                    z=data(c,p,indexg(n),s).gain2;
                    if length(z)>0
                        tempgain2{s}(:,n)=z;
                    end
                    z=data(c,p,indexg(n),s).Q;
                    if length(z)>0
                        tempQ{s}(:,n)=z;
                    end
                    z=data(c,p,indexg(n),s).BF;
                    if length(z)>0
                        tempBF{s}(:,n)=z;
                    end
                    z=data(c,p,indexg(n),s).diameter;
                    if length(z)>0
                        tempdiameter{s}(n)=z;
                    end

                    % now fill up a table to statistically compare small vs large pupils
                    % for each genotype using a multiple linear regression, paired testing approach
                    % (linear mixed model in R) - create a huge spreadsheet of average data from each mouse:
                    % experiment	condition	protocol	genotype	pupilSize	freq	level	
                    % mag	phi	gain1	gain2 q	bf	diameter

                    % mag/phase
                    temp=reshape(tempmag{s}(:,:,n),[],1);
                    npts=length(temp);
                    v1Temp=table('Size',[npts 14], ...
                        'VariableTypes',{'int32','string','string','string','string','double', ...
                        'double','double','double','double', 'double', ...
                        'double','double','double'},'VariableNames', ...
                        {'experiment','condition','protocol','genotype','pupilSize', ...
                        'freq','level','mag','phi','gain1','gain2','q','bf','diameter'});
                    v1Temp.experiment=repmat(indexg(n),[npts 1]);
                    v1Temp.condition=repmat(conditions(c),[npts 1]);
                    v1Temp.protocol=repmat(protocols(p),[npts 1]);
                    v1Temp.genotype=repmat(genotypes(g),[npts 1]);
                    v1Temp.pupilSize=repmat(pSize(s),[npts 1]);
                    v1Temp.freq=reshape(repmat(F,[length(L) 1])',[],1);
                    v1Temp.level=reshape(repmat(L,[length(F) 1]),[],1);
                    v1Temp.mag=temp;
                    v1Temp.phi=reshape(tempphi{s}(:,:,n),[],1);

                    % gain1/2
                    temp=reshape(tempgain1{s}(:,n),[],1);
                    npts=length(temp);
                    v2Temp=table('Size',[npts 14], ...
                        'VariableTypes',{'int32','string','string','string','string','double', ...
                        'double','double','double','double', 'double', ...
                        'double','double','double'},'VariableNames', ...
                        {'experiment','condition','protocol','genotype','pupilSize', ...
                        'freq','level','mag','phi','gain1','gain2','q','bf','diameter'});
                    v2Temp.experiment=repmat(indexg(n),[npts 1]);
                    v2Temp.condition=repmat(conditions(c),[npts 1]);
                    v2Temp.protocol=repmat(protocols(p),[npts 1]);
                    v2Temp.genotype=repmat(genotypes(g),[npts 1]);
                    v2Temp.pupilSize=repmat(pSize(s),[npts 1]);
                    v2Temp.freq=F';
                    v2Temp.gain1=temp;
                    v2Temp.gain2=reshape(tempgain2{s}(:,n),[],1);

                    % q, bf
                    temp=reshape(tempQ{s}(:,n),[],1);
                    npts=length(temp);
                    v3Temp=table('Size',[npts 14], ...
                        'VariableTypes',{'int32','string','string','string','string','double', ...
                        'double','double','double','double', 'double', ...
                        'double','double','double'},'VariableNames', ...
                        {'experiment','condition','protocol','genotype','pupilSize', ...
                        'freq','level','mag','phi','gain1','gain2','q','bf','diameter'});
                    v3Temp.experiment=repmat(indexg(n),[npts 1]);
                    v3Temp.condition=repmat(conditions(c),[npts 1]);
                    v3Temp.protocol=repmat(protocols(p),[npts 1]);
                    v3Temp.genotype=repmat(genotypes(g),[npts 1]);
                    v3Temp.pupilSize=repmat(pSize(s),[npts 1]);
                    v3Temp.level=L';
                    v3Temp.q=temp;
                    v3Temp.bf=reshape(tempBF{s}(:,n),[],1);

                    % diameter
                    temp=tempdiameter{s}(n);
                    npts=length(temp);
                    v4Temp=table('Size',[npts 14], ...
                        'VariableTypes',{'int32','string','string','string','string','double', ...
                        'double','double','double','double', 'double', ...
                        'double','double','double'},'VariableNames', ...
                        {'experiment','condition','protocol','genotype','pupilSize', ...
                        'freq','level','mag','phi','gain1','gain2','q','bf','diameter'});
                    v4Temp.experiment=repmat(indexg(n),[npts 1]);
                    v4Temp.condition=repmat(conditions(c),[npts 1]);
                    v4Temp.protocol=repmat(protocols(p),[npts 1]);
                    v4Temp.genotype=repmat(genotypes(g),[npts 1]);
                    v4Temp.pupilSize=repmat(pSize(s),[npts 1]);
                    v4Temp.diameter=temp;

                    v1=[v1;v1Temp];
                    v2=[v2;v2Temp];
                    v3=[v3;v3Temp];
                    v4=[v4;v4Temp];

                end
                dataAve1(c,p,g,s).mag=mean(tempmag{s},3,'omitnan');
                dataAve1(c,p,g,s).magSD=std(tempmag{s},0,3,'omitnan');
                dataAve1(c,p,g,s).magN=size(tempmag{s},3) - sum(isnan(tempmag{s}),3);
                dataAve1(c,p,g,s).phi=mean(tempphi{s},3,'omitnan');
                dataAve1(c,p,g,s).phiSD=std(tempphi{s},0,3,'omitnan');
                dataAve1(c,p,g,s).phiN=size(tempphi{s},3) - sum(isnan(tempphi{s}),3);
                dataAve1(c,p,g,s).gain1=mean(tempgain1{s},2,'omitnan');
                dataAve1(c,p,g,s).gain1SD=std(tempgain1{s},0,2,'omitnan');
                dataAve1(c,p,g,s).gain1N=size(tempgain1{s},2) - sum(isnan(tempgain1{s}),2);
                dataAve1(c,p,g,s).gain2=mean(tempgain2{s},2,'omitnan');
                dataAve1(c,p,g,s).gain2SD=std(tempgain2{s},0,2,'omitnan');
                dataAve1(c,p,g,s).gain2N=size(tempgain2{s},2) - sum(isnan(tempgain2{s}),2);
                dataAve1(c,p,g,s).Q=mean(tempQ{s},2,'omitnan');
                dataAve1(c,p,g,s).QSD=std(tempQ{s},0,2,'omitnan');
                dataAve1(c,p,g,s).QN=size(tempQ{s},2) - sum(isnan(tempQ{s}),2);
                dataAve1(c,p,g,s).BF=mean(tempBF{s},2,'omitnan');
                dataAve1(c,p,g,s).BFSD=std(tempBF{s},0,2,'omitnan');
                dataAve1(c,p,g,s).BFN=size(tempBF{s},2) - sum(isnan(tempBF{s}),2);
                dataAve1(c,p,g,s).diameter=mean(tempdiameter{s},1,'omitnan');
                dataAve1(c,p,g,s).diameterSD=std(tempdiameter{s},0,1,'omitnan');
                dataAve1(c,p,g,s).diameterN=size(tempdiameter{s},1) - sum(isnan(tempdiameter{s}),1);
            end
        end

        % Ave data for awake/anesth comparison
        if (c==1|c==2)&(p==1)
            for g=1:2
                logicalg=strcmp({dataAll.genotype}, genotypes{g});
                indexg=find(logicalg);
                nAve=sum(logicalg);
                tempmag=nan(20,8,nAve);
                tempphi=nan(20,8,nAve);
                tempgain1=nan(20,nAve);
                tempgain2=nan(20,nAve);
                tempQ=nan(8,nAve);
                tempBF=nan(8,nAve);
                for n=1:nAve
                    z=dataAA(c,indexg(n)).mag;
                    if length(z)>0
                        tempmag(:,:,n)=z;
                    end
                    z=dataAA(c,indexg(n)).phi;
                    if length(z)>0
                        tempphi(:,:,n)=z;
                    end
                    z=dataAA(c,indexg(n)).gain1;
                    if length(z)>0
                        tempgain1(:,n)=z;
                    end
                    z=dataAA(c,indexg(n)).gain2;
                    if length(z)>0
                        tempgain2(:,n)=z;
                    end
                    z=dataAA(c,indexg(n)).Q;
                    if length(z)>0
                        tempQ(:,n)=z;
                    end
                    z=dataAA(c,indexg(n)).BF;
                    if length(z)>0
                        tempBF(:,n)=z;
                    end
    
                    % now fill up a table to statistically compare awake vs anesth mice 
                    % for each genotype using a multiple linear regression, paired testing approach
                    % (linear mixed model in R) 
    
                    % mag/phase
                    temp=reshape(tempmag(:,:,n),[],1);
                    npts=length(temp);
                    v1Temp=table('Size',[npts 14], ...
                        'VariableTypes',{'int32','string','string','string','string','double', ...
                        'double','double','double','double', 'double', ...
                        'double','double','double'},'VariableNames', ...
                        {'experiment','condition','protocol','genotype','pupilSize', ...
                        'freq','level','mag','phi','gain1','gain2','q','bf','diameter'});
                    v1Temp.experiment=repmat(indexg(n),[npts 1]);
                    v1Temp.condition=repmat(conditions(c),[npts 1]);
                    v1Temp.protocol=repmat(protocols(p),[npts 1]);
                    v1Temp.genotype=repmat(genotypes(g),[npts 1]);
                    v1Temp.pupilSize=repmat('na',[npts 1]);
                    v1Temp.freq=reshape(repmat(F,[length(L) 1])',[],1);
                    v1Temp.level=reshape(repmat(L,[length(F) 1]),[],1);
                    v1Temp.mag=temp;
                    v1Temp.phi=reshape(tempphi(:,:,n),[],1);
    
                    % gain1/2
                    temp=reshape(tempgain1(:,n),[],1);
                    npts=length(temp);
                    v2Temp=table('Size',[npts 14], ...
                        'VariableTypes',{'int32','string','string','string','string','double', ...
                        'double','double','double','double', 'double', ...
                        'double','double','double'},'VariableNames', ...
                        {'experiment','condition','protocol','genotype','pupilSize', ...
                        'freq','level','mag','phi','gain1','gain2','q','bf','diameter'});
                    v2Temp.experiment=repmat(indexg(n),[npts 1]);
                    v2Temp.condition=repmat(conditions(c),[npts 1]);
                    v2Temp.protocol=repmat(protocols(p),[npts 1]);
                    v2Temp.genotype=repmat(genotypes(g),[npts 1]);
                    v2Temp.pupilSize=repmat('na',[npts 1]);
                    v2Temp.freq=F';
                    v2Temp.gain1=temp;
                    v2Temp.gain2=reshape(tempgain2(:,n),[],1);
    
                    % q, bf
                    temp=reshape(tempQ(:,n),[],1);
                    npts=length(temp);
                    v3Temp=table('Size',[npts 14], ...
                        'VariableTypes',{'int32','string','string','string','string','double', ...
                        'double','double','double','double', 'double', ...
                        'double','double','double'},'VariableNames', ...
                        {'experiment','condition','protocol','genotype','pupilSize', ...
                        'freq','level','mag','phi','gain1','gain2','q','bf','diameter'});
                    v3Temp.experiment=repmat(indexg(n),[npts 1]);
                    v3Temp.condition=repmat(conditions(c),[npts 1]);
                    v3Temp.protocol=repmat(protocols(p),[npts 1]);
                    v3Temp.genotype=repmat(genotypes(g),[npts 1]);
                    v3Temp.pupilSize=repmat('na',[npts 1]);
                    v3Temp.level=L';
                    v3Temp.q=temp;
                    v3Temp.bf=reshape(tempBF(:,n),[],1);
        
                    v1AA=[v1AA;v1Temp];
                    v2AA=[v2AA;v2Temp];
                    v3AA=[v3AA;v3Temp];
    
                end
                dataAAAve1(c,g).mag=mean(tempmag,3,'omitnan');
                dataAAAve1(c,g).magSD=std(tempmag,0,3,'omitnan');
                dataAAAve1(c,g).magN=size(tempmag,3) - sum(isnan(tempmag),3);
                dataAAAve1(c,g).phi=mean(tempphi,3,'omitnan');
                dataAAAve1(c,g).phiSD=std(tempphi,0,3,'omitnan');
                dataAAAve1(c,g).phiN=size(tempphi,3) - sum(isnan(tempphi),3);
                dataAAAve1(c,g).gain1=mean(tempgain1,2,'omitnan');
                dataAAAve1(c,g).gain1SD=std(tempgain1,0,2,'omitnan');
                dataAAAve1(c,g).gain1N=size(tempgain1,2) - sum(isnan(tempgain1),2);
                dataAAAve1(c,g).gain2=mean(tempgain2,2,'omitnan');
                dataAAAve1(c,g).gain2SD=std(tempgain2,0,2,'omitnan');
                dataAAAve1(c,g).gain2N=size(tempgain2,2) - sum(isnan(tempgain2),2);
                dataAAAve1(c,g).Q=mean(tempQ,2,'omitnan');
                dataAAAve1(c,g).QSD=std(tempQ,0,2,'omitnan');
                dataAAAve1(c,g).QN=size(tempQ,2) - sum(isnan(tempQ),2);
                dataAAAve1(c,g).BF=mean(tempBF,2,'omitnan');
                dataAAAve1(c,g).BFSD=std(tempBF,0,2,'omitnan');
                dataAAAve1(c,g).BFN=size(tempBF,2) - sum(isnan(tempBF),2);
            end
        end



    end
end

%% Plot summary figure about pupil sizes
c=1;
p=1;
f1=figure(1);
f1.Position=[20,70,750,1000];

t1=tiledlayout(8,6);
TOLC_base=viridis(7); % get distinct colors from a palette that works for colorblind
NF=0.8;     %noise floor
capsize=3;
xlims = [4 15];
ylims = [0.8 16];
ylimsP = [-3 .5];
x80=[4.3, 4.2];
y80=[9.1, 9.5];
x10=[11.0, 10.0];
y10=[0.8,0.8];

% Plot representative tuning curve
TOLC=[TOLC_base(2,:);TOLC_base(4,:);TOLC_base(6,:)];
e=3; % WT
nexttile(1,[2,2]);
x=data(c,p,e,1).Fmat/1000;
y=(data(c,p,e,1).mag);
y(y<NF)=nan;
yP1=(data(c,p,e,1).phi);
yP1(y<NF)=nan;
z=(data(c,p,e,1).magSD./sqrt(data(c,p,e,1).magN));
errorbar(x,y,z,'-', 'LineWidth', 1,'Color',TOLC(1,:),'CapSize',capsize, 'LineWidth', 2);
hold on;
y=(data(c,p,e,2).mag);
y(y<NF)=nan;
yP2=(data(c,p,e,2).phi);
yP2(y<NF)=nan;
z=(data(c,p,e,2).magSD./sqrt(data(c,p,e,2).magN));
errorbar(x,y,z,'-', 'LineWidth', 1,'Color',TOLC(2,:),'CapSize',capsize, 'LineWidth', 2);
y=(data(c,p,e,3).mag);
y(y<NF)=nan;
yP3=(data(c,p,e,3).phi);
yP3(y<NF)=nan;
z=(data(c,p,e,3).magSD./sqrt(data(c,p,e,3).magN));
errorbar(x,y,z,'-', 'LineWidth', 1,'Color',TOLC(3,:),'CapSize',capsize, 'LineWidth', 2);
xlim(xlims); set(gca,'Xscale','log','Xtick',[1 2 4 5 6 7 8 10 12 15 20 50]);
ylim(ylims); set(gca,'Yscale','log','Ytick',[0.1 0.2 0.5 1 2 5 10 20 50]);
xlabel('Freq (kHz)');ylabel('Magnitude (nm)');
title('WT-10');
text(x80(1),y80(1),'80','Color','k','FontSize',9)
% text(x10(1),y10(1),'10','Color','k','FontSize',9)
xtickangle(0);
hold off;

e=8; %alpha9
nexttile(13,[2,2]);
x=data(c,p,e,1).Fmat/1000;
y=(data(c,p,e,1).mag);
y(y<NF)=nan;
yP1=(data(c,p,e,1).phi);
yP1(y<NF)=nan;
z=(data(c,p,e,1).magSD./sqrt(data(c,p,e,1).magN));
errorbar(x,y,z,'-', 'LineWidth', 1,'Color',TOLC(1,:),'CapSize',capsize, 'LineWidth', 2)
hold on;
y=(data(c,p,e,2).mag);
y(y<NF)=nan;
yP2=(data(c,p,e,2).phi);
yP2(y<NF)=nan;
z=(data(c,p,e,2).magSD./sqrt(data(c,p,e,2).magN));
errorbar(x,y,z,'g-', 'LineWidth', 1,'Color',TOLC(2,:),'CapSize',capsize, 'LineWidth', 2)
y=(data(c,p,e,3).mag);
y(y<NF)=nan;
yP3=(data(c,p,e,3).phi);
yP3(y<NF)=nan;
z=(data(c,p,e,3).magSD./sqrt(data(c,p,e,3).magN));
errorbar(x,y,z,'r-', 'LineWidth', 1,'Color',TOLC(3,:),'CapSize',capsize, 'LineWidth', 2)
xlim(xlims); set(gca,'Xscale','log','Xtick',[1 2 4 5 6 7 8 10 12 15 20 50]);
ylim(ylims); set(gca,'Yscale','log','Ytick',[0.1 0.2 0.5 1 2 5 10 20 50]);
xlabel('Freq (kHz)');ylabel('Magnitude (nm)');
title('Alpha9^{-/-}-04');
text(x80(2),y80(2),'80','Color','k','FontSize',9)
% text(x10(2),y10(2),'10','Color','k','FontSize',9)
xtickangle(0);
hold off;

% Plot averaged tuning curves
g=1; % WT
nexttile(3,[2,2]);
x=data(c,p,e,1).Fmat/1000;
y=(dataAve1(c,p,g,1).mag);
y(y<NF)=nan;
yP1=(dataAve1(c,p,g,1).phi);
yP1(y<NF)=nan;
z=(dataAve1(c,p,g,1).magSD./sqrt(dataAve1(c,p,g,1).magN));
errorbar(x,y,z,'-', 'LineWidth', 1,'Color',TOLC(1,:),'CapSize',capsize, 'LineWidth', 2)
hold on;
y=(dataAve1(c,p,g,2).mag);
y(y<NF)=nan;
yP2=(dataAve1(c,p,g,2).phi);
yP2(y<NF)=nan;
z=(dataAve1(c,p,g,2).magSD./sqrt(dataAve1(c,p,g,2).magN));
errorbar(x,y,z,'-', 'LineWidth', 1,'Color',TOLC(2,:),'CapSize',capsize, 'LineWidth', 2)
y=(dataAve1(c,p,g,3).mag);
y(y<NF)=nan;
yP3=(dataAve1(c,p,g,3).phi);
yP3(y<NF)=nan;
z=(dataAve1(c,p,g,3).magSD./sqrt(dataAve1(c,p,g,3).magN));
errorbar(x,y,z,'-', 'LineWidth', 1,'Color',TOLC(3,:),'CapSize',capsize, 'LineWidth', 2)
xlim(xlims); set(gca,'Xscale','log','Xtick',[1 2 4 5 6 7 8 10 12 15 20 50]);
ylim(ylims); set(gca,'Yscale','log','Ytick',[0.1 0.2 0.5 1 2 5 10 20 50]);
xlabel('Freq (kHz)');ylabel('Magnitude (nm)');
G=genotypes{g};
if G=="Alpha9KO"
    G="Alpha9^{-/-}";
end
title(G);
text(11,13.5,'n=6','Color','k','FontSize',9)
text(x80(1),y80(1),'80','Color','k','FontSize',9)
xtickangle(0);
hold off;

nexttile(5,[2,2]);
z=(dataAve1(c,p,g,1).phiSD./sqrt(dataAve1(c,p,g,1).phiN));
errorbar(x,yP1,z,'-', 'LineWidth', 1,'Color',TOLC(1,:),'CapSize',capsize, 'LineWidth', 2);
hold on;
z=(dataAve1(c,p,g,2).phiSD./sqrt(dataAve1(c,p,g,2).phiN));
errorbar(x,yP2,z,'-', 'LineWidth', 1,'Color',TOLC(2,:),'CapSize',capsize, 'LineWidth', 2)
z=(dataAve1(c,p,g,3).phiSD./sqrt(dataAve1(c,p,g,3).phiN));
errorbar(x,yP3,z,'-', 'LineWidth', 1,'Color',TOLC(3,:),'CapSize',capsize, 'LineWidth', 2)
hold off;
xlim(xlims); set(gca,'Xscale','log','Xtick',[1 2 4 5 6 7 8 10 12 15 20 50]);
ylim(ylimsP); 
xlabel('Freq (kHz)');ylabel('Phase (cycles)');
text(8.3,0.2,'Small','Color',TOLC(1,:))
text(8.3,-0.3,'Medium','Color',TOLC(2,:))
text(8.3,-0.8,'Large','Color',TOLC(3,:))
G=genotypes{g};
if G=="Alpha9KO"
    G="Alpha9^{-/-}";
end
xtickangle(0);
title(G);

g=2; %alpha9
nexttile(15,[2,2]);
x=data(c,p,e,1).Fmat/1000;
y=(dataAve1(c,p,g,1).mag);
y(y<NF)=nan;
yP1=(dataAve1(c,p,g,1).phi);
yP1(y<NF)=nan;
z=(dataAve1(c,p,g,1).magSD./sqrt(dataAve1(c,p,g,1).magN));
errorbar(x,y,z,'-', 'LineWidth', 1,'Color',TOLC(1,:),'CapSize',capsize, 'LineWidth', 2)
hold on;
y=(dataAve1(c,p,g,2).mag);
y(y<NF)=nan;
yP2=(dataAve1(c,p,g,2).phi);
yP2(y<NF)=nan;
z=(dataAve1(c,p,g,2).magSD./sqrt(dataAve1(c,p,g,2).magN));
errorbar(x,y,z,'-', 'LineWidth', 1,'Color',TOLC(2,:),'CapSize',capsize, 'LineWidth', 2)
y=(dataAve1(c,p,g,3).mag);
y(y<NF)=nan;
yP3=(dataAve1(c,p,g,3).phi);
yP3(y<NF)=nan;
z=(dataAve1(c,p,g,3).magSD./sqrt(dataAve1(c,p,g,3).magN));
errorbar(x,y,z,'-', 'LineWidth', 1,'Color',TOLC(3,:),'CapSize',capsize, 'LineWidth', 2)
xlim(xlims); set(gca,'Xscale','log','Xtick',[1 2 4 5 6 7 8 10 12 15 20 50]);
ylim(ylims); set(gca,'Yscale','log','Ytick',[0.1 0.2 0.5 1 2 5 10 20 50]);
xlabel('Freq (kHz)');ylabel('Magnitude (nm)');
G=genotypes{g};
if G=="Alpha9KO"
    G="Alpha9^{-/-}";
end
title(G);
text(11,13.5,'n=6','Color','k','FontSize',9)
text(x80(2),y80(2),'80','Color','k','FontSize',9)
xtickangle(0);
hold off;

nexttile(17,[2,2]);
z=(dataAve1(c,p,g,1).phiSD./sqrt(dataAve1(c,p,g,1).phiN));
errorbar(x,yP1,z,'-', 'LineWidth', 1,'Color',TOLC(1,:),'CapSize',capsize, 'LineWidth', 2);
hold on;
z=(dataAve1(c,p,g,2).phiSD./sqrt(dataAve1(c,p,g,2).phiN));
errorbar(x,yP2,z,'-', 'LineWidth', 1,'Color',TOLC(2,:),'CapSize',capsize, 'LineWidth', 2)
z=(dataAve1(c,p,g,3).phiSD./sqrt(dataAve1(c,p,g,3).phiN));
errorbar(x,yP3,z,'-', 'LineWidth', 1,'Color',TOLC(3,:),'CapSize',capsize, 'LineWidth', 2)
hold off;
xlim(xlims); set(gca,'Xscale','log','Xtick',[1 2 4 5 6 7 8 10 12 15 20 50]);
ylim(ylimsP); 
xlabel('Freq (kHz)');ylabel('Phase (cycles)');
G=genotypes{g};
if G=="Alpha9KO"
    G="Alpha9^{-/-}";
end
xtickangle(0);
title(G);

% Plot level dependence of BF/Q10dB
nexttile(25,[2,2]);
g=1; % WT
x=F/1000;
y=(dataAve1(c,p,g,1).gain1);
z=(dataAve1(c,p,g,1).gain1SD./sqrt(dataAve1(c,p,g,1).gain1N));
errorbar(x,y,z,'o-', 'LineWidth', 1,'Color',TOLC(1,:),'CapSize',capsize, 'LineWidth', 2)
hold on;
y=(dataAve1(c,p,g,2).gain1);
z=(dataAve1(c,p,g,2).gain1SD./sqrt(dataAve1(c,p,g,2).gain1N));
errorbar(x,y,z,'o-', 'LineWidth', 1,'Color',TOLC(2,:),'CapSize',capsize, 'LineWidth', 2)
y=(dataAve1(c,p,g,3).gain1);
z=(dataAve1(c,p,g,3).gain1SD./sqrt(dataAve1(c,p,g,3).gain1N));
errorbar(x,y,z,'o-', 'LineWidth', 1,'Color',TOLC(3,:),'CapSize',capsize, 'LineWidth', 2)
xlim([8 13]); set(gca,'Xscale','log','Xtick',[1 2 4 5 6 7 8 9 10 11 12 13 14 15 20 50]);
ylim([10 40]); set(gca,'Yscale','linear','Ytick',[0 10 20 30 40 50]);
xlabel('Freq (kHz)');
ylabel('Gain (dB)');
text(8.2,37,'20-80 dB SPL','Color','k','FontSize',9)
G=genotypes{g};
if G=="Alpha9KO"
    G="Alpha9^{-/-}";
end
title(G);
hold off;

nexttile(37,[2,2]);
g=2; % Alpha9
x=F/1000;
y=(dataAve1(c,p,g,1).gain1);
z=(dataAve1(c,p,g,1).gain1SD./sqrt(dataAve1(c,p,g,1).gain1N));
errorbar(x,y,z,'o-', 'LineWidth', 1,'Color',TOLC(1,:),'CapSize',capsize, 'LineWidth', 2)
hold on;
y=(dataAve1(c,p,g,2).gain1);
z=(dataAve1(c,p,g,2).gain1SD./sqrt(dataAve1(c,p,g,2).gain1N));
errorbar(x,y,z,'o-', 'LineWidth', 1,'Color',TOLC(2,:),'CapSize',capsize, 'LineWidth', 2)
y=(dataAve1(c,p,g,3).gain1);
z=(dataAve1(c,p,g,3).gain1SD./sqrt(dataAve1(c,p,g,3).gain1N));
errorbar(x,y,z,'o-', 'LineWidth', 1,'Color',TOLC(3,:),'CapSize',capsize, 'LineWidth', 2)
xlim([8 13]); set(gca,'Xscale','log','Xtick',[1 2 4 5 6 7 8 9 10 11 12 13 14 15 20 50]);
ylim([10 40]); set(gca,'Yscale','linear','Ytick',[0 10 20 30 40 50]);
xlabel('Freq (kHz)');
ylabel('Gain (dB)');
text(8.2,37,'20-80 dB SPL','Color','k','FontSize',9)
G=genotypes{g};
if G=="Alpha9KO"
    G="Alpha9^{-/-}";
end
title(G);
hold off;

nexttile(27,[2,2]);
g=1; % WT
x=L;
y=(dataAve1(c,p,g,1).BF)/1000;
z=(dataAve1(c,p,g,1).BFSD./sqrt(dataAve1(c,p,g,1).BFN))/1000;
errorbar(x,y,z,'o-', 'LineWidth', 1,'Color',TOLC(1,:),'CapSize',capsize, 'LineWidth', 2)
hold on;
y=(dataAve1(c,p,g,2).BF)/1000;
z=(dataAve1(c,p,g,2).BFSD./sqrt(dataAve1(c,p,g,2).BFN))/1000;
errorbar(x,y,z,'o-', 'LineWidth', 1,'Color',TOLC(2,:),'CapSize',capsize, 'LineWidth', 2)
y=(dataAve1(c,p,g,3).BF)/1000;
z=(dataAve1(c,p,g,3).BFSD./sqrt(dataAve1(c,p,g,3).BFN))/1000;
errorbar(x,y,z,'o-', 'LineWidth', 1,'Color',TOLC(3,:),'CapSize',capsize, 'LineWidth', 2)
xlim([10 80]); set(gca,'Xscale','linear','Xtick',[10 30 50 70]);
ylim([8 12]); set(gca,'Yscale','linear','Ytick',[6 8 9 10 11 12 14 16]);
xlabel('Level (dB SPL)');ylabel('BF (kHz)');
G=genotypes{g};
if G=="Alpha9KO"
    G="Alpha9^{-/-}";
end
title(G);
hold off;

nexttile(39,[2,2]);
g=2; % alpha9
x=L;
y=(dataAve1(c,p,g,1).BF)/1000;
z=(dataAve1(c,p,g,1).BFSD./sqrt(dataAve1(c,p,g,1).BFN))/1000;
errorbar(x,y,z,'o-', 'LineWidth', 1,'Color',TOLC(1,:),'CapSize',capsize, 'LineWidth', 2)
hold on;
y=(dataAve1(c,p,g,2).BF)/1000;
z=(dataAve1(c,p,g,2).BFSD./sqrt(dataAve1(c,p,g,2).BFN))/1000;
errorbar(x,y,z,'o-', 'LineWidth', 1,'Color',TOLC(2,:),'CapSize',capsize, 'LineWidth', 2)
y=(dataAve1(c,p,g,3).BF)/1000;
z=(dataAve1(c,p,g,3).BFSD./sqrt(dataAve1(c,p,g,3).BFN))/1000;
errorbar(x,y,z,'o-', 'LineWidth', 1,'Color',TOLC(3,:),'CapSize',capsize, 'LineWidth', 2)
xlim([10 80]); set(gca,'Xscale','linear','Xtick',[10 30 50 70]);
ylim([8 12]); set(gca,'Yscale','linear','Ytick',[6 8 9 10 11 12 14 16]);
xlabel('Level (dB SPL)');ylabel('BF (kHz)');
G=genotypes{g};
if G=="Alpha9KO"
    G="Alpha9^{-/-}";
end
title(G);
hold off;

nexttile(29,[2,2]);
g=1; % WT
x=L;
y=(dataAve1(c,p,g,1).Q);
z=(dataAve1(c,p,g,1).QSD./sqrt(dataAve1(c,p,g,1).QN));
errorbar(x,y,z,'o-', 'LineWidth', 1,'Color',TOLC(1,:),'CapSize',capsize, 'LineWidth', 2)
hold on;
y=(dataAve1(c,p,g,2).Q);
z=(dataAve1(c,p,g,2).QSD./sqrt(dataAve1(c,p,g,2).QN));
errorbar(x,y,z,'o-', 'LineWidth', 1,'Color',TOLC(2,:),'CapSize',capsize, 'LineWidth', 2)
y=(dataAve1(c,p,g,3).Q);
z=(dataAve1(c,p,g,3).QSD./sqrt(dataAve1(c,p,g,3).QN));
errorbar(x,y,z,'o-', 'LineWidth', 1,'Color',TOLC(3,:),'CapSize',capsize, 'LineWidth', 2)
xlim([10 80]); set(gca,'Xscale','linear','Xtick',[10 30 50 70]);
ylim([1 6]); set(gca,'Yscale','linear','Ytick',[0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 16]);
xlabel('Level (dB SPL)');ylabel('Q10dB');
G=genotypes{g};
if G=="Alpha9KO"
    G="Alpha9^{-/-}";
end
title(G);
hold off;

nexttile(41,[2,2]);
g=2; % Alpha9
x=L;
y=(dataAve1(c,p,g,1).Q);
z=(dataAve1(c,p,g,1).QSD./sqrt(dataAve1(c,p,g,1).QN));
errorbar(x,y,z,'o-', 'LineWidth', 1,'Color',TOLC(1,:),'CapSize',capsize, 'LineWidth', 2)
hold on;
y=(dataAve1(c,p,g,2).Q);
z=(dataAve1(c,p,g,2).QSD./sqrt(dataAve1(c,p,g,2).QN));
errorbar(x,y,z,'o-', 'LineWidth', 1,'Color',TOLC(2,:),'CapSize',capsize, 'LineWidth', 2)
y=(dataAve1(c,p,g,3).Q);
z=(dataAve1(c,p,g,3).QSD./sqrt(dataAve1(c,p,g,3).QN));
errorbar(x,y,z,'o-', 'LineWidth', 1,'Color',TOLC(3,:),'CapSize',capsize, 'LineWidth', 2)
xlim([10 80]); set(gca,'Xscale','linear','Xtick',[10 30 50 70]);
ylim([1 6]); set(gca,'Yscale','linear','Ytick',[0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 16]);
xlabel('Level (dB SPL)');ylabel('Q10dB');
G=genotypes{g};
if G=="Alpha9KO"
    G="Alpha9^{-/-}";
end
title(G);
hold off;

fontsize(f1,scale=1.5);
t1.TileSpacing = 'compact';
t1.Padding = 'compact';
shg
%% Plot data binned by % pupil diameter

f2=figure(2);
f2.Position=[800,70,750,500];
t2=tiledlayout(4,6);
TOLC=[TOLC_base(1,:);TOLC_base(5,:)];

% plot binned mag/gain/BF/Q data
xBin=dataScatter(c,p,1).xBin;
yMagAve1=mean(yMagAveData(1:length(dataAll)/2,:),'omitnan');
yMagAve2=mean(yMagAveData((length(dataAll)/2+1):end,:),'omitnan');
yMagSEM1=std(yMagAveData(1:length(dataAll)/2,:),'omitnan')/sqrt(length(dataAll)/2);
yMagSEM2=std(yMagAveData((length(dataAll)/2+1):end,:),'omitnan')/sqrt(length(dataAll)/2);
yMagLowAve1=mean(yMagLowAveData(1:length(dataAll)/2,:),'omitnan');
yMagLowAve2=mean(yMagLowAveData((length(dataAll)/2+1):end,:),'omitnan');
yMagLowSEM1=std(yMagLowAveData(1:length(dataAll)/2,:),'omitnan')/sqrt(length(dataAll)/2);
yMagLowSEM2=std(yMagLowAveData((length(dataAll)/2+1):end,:),'omitnan')/sqrt(length(dataAll)/2);
yPhaseAve1=mean(yPhaseAveData(1:length(dataAll)/2,:),'omitnan');
yPhaseAve2=mean(yPhaseAveData((length(dataAll)/2+1):end,:),'omitnan');
yPhaseSEM1=std(yPhaseAveData(1:length(dataAll)/2,:),'omitnan')/sqrt(length(dataAll)/2);
yPhaseSEM2=std(yPhaseAveData((length(dataAll)/2+1):end,:),'omitnan')/sqrt(length(dataAll)/2);
yGainAve1=mean(yGainAveData(1:length(dataAll)/2,:),'omitnan');
yGainAve2=mean(yGainAveData((length(dataAll)/2+1):end,:),'omitnan');
yGainSEM1=std(yGainAveData(1:length(dataAll)/2,:),'omitnan')/sqrt(length(dataAll)/2);
yGainSEM2=std(yGainAveData((length(dataAll)/2+1):end,:),'omitnan')/sqrt(length(dataAll)/2);
yBFAve1=mean(yBFAveData(1:length(dataAll)/2,:),'omitnan')/1000;
yBFAve2=mean(yBFAveData((length(dataAll)/2+1):end,:),'omitnan')/1000;
yBFSEM1=std(yBFAveData(1:length(dataAll)/2,:),'omitnan')/sqrt(length(dataAll)/2)/1000;
yBFSEM2=std(yBFAveData((length(dataAll)/2+1):end,:),'omitnan')/sqrt(length(dataAll)/2)/1000;
yQAve1=mean(yQAveData(1:length(dataAll)/2,:),'omitnan');
yQAve2=mean(yQAveData((length(dataAll)/2+1):end,:),'omitnan');
yQSEM1=std(yQAveData(1:length(dataAll)/2,:),'omitnan')/sqrt(length(dataAll)/2);
yQSEM2=std(yQAveData((length(dataAll)/2+1):end,:),'omitnan')/sqrt(length(dataAll)/2);

nexttile(1,[2,2]);
errorbar(xBin,yGainAve1,yGainSEM1,'o-','Color',TOLC(1,:),'CapSize',capsize, 'LineWidth', 2);
hold on;
errorbar(xBin,yGainAve2,yGainSEM2,'o-','Color',TOLC(2,:),'CapSize',capsize, 'LineWidth', 2);
hold off;
text(2,23.8,'WT (n=6)','Color',TOLC(1,:))
text(2,22.2,'Alpha9^{-/-} (n=6)','Color',TOLC(2,:))
xlim([0 xBin(end)]);
ylim([21 37]);
xticks([0 25 50 75 100]);
xtickangle(0);
xlabel('Pupil (% max)');
ylabel('Gain (dB)');
text(2,36,'20-80 dB SPL','Color','k','FontSize',9)

nexttile(3,[2,2]);
errorbar(xBin,yBFAve1,yBFSEM1,'o-','Color',TOLC(1,:),'CapSize',capsize, 'LineWidth', 2);
hold on;
errorbar(xBin,yBFAve2,yBFSEM2,'o-','Color',TOLC(2,:),'CapSize',capsize, 'LineWidth', 2);
ylim([9.5 12]);
hold off;
xlim([0 xBin(end)]);
xticks([0 25 50 75 100]);
xtickangle(0);
xlabel('Pupil (% max)');
ylabel('CF (kHz)');

nexttile(5,[2,2]);
errorbar(xBin,yQAve1,yQSEM1,'o-','Color',TOLC(1,:),'CapSize',capsize, 'LineWidth', 2);
hold on;
errorbar(xBin,yQAve2,yQSEM2,'o-','Color',TOLC(2,:),'CapSize',capsize, 'LineWidth', 2);
ylim([3.5 5.0]);
hold off;
xlim([0 xBin(end)]);
xticks([0 25 50 75 100]);
xtickangle(0);
xlabel('Pupil (% max)');
ylabel('Q10dB');

nexttile(13,[2,2]);
errorbar(xBin,yMagAve1,yMagSEM1,'o-','Color',TOLC(1,:),'CapSize',capsize, 'LineWidth', 2);
hold on;
errorbar(xBin,yMagAve2,yMagSEM2,'o-','Color',TOLC(2,:),'CapSize',capsize, 'LineWidth', 2);
hold off;
xlim([0 xBin(end)]);
ylim([1.5 3.5]);
xticks([0 25 50 75 100]);
xtickangle(0);
xlabel('Pupil (% max)');
ylabel('Mag at BF(nm)');
text(2,3.3,'20 dB SPL','Color','k','FontSize',9)

nexttile(15,[2,2]);
errorbar(xBin,yMagLowAve1,yMagLowSEM1,'o-','Color',TOLC(1,:),'CapSize',capsize, 'LineWidth', 2);
hold on;
errorbar(xBin,yMagLowAve2,yMagLowSEM2,'o-','Color',TOLC(2,:),'CapSize',capsize, 'LineWidth', 2);
hold off;
xlim([0 xBin(end)]);
xticks([0 25 50 75 100]);
xtickangle(0);
ylim([1.5 3.5]);
xlabel('Pupil (% max)');
ylabel('Mag at 0.5*BF (nm)');
text(2,3.3,'50 dB SPL','Color','k','FontSize',9)

nexttile(17,[2,2]);
errorbar(xBin,yPhaseAve1,yPhaseSEM1,'o-','Color',TOLC(1,:),'CapSize',capsize, 'LineWidth', 2);
hold on;
errorbar(xBin,yPhaseAve2,yPhaseSEM2,'o-','Color',TOLC(2,:),'CapSize',capsize, 'LineWidth', 2);
hold off;
xlim([0 xBin(end)]);
ylim([-3.5 -1.5]);
xticks([0 25 50 75 100]);
xtickangle(0);
xlabel('Pupil (% max)');
ylabel('Phase at CF (cycles)');

fontsize(gcf,scale=1.5);
t2.TileSpacing = 'compact';
t2.Padding = 'compact';
shg


%% Plot summary figure about awake vs anesth
f3=figure(3);
f3.Position=[2200,70,750, 1000];
t3=tiledlayout(8,6);
TOLC=[TOLC_base(1,:);TOLC_base(5,:)];

% Plot representative tuning curve
e=3; % WT
c=1; % awake
nexttile(1,[2,2]);
x=dataAA(c,e).Fmat/1000;
y=(dataAA(c,e).mag);
y(y<NF)=nan;
yP1=(dataAA(c,e).phi);
yP1(y<NF)=nan;
z=(dataAA(c,e).magSD./sqrt(dataAA(c,e).magN));
errorbar(x,y,z,'-', 'LineWidth', 1,'Color',TOLC(1,:),'CapSize',capsize, 'LineWidth', 2);
hold on;
c=2; % anesth
y=(dataAA(c,e).mag);
y(y<NF)=nan;
yP1=(dataAA(c,e).phi);
yP2(y<NF)=nan;
z=(dataAA(c,e).magSD./sqrt(dataAA(c,e).magN));
errorbar(x,y,z,':', 'LineWidth', 1,'Color',TOLC(1,:),'CapSize',capsize, 'LineWidth', 2);
xlim(xlims); set(gca,'Xscale','log','Xtick',[1 2 4 5 6 7 8 10 12 15 20 50]);
ylim(ylims); set(gca,'Yscale','log','Ytick',[0.1 0.2 0.5 1 2 5 10 20 50]);
xlabel('Freq (kHz)');ylabel('Magnitude (nm)');
title('WT-10');
text(x80(1),y80(1),'80','Color','k','FontSize',9)
% text(x10(1),y10(1),'10','Color','k','FontSize',9);
xtickangle(0);
hold off;

e=8; %alpha9
c=1; % awake
nexttile(13,[2,2]);
x=dataAA(c,e).Fmat/1000;
y=(dataAA(c,e).mag);
y(y<NF)=nan;
yP1=(dataAA(c,e).phi);
yP1(y<NF)=nan;
z=(dataAA(c,e).magSD./sqrt(dataAA(c,e).magN));
errorbar(x,y,z,'-', 'LineWidth', 1,'Color',TOLC(2,:),'CapSize',capsize, 'LineWidth', 2);
hold on;
c=2; % anesth
y=(dataAA(c,e).mag);
y(y<NF)=nan;
yP1=(dataAA(c,e).phi);
yP2(y<NF)=nan;
z=(dataAA(c,e).magSD./sqrt(dataAA(c,e).magN));
errorbar(x,y,z,':', 'LineWidth', 1,'Color',TOLC(2,:),'CapSize',capsize, 'LineWidth', 2);
xlim(xlims); set(gca,'Xscale','log','Xtick',[1 2 4 5 6 7 8 10 12 15 20 50]);
ylim(ylims); set(gca,'Yscale','log','Ytick',[0.1 0.2 0.5 1 2 5 10 20 50]);
xlabel('Freq (kHz)');ylabel('Magnitude (nm)');
title('Alpha9^{-/-}-04');
text(x80(2),y80(2),'80','Color','k','FontSize',9)
% text(x10(2),y10(2),'10','Color','k','FontSize',9)
xtickangle(0);
hold off;

% Plot averaged tuning curves
g=1; % WT
nexttile(3,[2,2]);
x=dataAA(c,e).Fmat/1000;
c=1;
y=(dataAAAve1(c,g).mag);
y(y<NF)=nan;
yP1=(dataAAAve1(c,g).phi);
yP1(y<NF)=nan;
z=(dataAAAve1(c,g).magSD./sqrt(dataAAAve1(c,g).magN));
errorbar(x,y,z,'-', 'LineWidth', 1,'Color',TOLC(1,:),'CapSize',capsize, 'LineWidth', 2)
hold on;
c=2;
y=(dataAAAve1(c,g).mag);
y(y<NF)=nan;
yP1=(dataAAAve1(c,g).phi);
yP1(y<NF)=nan;
z=(dataAAAve1(c,g).magSD./sqrt(dataAAAve1(c,g).magN));
errorbar(x,y,z,':', 'LineWidth', 1,'Color',TOLC(1,:),'CapSize',capsize, 'LineWidth', 2)
xlim(xlims); set(gca,'Xscale','log','Xtick',[1 2 4 5 6 7 8 10 12 15 20 50]);
ylim(ylims); set(gca,'Yscale','log','Ytick',[0.1 0.2 0.5 1 2 5 10 20 50]);
xlabel('Freq (kHz)');ylabel('Magnitude (nm)');
G=genotypes{g};
if G=="Alpha9KO"
    G="Alpha9^{-/-}";
end
title(G);
text(11,13.5,'n=6','Color','k','FontSize',9)
text(x80(1),y80(1),'80','Color','k','FontSize',9);

xtickangle(0);
hold off;

nexttile(5,[2,2]);
c=1;
z=(dataAAAve1(c,g).phiSD./sqrt(dataAAAve1(c,g).phiN));
errorbar(x,yP1,z,'-', 'LineWidth', 1,'Color',TOLC(1,:),'CapSize',capsize, 'LineWidth', 2);
hold on;
c=2;
z=(dataAAAve1(c,g).phiSD./sqrt(dataAAAve1(c,g).phiN));
errorbar(x,yP2,z,':', 'LineWidth', 1,'Color',TOLC(1,:),'CapSize',capsize, 'LineWidth', 2)
hold off;
xlim(xlims); set(gca,'Xscale','log','Xtick',[1 2 4 5 6 7 8 10 12 15 20 50]);
ylim(ylimsP); 
xlabel('Freq (kHz)');ylabel('Phase (cycles)');
G=genotypes{g};
if G=="Alpha9KO"
    G="Alpha9^{-/-}";
end
title(G);
text(4.3,0.3,'Awake','Color',TOLC(1,:),'FontSize',9)
text(4.2,0.0,'\it Anesthetized','Color',TOLC(1,:),'FontSize',9)
xtickangle(0);

g=2; %alpha9
nexttile(15,[2,2]);
x=dataAA(c,e).Fmat/1000;
c=1;
y=(dataAAAve1(c,g).mag);
y(y<NF)=nan;
yP1=(dataAAAve1(c,g).phi);
yP1(y<NF)=nan;
z=(dataAAAve1(c,g).magSD./sqrt(dataAAAve1(c,g).magN));
errorbar(x,y,z,'-', 'LineWidth', 1,'Color',TOLC(2,:),'CapSize',capsize, 'LineWidth', 2)
hold on;
c=2;
y=(dataAAAve1(c,g).mag);
y(y<NF)=nan;
yP1=(dataAAAve1(c,g).phi);
yP1(y<NF)=nan;
z=(dataAAAve1(c,g).magSD./sqrt(dataAAAve1(c,g).magN));
errorbar(x,y,z,':', 'LineWidth', 1,'Color',TOLC(2,:),'CapSize',capsize, 'LineWidth', 2)
xlim(xlims); set(gca,'Xscale','log','Xtick',[1 2 4 5 6 7 8 10 12 15 20 50]);
ylim(ylims); set(gca,'Yscale','log','Ytick',[0.1 0.2 0.5 1 2 5 10 20 50]);
xlabel('Freq (kHz)');ylabel('Magnitude (nm)');
G=genotypes{g};
if G=="Alpha9KO"
    G="Alpha9^{-/-}";
end
title(G);
text(11,13.5,'n=6','Color','k','FontSize',9)
text(x80(2),y80(2),'80','Color','k','FontSize',9)
xtickangle(0);
hold off;

nexttile(17,[2,2]);
c=1;
z=(dataAAAve1(c,g).phiSD./sqrt(dataAAAve1(c,g).phiN));
errorbar(x,yP1,z,'-', 'LineWidth', 1,'Color',TOLC(2,:),'CapSize',capsize, 'LineWidth', 2);
hold on;
c=2;
z=(dataAAAve1(c,g).phiSD./sqrt(dataAAAve1(c,g).phiN));
errorbar(x,yP2,z,':', 'LineWidth', 1,'Color',TOLC(2,:),'CapSize',capsize, 'LineWidth', 2)
hold off;
xlim(xlims); set(gca,'Xscale','log','Xtick',[1 2 4 5 6 7 8 10 12 15 20 50]);
ylim(ylimsP); 
xlabel('Freq (kHz)');ylabel('Phase (cycles)');
G=genotypes{g};
if G=="Alpha9KO"
    G="Alpha9^{-/-}";
end
title(G);
text(4.3,0.3,'Awake','Color',TOLC(2,:),'FontSize',9)
text(4.2,0.0,'\it Anesthetized','Color',TOLC(2,:),'FontSize',9)
xtickangle(0);

% Plot level dependence of BF/Q10dB
nexttile(25,[2,2]);
g=1; % WT
x=F/1000;
c=1;
y=(dataAAAve1(c,g).gain1);
z=(dataAAAve1(c,g).gain1SD./sqrt(dataAAAve1(c,g).gain1N));
errorbar(x,y,z,'o-', 'LineWidth', 1,'Color',TOLC(1,:),'CapSize',capsize, 'LineWidth', 2)
hold on;
c=2;
y=(dataAAAve1(c,g).gain1);
z=(dataAAAve1(c,g).gain1SD./sqrt(dataAAAve1(c,g).gain1N));
errorbar(x,y,z,'o:', 'LineWidth', 1,'Color',TOLC(1,:),'CapSize',capsize, 'LineWidth', 2)
xlim([8 13]); set(gca,'Xscale','log','Xtick',[1 2 4 5 6 7 8 9 10 11 12 13 14 15 20 50]);
ylim([10 40]); set(gca,'Yscale','linear','Ytick',[0 10 20 30 40 50]);
xlabel('Freq (kHz)');
ylabel('Gain (dB)');
text(8.2,37,'20-80 dB SPL','Color','k','FontSize',9)
G=genotypes{g};
if G=="Alpha9KO"
    G="Alpha9^{-/-}";
end
title(G);
hold off;
xtickangle(0);

nexttile(37,[2,2]);
g=2; % Alpha9
x=F/1000;
c=1;
y=(dataAAAve1(c,g).gain1);
z=(dataAAAve1(c,g).gain1SD./sqrt(dataAAAve1(c,g).gain1N));
errorbar(x,y,z,'o-', 'LineWidth', 1,'Color',TOLC(2,:),'CapSize',capsize, 'LineWidth', 2)
hold on;
c=2;
y=(dataAAAve1(c,g).gain1);
z=(dataAAAve1(c,g).gain1SD./sqrt(dataAAAve1(c,g).gain1N));
errorbar(x,y,z,'o:', 'LineWidth', 1,'Color',TOLC(2,:),'CapSize',capsize, 'LineWidth', 2)
xlim([8 13]); set(gca,'Xscale','log','Xtick',[1 2 4 5 6 7 8 9 10 11 12 13 14 15 20 50]);
ylim([10 40]); set(gca,'Yscale','linear','Ytick',[0 10 20 30 40 50]);
xlabel('Freq (kHz)');
ylabel('Gain (dB)');
text(8.2,37,'20-80 dB SPL','Color','k','FontSize',9)
G=genotypes{g};
if G=="Alpha9KO"
    G="Alpha9^{-/-}";
end
title(G);
hold off;
xtickangle(0);

nexttile(27,[2,2]);
g=1; % WT
x=L;
c=1;
y=(dataAAAve1(c,g).BF)/1000;
z=(dataAAAve1(c,g).BFSD./sqrt(dataAAAve1(c,g).BFN))/1000;
errorbar(x,y,z,'o-', 'LineWidth', 1,'Color',TOLC(1,:),'CapSize',capsize, 'LineWidth', 2)
hold on;
c=2;
y=(dataAAAve1(c,g).BF)/1000;
z=(dataAAAve1(c,g).BFSD./sqrt(dataAAAve1(c,g).BFN))/1000;
errorbar(x,y,z,'o:', 'LineWidth', 1,'Color',TOLC(1,:),'CapSize',capsize, 'LineWidth', 2)
xlim([10 80]); set(gca,'Xscale','linear','Xtick',[10 30 50 70]);
ylim([8 12]); set(gca,'Yscale','linear','Ytick',[6 8 9 10 11 12 14 16]);
xlabel('Level (dB SPL)');ylabel('BF (kHz)');
G=genotypes{g};
title(G);
hold off;

nexttile(39,[2,2]);
g=2; % alpha9
x=L;
c=1;
y=(dataAAAve1(c,g).BF)/1000;
z=(dataAAAve1(c,g).BFSD./sqrt(dataAAAve1(c,g).BFN))/1000;
errorbar(x,y,z,'o-', 'LineWidth', 1,'Color',TOLC(2,:),'CapSize',capsize, 'LineWidth', 2)
hold on;
c=2;
y=(dataAAAve1(c,g).BF)/1000;
z=(dataAAAve1(c,g).BFSD./sqrt(dataAAAve1(c,g).BFN))/1000;
errorbar(x,y,z,'o:', 'LineWidth', 1,'Color',TOLC(2,:),'CapSize',capsize, 'LineWidth', 2)
xlim([10 80]); set(gca,'Xscale','linear','Xtick',[10 30 50 70]);
ylim([8 12]); set(gca,'Yscale','linear','Ytick',[6 8 9 10 11 12 14 16]);
xlabel('Level (dB SPL)');ylabel('BF (kHz)');
G=genotypes{g};
if G=="Alpha9KO"
    G="Alpha9^{-/-}";
end
title(G);
hold off;

nexttile(29,[2,2]);
g=1; % WT
x=L;
c=1;
y=(dataAAAve1(c,g).Q);
z=(dataAAAve1(c,g).QSD./sqrt(dataAAAve1(c,g).QN));
errorbar(x,y,z,'o-', 'LineWidth', 1,'Color',TOLC(1,:),'CapSize',capsize, 'LineWidth', 2)
hold on;
c=2;
y=(dataAAAve1(c,g).Q);
z=(dataAAAve1(c,g).QSD./sqrt(dataAAAve1(c,g).QN));
errorbar(x,y,z,'o:', 'LineWidth', 1,'Color',TOLC(1,:),'CapSize',capsize, 'LineWidth', 2)
xlim([10 80]); set(gca,'Xscale','linear','Xtick',[10 30 50 70]);
ylim([1 6]); set(gca,'Yscale','linear','Ytick',[0 1 2 3 4 5 6 8 10 12 14 16]);
xlabel('Level (dB SPL)');ylabel('Q10dB');
G=genotypes{g};
if G=="Alpha9KO"
    G="Alpha9^{-/-}";
end
title(G);
hold off;

nexttile(41,[2,2]);
g=2; % Alpha9
x=L;
c=1;
y=(dataAAAve1(c,g).Q);
z=(dataAAAve1(c,g).QSD./sqrt(dataAAAve1(c,g).QN));
errorbar(x,y,z,'o-', 'LineWidth', 1,'Color',TOLC(2,:),'CapSize',capsize, 'LineWidth', 2)
hold on;
c=2;
y=(dataAAAve1(c,g).Q);
z=(dataAAAve1(c,g).QSD./sqrt(dataAAAve1(c,g).QN));
errorbar(x,y,z,'o:', 'LineWidth', 1,'Color',TOLC(2,:),'CapSize',capsize, 'LineWidth', 2)
xlim([10 80]); set(gca,'Xscale','linear','Xtick',[10 30 50 70]);
ylim([1 6]); set(gca,'Yscale','linear','Ytick',[0 1 2 3 4 5 6 8 10 12 14 16]);
xlabel('Level (dB SPL)');ylabel('Q10dB');
G=genotypes{g};
if G=="Alpha9KO"
    G="Alpha9^{-/-}";
end
title(G);
hold off;

fontsize(f3,scale=1.5);
t3.TileSpacing = 'compact';
t3.Padding = 'compact';
shg


%% Save the data for R analyses

% save data organized by small, med, large pupils
v1(1,:)=[];
v2(1,:)=[];
v3(1,:)=[];
v4(1,:)=[];
v1=removevars(v1,{'gain1','gain2','q','bf','diameter'});
v2=removevars(v2,{'level','mag','phi','q','bf','diameter'});
v3=removevars(v3,{'freq','mag','phi','gain1','gain2','diameter'});
v4=removevars(v4,{'freq','level','mag','phi','gain1','gain2','q','bf'});
writetable(v1,'v1');
writetable(v2,'v2');
writetable(v3,'v3');
writetable(v4,'v4');

% save data organized by awake/anesth
v1AA(1,:)=[];
v2AA(1,:)=[];
v3AA(1,:)=[];
v1AA=removevars(v1AA,{'gain1','gain2','q','bf','diameter'});
v2AA=removevars(v2AA,{'level','mag','phi','q','bf','diameter'});
v3AA=removevars(v3AA,{'freq','mag','phi','gain1','gain2','diameter'});
writetable(v1AA,'v1AA');
writetable(v2AA,'v2AA');
writetable(v3AA,'v3AA');

% save data binned by pupil diameter for stats
% save('yMagAveData.txt','yMagAveData','-ascii');
% save('yGainAveData.txt','yGainAveData','-ascii');
% save('yBFAveData.txt','yBFAveData','-ascii');
% save('yQAveData.txt','yQAveData','-ascii');

%% Put all binned data into one table for R stats analysis and save it
vBin=table('Size',[nBins*length(dataAll) 9], ...
    'VariableTypes',{'int32','string','int32', ...
    'double','double','double','double', 'double', ...
    'double'}, 'VariableNames', ...
    {'experiment','genotype','pupilBin', ...
    'gain','cf','q','mag','maglow','phase'});
vBin.mag=reshape(yMagAveData',[],1);
vBin.maglow=reshape(yMagLowAveData',[],1);
vBin.phase=reshape(yPhaseAveData',[],1);
vBin.gain=reshape(yGainAveData',[],1);
vBin.cf=reshape(yBFAveData',[],1);
vBin.q=reshape(yQAveData',[],1);
vBin.pupilBin=repmat((1:nBins)',[length(dataAll) 1]);
bin=1;
for e=1:length(dataAll)
    for b=1:nBins
        vBin.genotype(bin)=dataAll(e).genotype;
        vBin.experiment(bin)=dataAll(e).experiment;
        bin=bin+1;
    end
end
writetable(vBin,'vBin');

%%
function pupilSizeAll=calc_pupilSizes(dataAll,percent1,percent2)
    for e=1:length(dataAll)
        d1=sort(dataAll(e).pupil.diameter);
        d2=sort(dataAll(e).pupil.derivative);
        n1=length(d1);
        n11=int16(n1*(percent1/100));
        n2=length(d2);
        n21=int16(n2*(percent2/100));
        pupilSizeAll{e}=[d1(n11) d1(n1-n11) d2(n21) d2(n2-n21)];
    end
end


%%
function dataAll=loadData(conditions,protocolLabels,protocols,genotypes,pupilSmooth)
    % load in all the data and organize it into structures
    % these arrays contain the output data
    expStep=1;
    expFolders1=dir(fullfile('Ari_data_backup','WT*'));
    expFolders2=dir(fullfile('Ari_data_backup','Alpha*'));
    expFolders=[expFolders1;expFolders2];
    %expFolders=expFolders1(1);
    for e=1:length(expFolders)
        filename=expFolders(e).name;
        path1=fullfile(pwd,'Ari_data_backup',filename);
        pupilData=[];
        if filename(1:2)==genotypes{1}(1:2)
            g=1;
        elseif filename(1:2)==genotypes{2}(1:2)
            g=2;
        else
            g=0;
        end

        %first load in the pupillometry data
        pupFolder=dir(fullfile(path1,'*Pupillometry*'));
        path2=fullfile(path1,pupFolder.name);
        pupilData=loadPupillometry(path2);

        tableIndex=1;    
        %now load in the vibrometry data
        for c=1:length(conditions)
            path3=fullfile(path1,conditions{c});
            for p=1:length(protocols)         
                path4=dir(fullfile(path3,protocolLabels{p}));
                path5=dir(fullfile(path3,path4.name,'**','2021-*'));  %find all data files in subfolders
                toRemove = contains({path5.name}, '.xlsx'); 
                path5(toRemove)=[];  %remove excel files, just leaving the folder names
                                
                if length(path5)>0
                    for m=1:length(path5)
                        path6=fullfile(path5(m).folder,path5(m).name);
                        path7=dir(fullfile(path6,'analysis *'));
                        if isempty(path7)
                            disp(path6)
                        else
                            %get the time the data was collected
                            x1=split(path7.folder,filesep);
                            x2=split(x1(length(x1)),' ');
                            x3=split(x2(1),'-');
                            x4=split(x2(2),'_');
                            tVector=[str2double(x3{1}),str2double(x3{2}),str2double(x3{3}),str2double(x4{1}),str2double(x4{2}),str2double(x4{3})];

                            filename1=fullfile(path6,path7.name,'vibNal.mat');
                            load(filename1);
                            filename2=fullfile(path6,path7.name,'params.mat');
                            load(filename2);

                            %now analyze the vibrometry data and put the results into a table:
                            % datetime,CF,Q10dB,gain at CF, gain at ~1/2 CF, phase at CF,
                            nal=analyze_data(params,vib,e,c,p,m);
                            nal.condition=conditions{c};
                            nal.protocol=protocols{p};
                            nal.measurement=m;
                            nal.time=datetime(tVector);
                            nal.diameter=nan;
                            z=1;
                            vibData(tableIndex,:)=nal;
                            tableIndex=tableIndex+1;

                            %store start and stop time for each part of the experiment to
                            %sync with pupil data
                            if m==1
                                times(expStep).condition=conditions{c};
                                times(expStep).protocol=protocols{p};
                                times(expStep).start=datetime(tVector);
                            end
                            if m==length(path5)
                                times(expStep).end=datetime(tVector);
                            end 
                        end
                    end
                    expStep=expStep+1;
                end
            end
        end

        % remove NaNs, then smooth and downsample pupil data 
        v2=pupilData;
        v2=v2(~(isnan(v2.diameter)),:);
        v2{:,2}=movmean(v2{:,2},pupilSmooth);
        v2=downsample(v2,pupilSmooth);

        % fix the date of the pupil data (since only the time was saved correctly)
        nPupil=height(v2);    
        dVec=datevec(v2{:,1});
        dVec(:,1)=repmat(tVector(1),[1,nPupil]);
        dVec(:,2)=repmat(tVector(2),[1,nPupil]);
        dVec(:,3)=repmat(tVector(3),[1,nPupil]);
        dVec1=datetime(dVec);
        v2{:,1}=dVec1;
        x=seconds(v2{:,1}-v2{1,1});
        y=v2{:,2};
        y1=get_derivative(x,y);
        v2=addvars(v2,y1,'NewVariableNames','derivative');

        % now append this experiment onto the overall data table
        dataAll(e).experiment=e;
        dataAll(e).filename=filename;
        dataAll(e).genotype=genotypes{g};
        dataAll(e).vib=vibData;
        dataAll(e).pupil=v2;
        dataAll(e).times=times;
    end
end

function y1=get_derivative(x,y)
    f = spline(x,y);
    deriv = fnder(f,1);
    y1=ppval(deriv,x);
    plot(x,y,'b',x,y1,'r');
end

function out=createVarNames(letters,num)
    out={};
    for i=1:num
        out{i}=[letters num2str(i)];
    end
end

function pupilData=loadPupillometry(path2)
    path3=dir(fullfile(path2,'pupil_diameters_2*'));
    for file=1:length(path3)
        filename=fullfile(path2,path3(file).name);
        opts = detectImportOptions(filename);
        opts=setvartype(opts,'Date_Time','datetime');
        opts=setvartype(opts,'Diameter','double');
        T=readtable(filename,opts);
        if file==1
            pupilData=T;
        else
            pupilData=cat(1,pupilData,T);
        end
    end
    pupilData=sortrows(pupilData);
    pupilData=renamevars(pupilData,['Date_Time'],['datetime']);
    pupilData=renamevars(pupilData,['Diameter'],['diameter']);
end

function nal=analyze_data(params,vib,e,c,p,m)
    plotData=0;  % 0=do not create plots

    %% Generic q analysis
    qdB = 10; % dB down from maximum to calculate Q
    nsd = 2; % # SD above mean noise floor required to be considered 'clean' data
    use_clean_data = 1; % 0=use unclean data; 1=use only clean data in Q calculation
    min_clean_pointN = 3; % minimum # of clean data points required to attempt analysis
    CFlevel=20; %level to measure CF at
    gain1=[20 80]; %levels to measure gain1 at
    gain2=[50 70]; %levels to measure gain1 at
    
    levelArray=[10 20 30 40 50 60 70 80];

    xlims = [.98 15.2];

    % Parameters
    F1s = params.f1.f;
    F1N = length(F1s);
    L1s = params.f1.L;
    L1N = length(L1s);
    
    nans = nan(1,L1N); % dummy
    nansF= nan(1,F1N);
    nansVib=nan(F1N,L1N);
    
    % set up nal structure
    nal.L = L1s; % level (dB SPL)
    nal.F = F1s; % frequencies
    nal.bf = nans; % BF (Hz)
    nal.bf_mag = nans; % Peak magnitude (nm)
    nal.bfIndex=nans; %frequency index for the BF
    nal.CF=nan; % CF (BF at a CFlevel)
    nal.CFIndex=nan; % frequency index for the CF
    nal.bw=nans; %bandwidth
    nal.q = nans; % Q
    nal.mag=nansVib;
    nal.phi=nansVib;
    nal.gain1=nansF;
    nal.gain2=nansF;
   
    % Clean magnitudes
    mags = vib.f1.mag;
    nfs = vib.f1.nf + nsd.* vib.f1.nfsd;
    magsC = mags;
    magsC(mags<nfs)=NaN;

    % Clean and adjust phases
    phiTemp=[zeros([1 L1N]);vib.f1.phi]; % put in a row of zeros first to make sure phases start off near zero
    phisTemp = unwrap(phiTemp); % phase unwrap along freq axis first.
    phis=phisTemp(2:end,:); % remove the zeros
    phis=flip(unwrap(flip(phis,2),[],2),2); % Unwrap phases along stim intensity (backwards from high to low)   
    phisC = phis; 
    phisC(mags<nfs)=NaN;
    
    %convert to cycles
    phisC=phisC/(2*pi); % phase (cycles)
    phis=phis/(2*pi); % phase (cycles)

    % some experiments didn't use standard parameters. Change them so they
    % all fit into one big table
    if F1N==19
        F1N=20;
        F1s(20)=14375;
        nal.F = F1s; % frequencies
        mags(20,:)=nan;
        magsC(20,:)=nan;
        phis(20,:)=nan;
        phisC(20,:)=nan;
    end

    if L1N<8 || L1N>8
        magsTemp=ones(20,8)*nan;
        magsCTemp=magsTemp;
        phisTemp=magsTemp;
        phisCTemp=magsTemp;
        [Lia,Locb]=ismember(L1s,levelArray);
        for i=1:length(Lia)
            if Lia(i)
                magsTemp(:,Locb(i))=mags(:,i);
                magsCTemp(:,Locb(i))=magsC(:,i);
                phisTemp(:,Locb(i))=phis(:,i);
                phisCTemp(:,Locb(i))=phisC(:,i);
            end
        end
        clear mags magsC phis phisC L1s L1N
        mags=magsTemp;
        magsC=magsCTemp;
        phis=phisTemp;
        phisC=phisCTemp;
        L1s = levelArray;
        L1N = length(L1s);
        nal.L=L1s;
        nans=nan(1,L1N);
        nal.bf = nans; % BF (Hz)
        nal.bf_mag = nans; % Peak magnitude (nm)
        nal.bfIndex=nans; %frequency index for the BF
        nal.bw=nans; %bandwidth
        nal.q = nans; % Q
    end

    % set frequency to work with
    f = F1s;
    fN = length(f);

    if plotData==1
        %% Plot
        h=figure('units','normalized','position',[.2 .2 .5 .7]);
        ax1 = axes('position',[.1 .55 .35 .35], 'box','off','LineWidth',1.8, 'FontSize',14); hold on; % Displacement mag. 
        ax2 = axes('position',[.1 .1 .35 .35], 'box','off','LineWidth',1.8, 'FontSize',14); hold on; % Displacement phase (cycles)
        ax3 = axes('position',[.6 .55 .35 .35], 'box','off','LineWidth',1.8, 'FontSize',14); hold on; % Displacement delay (ms)

        axes(ax1); % Displacement (dB re 1 nm)
        xlim(xlims); set(gca,'Xscale','log','Xtick',[1 2 5 10 20 50]);
        ylim([-40 45]); set(gca,'Ytick',-100:20:100);
        xlabel('Frequency (kHz)');ylabel('Displacement (dB re 1 nm)');

        axes(ax2); % Phase (cycles)
        xlim(xlims); set(gca,'Xscale','log','Xtick',[1 2 5 10 20 50]);
        ylim([-5 .5]);
        xlabel('Frequency (kHz)');ylabel('Phase (cycles)');

        axes(ax3); % Q
        xlim([10 100]); set(gca,'Xtick',[0:20:100]);
        ylim([0 5]);
        xlabel('Stimulus level (dB SPL)');ylabel('Q10dB');
    end
    
    %% Calculate params for each level
    for L_i = 1:L1N
        L1 = L1s(L_i);
%
        mag = squeeze(mags(:,L_i));
        magC = squeeze(magsC(:,L_i));
        phi = squeeze(phis(:,L_i));
        phiC = squeeze(phisC(:,L_i));

        %% Calculate Q
        if use_clean_data
            for f_i = 1:fN    
                if f_i == 1
                    if isnan(magC(f_i+1)) % if first point is clean but second is not
                        magC(f_i) = NaN;
                    end    
                elseif f_i == fN
                    if isnan(magC(f_i-1)) % if last point is clean but previous is not
                        magC(f_i) = NaN;
                    end
                else
                    if isnan(magC(f_i-1)) && isnan(magC(f_i+1)) % if points on either side are not clean
                        magC(f_i) = NaN;
                    end
                end
            end
        end

        % Calculate max
        peak_i = []; % re-set max index

        if use_clean_data
            magx = magC;
        else
            magx = mag;
        end
        
        clean_i = find(~isnan(magx));
        [peak_mag, peak_clean_i] = max(magx(clean_i)); % maximum in clean amplitudes
        peak_i = clean_i(peak_clean_i); % actual index in data (peak_i is index in non-nan values)
        if length(magx(clean_i)) < min_clean_pointN % if too few points are clean, re-set max index to empty
            peak_i = [];
        elseif peak_i==1||peak_i==F1N
            peak_i = [];
        elseif isnan(magx(peak_i-1)) && isnan(magx(peak_i+1)) % if points on either side aren't clean, re-set to empty
            peak_i = [];
        end

        if ~isempty(peak_i)
            nal.bf_mag(L_i) = peak_mag; % Peak mag. (linear units?)
            peak_f = f(peak_i);
            nal.bf(L_i) = peak_f; % Peak frequency
            nal.bfIndex(L_i)=peak_i; %index at peak freq

            % Calculate Q
            intMag = 20*log10(peak_mag) - qdB; % intercept magnitude, in dB (magnitude 'qdB' down from peak mag.)
            int1 = NaN; % low-intercept
            int2 = NaN; % high-intercept

            % Intercept 1: Reset
            ctr = peak_i; % reset counter to be the peak magnitude
            int_not_found = 1; % reset intercept not found as default

            while ctr > 1 && int_not_found
                if 20*log10(magx(ctr)) > intMag && 20*log10(magx(ctr-1)) < intMag
                    int1 = interp1(20*log10(magx(ctr-1:ctr)), f(ctr-1:ctr), intMag, 'linear');
                    int_not_found = 0;
                else
                    ctr = ctr-1;
                end
            end

            % Intercept 2: Reset
            ctr = peak_i;
            int_not_found = 1;

            while ctr < fN && int_not_found
                if 20*log10(magx(ctr)) > intMag && 20*log10(magx(ctr+1)) < intMag
                    int2 = interp1(20*log10(magx(ctr:ctr+1)), f(ctr:ctr+1), intMag, 'linear');
                    int_not_found = 0;   
                else
                    ctr = ctr+1;
                end
            end

            if ~isnan(int1) && ~isnan(int2) % if intercepts found, calculate BW and Q
                nal.bw(L_i) = int2-int1;
                nal.q(L_i) = peak_f / abs(nal.bw(L_i));  
            else
                nal.bw(L_i) = nan;
                nal.q(L_i) = nan;  
            end         
            
            if plotData==1
                %% Plot tuning curves w/ BF and intercepts
                axes(ax1);
                plot(f/1000, 20*log10(mag), 'k--', 'LineWidth', .5);
                plot(f/1000, 20*log10(magC), 'k-', 'LineWidth', 1, 'MarkerFaceColor', 'k');

                plot(peak_f/1000, 20*log10(peak_mag), 'ro', 'MarkerSize', 5, 'LineWidth', 2);
                if ~isnan(int1)
                   plot(int1/1000, intMag, 'bo', 'MarkerSize', 5, 'LineWidth', 2);
                end
                if ~isnan(int2)
                    plot(int2/1000, intMag, 'bo', 'MarkerSize', 5, 'LineWidth', 2);
                end

                axes(ax2);
                plot(f/1000, phi, 'k--', 'LineWidth', .5);
                plot(f/1000, phiC, 'k-', 'LineWidth', 1);
            end
        else  %if there is no peak
            nal.bf_mag(L_i)=nan;
            nal.bf(L_i)=nan;
            nal.bfIndex(L_i)=nan;
            nal.bw(L_i)=nan;
            nal.q(L_i)=nan;
        end
    end

    % calculate specific values to save into the table
    if nal.bf(find(L1s==CFlevel))>0
        nal.CF=nal.bf(find(L1s==CFlevel)); % CF (BF at a specific intensity (i.e. 20 dB SPL))
        nal.CFIndex=find(nal.CF==nal.F);
    end
    
    nal.mag=magsC;
    nal.phi=phisC;   
    nal.gain1=20*log10(nal.mag(:,find(L1s==gain1(1)))./nal.mag(:,find(L1s==gain1(2))))+20*log10(gain1(2)-gain1(1));  
    nal.gain2=20*log10(nal.mag(:,find(L1s==gain2(1)))./nal.mag(:,find(L1s==gain2(2))))+20*log10(gain2(2)-gain2(1)); 
    if plotData==1
        %% Plot Q vs. level
        axes(ax3);
        plot(L1s, nal.q, 'o-', 'LineWidth', 1, 'Color', 'k', 'MarkerFaceColor','k'); hold on;    
    end
end


function distinctColors = ptc6(n, varargin)
% PTC6: 'Bright' color palettes from Paul Tol at SRON
% Creates a palette of n (between 1 and 6) distinctive
% colors that are:
%  - distinct for all people, including colour-blind readers
%  - distinct from black and white
%  - distinct on screen and paper
%  - still match well together
%
% Output is (n,3) RGB colors
% A light grey is included at n+1 (can indicate absence of data)
% 
% Calling ptc6(n,'check') will show a figure displaying 
% the color scheme that is selected
% Colors from ptc6 are slightly brighter than those from ptc12
%
% All colors and schemes created by Paul Tol at SRON:
% https://personal.sron.nl/~pault/
% Written by Nans Bujan, PhD, 02/03/2019. 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if n <0 || n >6
    error('N must be comprised between 0 and 6')
end    

% Colors definition in RGB
PAULTOLCOLORS =     [68,    119,    170;...                 % blue
                    102,    204,    238;...                 % cyan   
                     34,    136,     51;...                 % green
                    204,    187,     68;...                 % yellow
                    238,    102,    119;...                 % red
                    170,     51,    119;...                 % pink
                    187,    187,    187]...                 % light grey
                    / 255; 

% Colors schemes 
PAULTOLSCHEMES = { ...
                [1], ...                                       
                [1, 5], ... 
                [1, 5, 3], ... 
                [1, 5, 3, 4], ...                           
                [1, 2, 3, 4, 5], ...                             
                [1, 2, 3, 4, 5, 6]};

distinctColors = PAULTOLCOLORS([PAULTOLSCHEMES{n},7],:);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Displays a figure to show the selected colors:
    if strcmp(varargin,'check')
        X       = [0:0.1:5]; % a vector to plot some SIN functions
        LW      = 6; % Width of plot lines
        figure
        
        sp1     = subplot(2,1,1);   % white background
            for  i = 1:n+1
                plot(X,sin(X+i/2),'Color',distinctColors(i,:),'LineWidth',LW)
                hold on
            end
            grid on
            title(['Showing ',num2str(n),' colors (and one light grey)'])
        
        sp2     = subplot(2,1,2);   % black background
            copyobj(get(sp1,'Children'),sp2);
            set(gca, 'color', 'k')
            grid on
    end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

