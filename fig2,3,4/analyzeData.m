% @file x20230626_mscan_mouse_group_avg_tc.m
%
% Plots averaged tuning curves from the data saved from
% F:\deweyMscanAnalysis\2023-06-26\x20230626_mscan_extraction_across_mouse_groups.m
%
% @author mjs
%
% 2023-07-02
%
% Need to be in this directory for the program to run:
% E:\Work\2024\02-01--jso-lab-backup\Oghalai_vglut3xalpha9\code4John
%



close all;clear;clc;
warning('off','all')
p=pwd;  %get current directory
%ind=strfind(p,filesep); % find locations of the separators.
%baseDir=p(1:ind(end)-1);% move up one level
baseDir=p;

% this is the vglut3 vs cbatime domain data for the paper 
%td{1}=importdata(strcat(baseDir,filesep,"\awkbhvn4oghz\2023-06-27 20_03_51 MScan wtBMclk_td.mat"));
%td{2}=importdata(strcat(baseDir,filesep,"\awkbhvn4oghz\2023-06-14 22_47_25 MScan vglut3kohomoshftngp2BMclk_td.mat"));
td{1}=importdata(strcat(baseDir,filesep,"\2023-06-27 20_03_51 MScan wtBMclk_td.mat"));
td{2}=importdata(strcat(baseDir,filesep,"\2023-06-14 22_47_25 MScan vglut3kohomoshftngp2BMclk_td.mat"));

% another example of vglut3 vs cba with tuning curves and TD data, but not
% used in the paper
%exampleTD{1}=importdata(strcat(baseDir,filesep,"\clickdataCBAvsvglut3koHOMO\newsamedayrecs\092823\ms2\analysisPDFs\td\2023-09-28 21_01_06 MScan CBACaJBMclk_td.mat"));
%exampleTD{2}=importdata(strcat(baseDir,filesep,"\clickdataCBAvsvglut3koHOMO\newsamedayrecs\092823\ms1\analysisPDFs\td\2023-09-28 17_27_57 MScan vglut3koHOMOBMclk_td.mat"));
%exampleTC{1}.filename=strcat(baseDir,filesep,'\clickdataCBAvsvglut3koHOMO\newsamedayrecs\092823\ms2\2023-09-28 21_30_19 MScan CBACaJBM\analysis 21_30_19');
%exampleTC{2}.filename=strcat(baseDir,filesep,'\clickdataCBAvsvglut3koHOMO\newsamedayrecs\092823\ms1\2023-09-28 17_35_34 MScan vglut3koHOMOBM\analysis 17_35_34');
exampleTD{1}=importdata(strcat(baseDir,filesep,"\2023-09-28 21_01_06 MScan CBACaJBMclk_td.mat"));
exampleTD{2}=importdata(strcat(baseDir,filesep,"\2023-09-28 17_27_57 MScan vglut3koHOMOBMclk_td.mat"));
exampleTC{1}.filename=strcat(baseDir,filesep,'\analysis 21_30_19');
exampleTC{2}.filename=strcat(baseDir,filesep,'\analysis 17_35_34');
exampleTC{1}.title='WT 2023-09-28 21-30-19';
exampleTC{2}.title='VGLUT3^{-/-} 2023-09-28 17-35-34';

%just need to run these functions once:
%scc_mscanNal_jso(exampleTC{1}.filename);
%scc_mscanNal_jso(exampleTC{2}.filename);
% then run scc_tcNalAll_mjsMod_jso.m to analyze the example data (need to
% comment in beginning of program

for i=1:2
    % Load the mat files params.mat, tcNal.mat, vibNal.mat for the example
    % data figure
    params = importdata(fullfile(exampleTC{i}.filename, 'params.mat'));
    tcNal = importdata(fullfile(exampleTC{i}.filename, 'tcNal.mat'));
    vibNal = importdata(fullfile(exampleTC{i}.filename, 'vibNal.mat'));
    exampleTC{i}.params=params;
    exampleTC{i}.tcNal=tcNal;
    exampleTC{i}.vibNal=vibNal;
    [mags, nfs, phis, magsC, phisC, gains, gainsC] = vibNal_post_proc(params, vibNal);
    exampleTC{i}.mag=magsC;
    exampleTC{i}.phi=phisC+1;
    exampleTC{i}.sens=gainsC;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Path to the extracted data saved with the earlier script. This will have
% extracted data organized according to the user-defined mouse groups.
p=pwd;  %get current directory
%ind=strfind(p,filesep); % find locations of the separators.
%basedir=p(1:ind(end)-1);% move up one level
basedir=p
%extracted_data_path = strcat(basedir,filesep,'analysis',filesep,'extracted_data.mat');
extracted_data_path = strcat(basedir,filesep,'extracted_data.mat');

option__include_scatter = true;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d = importdata(extracted_data_path);
d_in = d.group_data;

mouse_groups = fieldnames(d_in);
mouse_groups1 = {"WT","Alpha9^{-/-}","VGLUT3^{-/-}","VGLUT3^{-/-}Alpha9^{-/-}"};
gain1=[20 80]; %levels to measure gain1 at
repMouse=[7,16,23,31]; % the ppt slide number of each mouse within each group to plot as a representative tuning curve
repMouse(2:4)=repMouse(2:4)-numel(d_in.(mouse_groups{1}).tcNal); %subtract out the number of mice in the previous groups
repMouse(3:4)=repMouse(3:4)-numel(d_in.(mouse_groups{2}).tcNal);
repMouse(4)=repMouse(4)-numel(d_in.(mouse_groups{3}).tcNal);
repTitles={'WT 2023-06-27','Alpha9^{-/-} ms2-2023-06-21','VGLUT3^{-/-} ms2-2023-06-14','VGLUT3^{-/-}Alpha9^{-/-} 2023-07-13'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construct a table of all data to make further data manipulation easier.
tab_all = [];
nGroups=numel(mouse_groups);
[nFreq,nLevel]=size(d_in.(mouse_groups{1}).vibNal(1).f1.mag);
freqs=d_in.(mouse_groups{1}).vibNal(1).f1.f';
levels=d_in.(mouse_groups{1}).tcNal(1).L1s;
freq_array=repmat(freqs,1,nLevel);
level_array=repmat(levels,nFreq,1);

tempmag=cell(1,nGroups);
tempphi=cell(1,nGroups);
tempsens=cell(1,nGroups);
tempgain1=cell(1,nGroups);
tempQ=cell(1,nGroups);
tempBF=cell(1,nGroups);
tempgainMax=cell(1,nGroups);
tempCF=cell(1,nGroups);
tempQ10db=cell(1,nGroups);
tempSensLow=cell(1,nGroups);
tempSensHigh=cell(1,nGroups);

% set up tables for saving the data for later stats in R
v1=table('Size',[1 6], ...
    'VariableTypes',{'int32','string','double', ...
    'double','double','double'},'VariableNames', ...
    {'experiment','genotype', ...
    'freq','level','mag','phi'});
v2=table('Size',[1 4], ...
    'VariableTypes',{'int32','string','double', ...
    'double'},'VariableNames', ...
    {'experiment','genotype','freq','gain'});           
v3=table('Size',[1 5], ...
    'VariableTypes',{'int32','string','double', ...
    'double','double'},'VariableNames', ...
    {'experiment','genotype','level','q','bf'});
v4=table('Size',[1 5], ...
    'VariableTypes',{'int32','string','double', ...
    'double','double'},'VariableNames', ...
    {'experiment','genotype','maxGain', 'sensCF','sens5k'});

%%
mIndex=0;
for g = 1 : nGroups
    % g: iterate over mouse groups
    nAve=numel(d_in.(mouse_groups{g}).tcNal);
    tempmag{g}=nan(nFreq,nLevel,nAve);
    tempphi{g}=nan(nFreq,nLevel,nAve);
    tempsens{g}=nan(nFreq,nLevel,nAve);
    tempgain1{g}=nan(nFreq,nAve);
    tempQ{g}=nan(nLevel,nAve);
    tempBF{g}=nan(nLevel,nAve);
    tempgainMax{g}=nan(nAve,1);
    tempCF{g}=nan(nAve,1);
    tempQ10dB{g}=nan(nAve,1);
    tempSensLow{g}=nan(nAve,1);
    tempSensHigh{g}=nan(nAve,1);
    grp={};
    for m = 1:nAve
        mIndex=mIndex+1;
        % m: iterate over each MScan from a mouse. 
        grp{m}=mouse_groups{g};
        % get the tuning curves.
        params = d_in.(mouse_groups{g}).params(m);
        vib = d_in.(mouse_groups{g}).vibNal(m);
        [mags, nfs, phis, magsC, phisC, gains, gainsC] = ...
            vibNal_post_proc(params, vib);        
        phis=phis+1; % for some reason the phase is one cycle too low
        phisC=phisC+1;
        % look for the representative mouse in each group and store their
        % tuning curve data
        if m==repMouse(g)
            repData(g).mag=magsC;
            repData(g).phi=phisC;
            repData(g).sens=gainsC;
        end
        
        % ave data for genotype comparison
        z=magsC;
        if length(z)>0
            tempmag{g}(:,:,m)=z;
        end
        z=phisC;
        if length(z)>0
            tempphi{g}(:,:,m)=z;
        end
        z=gainsC;
        if length(z)>0
            tempsens{g}(:,:,m)=z;
            tempSensLow{g}(m)=max(max(gainsC(14:18,1:4)));  % max sens around 8500 Hz
            tempSensHigh{g}(m)=z(10,11);  %5000 at highest stim level Hz
        end
        z=20*log10(magsC(:,find(levels==gain1(1)))./magsC(:,find(levels==gain1(2))))+20*log10(gain1(2)-gain1(1));  
        if length(z)>0
            tempgain1{g}(:,m)=z;
        end
        z=d_in.(mouse_groups{g}).tcNal(m).q;
        if length(z)>0
            tempQ{g}(:,m)=z;
        end
        z=d_in.(mouse_groups{g}).tcNal(m).bf;
        if length(z)>0
            tempBF{g}(:,m)=z;
        end
        z=d_in.(mouse_groups{g}).tcNal(m).bf_gain_max;
        if length(z)>0
            tempgainMax{g}(m)=z;
        end
        z=d_in.(mouse_groups{g}).tcNal(m).bf_freq;
        if length(z)>0
            tempCF{g}(m)=z;
        end
        z=d_in.(mouse_groups{g}).tcNal(m).bf_q10db;
        if length(z)>0
            tempQ10dB{g}(m)=z;
        end
        
        xTable1=table(cellstr(repmat(mouse_groups{g},nFreq,1)),repmat(m,nFreq,1),freqs,tempgain1{g}(:,m),'VariableNames',["Group","Mouse","Freq","Gain3070"]);
        xTable1.Group=categorical(xTable1.Group,mouse_groups);
        if exist('pltTable1')
            pltTable1=[pltTable1; xTable1];
        else
            pltTable1=xTable1;
            
        end    
        
        xTable2=table(cellstr(repmat(mouse_groups{g},nLevel,1)),repmat(m,nLevel,1),levels',tempBF{g}(:,m),tempQ{g}(:,m),'VariableNames',["Group","Mouse","Level","BF","Q10dB"]);
        xTable2.Group=categorical(xTable2.Group,mouse_groups);
        if exist('pltTable2')
            pltTable2=[pltTable2; xTable2];
        else
            pltTable2=xTable2;
        end  

        % now fill up the tables to save the data for R stats
        % mag/phase
        temp=reshape(tempmag{g}(:,:,m),[],1);
        npts=length(temp);
        v1Temp=table('Size',[npts 6], ...
            'VariableTypes',{'int32','string','double', ...
            'double','double','double'},'VariableNames', ...
            {'experiment','genotype', ...
            'freq','level','mag','phi'});
        v1Temp.experiment=repmat(mIndex,[npts 1]);
        v1Temp.genotype=repmat(grp{m},[npts 1]);
        v1Temp.freq=reshape(repmat(freqs,[length(levels) 1])',[],1);
        v1Temp.level=reshape(repmat(levels,[length(freqs) 1]),[],1);
        v1Temp.mag=temp;
        v1Temp.phi=reshape(tempphi{g}(:,:,m),[],1);

        % gain1/2
        temp=reshape(tempgain1{g}(:,m),[],1);
        npts=length(temp);
        v2Temp=table('Size',[npts 4], ...
            'VariableTypes',{'int32','string','double', ...
            'double'},'VariableNames', ...
            {'experiment','genotype','freq','gain'});           
        v2Temp.experiment=repmat(mIndex,[npts 1]);
        v2Temp.genotype=repmat(grp{m},[npts 1]);
        v2Temp.freq=freqs;
        v2Temp.gain=temp;

        % q, bf
        temp=reshape(tempQ{g}(:,m),[],1);
        npts=length(temp);
        v3Temp=table('Size',[npts 5], ...
            'VariableTypes',{'int32','string','double', ...
            'double','double'},'VariableNames', ...
            {'experiment','genotype','level','q','bf'});
        v3Temp.experiment=repmat(mIndex,[npts 1]);
        v3Temp.genotype=repmat(grp{m},[npts 1]);
        v3Temp.level=levels';
        v3Temp.q=temp;
        v3Temp.bf=reshape(tempBF{g}(:,m),[],1);

        % misc
        npts=1;
        v4Temp=table('Size',[npts 5], ...
            'VariableTypes',{'int32','string','double', ...
            'double','double'},'VariableNames', ...
            {'experiment','genotype','maxGain', 'sensCF','sens5k'});
        v4Temp.experiment=mIndex;
        v4Temp.genotype=grp{m};
        v4Temp.maxGain=tempgainMax{g}(m);
        v4Temp.sensCF=tempSensLow{g}(m);
        v4Temp.sens5k=tempSensHigh{g}(m);

        v1=[v1;v1Temp];
        v2=[v2;v2Temp];
        v3=[v3;v3Temp];
        v4=[v4;v4Temp];
    end
    
    dataAve1(g).mag=mean(tempmag{g},3,'omitnan');
    dataAve1(g).magSD=std(tempmag{g},0,3,'omitnan');
    dataAve1(g).magN=size(tempmag{g},3) - sum(isnan(tempmag{g}),3);
    dataAve1(g).phi=mean(tempphi{g},3,'omitnan');
    dataAve1(g).phiSD=std(tempphi{g},0,3,'omitnan');
    dataAve1(g).phiN=size(tempphi{g},3) - sum(isnan(tempphi{g}),3);
    dataAve1(g).sens=mean(tempsens{g},3,'omitnan');
    dataAve1(g).sensSD=std(tempsens{g},0,3,'omitnan');
    dataAve1(g).sensN=size(tempsens{g},3) - sum(isnan(tempsens{g}),3);
    dataAve1(g).gain1=mean(tempgain1{g},2,'omitnan');
    dataAve1(g).gain1SD=std(tempgain1{g},0,2,'omitnan');
    dataAve1(g).gain1N=size(tempgain1{g},2) - sum(isnan(tempgain1{g}),2);
    dataAve1(g).Q=mean(tempQ{g},2,'omitnan');
    dataAve1(g).QSD=std(tempQ{g},0,2,'omitnan');
    dataAve1(g).QN=size(tempQ{g},2) - sum(isnan(tempQ{g}),2);
    dataAve1(g).BF=mean(tempBF{g},2,'omitnan');
    dataAve1(g).BFSD=std(tempBF{g},0,2,'omitnan');
    dataAve1(g).BFN=size(tempBF{g},2) - sum(isnan(tempBF{g}),2);
    dataAve1(g).gainMax=mean(tempgainMax{g},1,'omitnan');
    dataAve1(g).gainMaxSD=std(tempgainMax{g},0,1,'omitnan');
    dataAve1(g).gainMaxN=size(tempgainMax{g},1) - sum(isnan(tempgainMax{g}),1);
    dataAve1(g).CF=mean(tempCF{g},1,'omitnan');
    dataAve1(g).CFSD=std(tempCF{g},0,1,'omitnan');
    dataAve1(g).CFN=size(tempCF{g},1) - sum(isnan(tempCF{g}),1);
    dataAve1(g).Q10dB=mean(tempQ10dB{g},1,'omitnan');
    dataAve1(g).Q10dBSD=std(tempQ10dB{g},0,1,'omitnan');
    dataAve1(g).Q10dBN=size(tempQ10dB{g},1) - sum(isnan(tempQ10dB{g}),1);
    dataAve1(g).sensLow=mean(tempSensLow{g},1,'omitnan');
    dataAve1(g).sensLowSD=std(tempSensLow{g},0,1,'omitnan');
    dataAve1(g).sensLowN=size(tempSensLow{g},1) - sum(isnan(tempSensLow{g}),1);
    dataAve1(g).sensHigh=mean(tempSensHigh{g},1,'omitnan');
    dataAve1(g).sensHighSD=std(tempSensHigh{g},0,1,'omitnan');
    dataAve1(g).sensHighN=size(tempSensHigh{g},1) - sum(isnan(tempSensHigh{g}),1);
   
    xTable=table(grp',tempgainMax{g},tempCF{g},tempQ10dB{g},tempSensLow{g},tempSensHigh{g},'VariableNames',["Group","GainMax","CF","Q10dB","SensLow","SensHigh"]);
    if exist('pltTable')
        pltTable=[pltTable; xTable];
    else
        pltTable=xTable;
    end
        
end

%% Save the data for R analyses

% save data 
v1(1,:)=[];
v2(1,:)=[];
v3(1,:)=[];
v4(1,:)=[];
%savedDataPath = strcat(basedir,filesep,'analysis');
savedDataPath = basedir;
writetable(v1,strcat(savedDataPath,filesep,'v1'));
writetable(v2,strcat(savedDataPath,filesep,'v2'));
writetable(v3,strcat(savedDataPath,filesep,'v3'));
writetable(v4,strcat(savedDataPath,filesep,'v4'));

%% Plot example1 WT and VGLUT3KO - may not use this plot
f3=figure(1);
f3.Position=[20,70,1000,500];
t=tiledlayout(4,8);
TOLC_base=viridis(7); % get distinct colors from a palette that works for colorblind
TOLC=[TOLC_base(1,:);TOLC_base(5,:);TOLC_base(3,:);TOLC_base(7,:)];
Ylow=[42 42];
Yhigh=[75 75];
capsize=3;

xlims = [4 15];
ylims = [0.08 300];
ylimsS = [0 90];
ylimsP = [-3 .5];
Ytick_num=[0.1 0.3 1 3 10 30 100 300 1000];
x90=[4.2 4.2 4.2 4.2];
x90_2=[7 7 7 7];
y90_1=[200 160 160 200];
y90_2=[28 28 31 28];
x10=[8.4 8.8 9.2 7.8];
y10=[0.6 0.35 1.3 0.22];
y10_2=[79 77 86 85];
NF=1;     %noise floor
examp_groups=[1,3];
TD_level=30;
for g1 = 1 : 2
    g=examp_groups(g1);
    % representative mag
    nexttile(1+(g1-1)*8*2,[2,2]);
    x=freq_array/1000;
    y=repData(g).mag;
%     y(y<NF)=nan;
    plot(x,y,'-', 'Color',TOLC(g,:),'LineWidth', 2); 
    xlim(xlims); set(gca,'Xscale','log','Xtick',[1 2 4 5 6 7 8 10 12 15 20 50]);
    ylim(ylims); set(gca,'Yscale','log','Ytick',Ytick_num);
    xlabel('Freq (kHz)');ylabel('Magnitude (nm)');
    title(repTitles{g},'FontSize',8);
    text(x90(g),y90_1(g),'90','Color','k','FontSize',8)
    text(x10(g),y10(g),'10','Color','k','FontSize',8)
    xtickangle(0);

    % representative sensitivity 
    nexttile(3+(g1-1)*8*2,[2,2]);
    y=repData(g).sens;
    plot(x,y,'-', 'Color',TOLC(g,:),'LineWidth', 2); 
    xlim(xlims); set(gca,'Xscale','log','Xtick',[1 2 4 5 6 7 8 10 12 15 20 50]);
    ylim(ylimsS); set(gca,'Yscale','lin','Ytick',[0 20 40 60 80 100]);
    xlabel('Freq (kHz)');ylabel('Sens (dB re:1nm/Pa)');
    text(x90_2(g),y90_2(g),'90','Color','k','FontSize',8)
    text(x10(g),y10_2(g),'10','Color','k','FontSize',8)
    xtickangle(0);
    hold on;
    plot(xlims,Yhigh,'--', 'LineWidth', 1.5,'Color','r');
    hold off;

    % phase 
    nexttile(5+(g1-1)*8*2,[2,2]);
    y=repData(g).phi;
    plot(x,y,'-', 'Color',TOLC(g,:),'LineWidth', 2); 
    xlim(xlims); set(gca,'Xscale','log','Xtick',[1 2 4 5 6 7 8 10 12 15 20 50]);
    ylim(ylimsP); 
    xlabel('Freq (kHz)');ylabel('Phase (cycles)');   
    xtickangle(0);

    % time domain data
    xlimits=[10 12];
    ylimits=[-4 4];
    nexttile(7+(g1-1)*8*2,[2,2]);
    y=td{1,g1}.td_data.x30_dB_SPL;
    dt_s=td{1,g1}.dt_s;
    x=1000*[0:dt_s:(dt_s*(length(y)-1))]';
    ind=find((x>xlimits(1)) & (x<xlimits(2)));
    ave=mean(y(ind));
    y=y-ave;
    plot(x,y,'-', 'Color',TOLC(g,:),'LineWidth', 2); 
    xlim(xlimits);set(gca,'Xscale','lin');
    ylim(ylimits);
    xlabel('Time (ms)');ylabel('Displacement (nm)');   
    text(10.1,-3,'30 dB SPL click','Color','k','FontSize',8)

    xtickangle(0);
  
end
fontsize(f3,scale=1.5);
t.TileSpacing = 'compact';
t.Padding = 'compact';
shg

%% Plot example2 WT and VGLUT3KO - another example
f4=figure(4);
f4.Position=[20,70,1000,500];
t=tiledlayout(4,8);
TOLC_base=viridis(7); % get distinct colors from a palette that works for colorblind
TOLC=[TOLC_base(1,:);TOLC_base(5,:);TOLC_base(3,:);TOLC_base(7,:)];
Ylow=[42 42];
Yhigh=[74 74];
capsize=3;

xlims = [4 15];
ylims = [0.08 300];
ylimsS = [0 90];
ylimsP = [-3 .5];
Ytick_num=[0.1 0.3 1 3 10 30 100 300 1000];
x90=[4.2 4.2 4.2 4.2];
x90_2=[6.8 7 7 7];
y90_1=[170 160 190 200];
y90_2=[25 28 31 28];
x10=[8.4 8.8 8.9 7.8];
x10_2=[10.5 8.8 10.5 7.8];
y10=[0.65 0.35 2.5 0.22];
y10_2=[70 77 86 85];

NF=1;     %noise floor
examp_groups=[1,3];
TD_level=30;
for g1 = 1 : 2
    g=examp_groups(g1);
    % representative mag
    nexttile(1+(g1-1)*8*2,[2,2]);
    x=freq_array/1000;
    y=exampleTC{g1}.mag;
%     y(y<NF)=nan;
    plot(x,y,'-', 'Color',TOLC(g,:),'LineWidth', 2); 
    xlim(xlims); set(gca,'Xscale','log','Xtick',[1 2 4 5 6 7 8 10 12 15 20 50]);
    ylim(ylims); set(gca,'Yscale','log','Ytick',Ytick_num);
    xlabel('Freq (kHz)');ylabel('Magnitude (nm)');
    [ttt,sss] = title(exampleTC{g1}.title);
    ttt.FontSize = 6;
    ax = gca;
    ax.TitleHorizontalAlignment = 'left';
    text(x90(g),y90_1(g),'90','Color','k','FontSize',8)
    text(x10(g),y10(g),'10','Color','k','FontSize',8)
    xtickangle(0);

    % representative sensitivity 
    nexttile(3+(g1-1)*8*2,[2,2]);
    y=exampleTC{g1}.sens;
    plot(x,y,'-', 'Color',TOLC(g,:),'LineWidth', 2); 
    xlim(xlims); set(gca,'Xscale','log','Xtick',[1 2 4 5 6 7 8 10 12 15 20 50]);
    ylim(ylimsS); set(gca,'Yscale','lin','Ytick',[0 20 40 60 80 100]);
    xlabel('Freq (kHz)');ylabel('Sens (dB re:1nm/Pa)');
    text(x90_2(g),y90_2(g),'90','Color','k','FontSize',8)
    text(x10_2(g),y10_2(g),'10','Color','k','FontSize',8)
    xtickangle(0);
    hold on;
    plot(xlims,Yhigh,'--', 'LineWidth', 1.5,'Color','r');
    hold off;

    % phase 
    nexttile(5+(g1-1)*8*2,[2,2]);
    y=exampleTC{g1}.phi;
    plot(x,y,'-', 'Color',TOLC(g,:),'LineWidth', 2); 
    xlim(xlims); set(gca,'Xscale','log','Xtick',[1 2 4 5 6 7 8 10 12 15 20 50]);
    ylim(ylimsP); 
    xlabel('Freq (kHz)');ylabel('Phase (cycles)');   
    xtickangle(0);

    % time domain data
    xlimits=[10 12];
    ylimits=[-4 4.3];
    nexttile(7+(g1-1)*8*2,[2,2]);
    y=exampleTD{1,g1}.td_data.x30_dB_SPL;
    dt_s=exampleTD{1,g1}.dt_s;
    x=1000*[0:dt_s:(dt_s*(length(y)-1))]';
    ind=find((x>xlimits(1)) & (x<xlimits(2)));
    ave=mean(y(ind));
    y=y-ave;
    plot(x,y,'-', 'Color',TOLC(g,:),'LineWidth', 2); 
    xlim(xlimits);set(gca,'Xscale','lin');
    ylim(ylimits);
    xlabel('Time (ms)');ylabel('Displacement (nm)');   
    text(10.05,3.8,'30 dB SPL click','Color','k','FontSize',8)

    xtickangle(0);
  
end
fontsize(f4,scale=1.5);
t.TileSpacing = 'compact';
t.Padding = 'compact';
shg

%% Plot summary figure
c=1;
p=1;
f1=figure(3);
f1.Position=[20,70,750,1000];
t=tiledlayout(8,6);
TOLC_base=viridis(7); % get distinct colors from a palette that works for colorblind
TOLC=[TOLC_base(1,:);TOLC_base(5,:);TOLC_base(3,:);TOLC_base(7,:)];
Ylow=[42 42];
Yhigh=[75 75];

capsize=3;
xlims = [4 15];
ylims = [0.08 300];
ylimsS = [0 90];
ylimsP = [-3 .5];
Ytick_num=[0.1 0.3 1 3 10 30 100 300 1000];
% x90=[4.2 4.2 4.2 4.2];
x90=[6,6,6,6];
x90_2=[7 7 7 7];
y90_1=[120 120 120 120];
y90_2=[34 34 34 34];
x10=[8.4 8.8 9.3 7.8];
y10=[0.5 0.35 1.3 0.22];
% Plot tuning curves
NF=1;     %noise floor
for g = 1 : nGroups
    x=freq_array/1000;
%     %representative mag
%     nexttile(1+(g-1)*8*2,[2,2]);
%     x=freq_array/1000;
%     y=repData(g).mag;
% %     y(y<NF)=nan;
%     plot(x,y,'-', 'LineWidth', 1,'Color',TOLC(g,:),'LineWidth', 2); 
%     xlim(xlims); set(gca,'Xscale','log','Xtick',[1 2 4 5 6 7 8 9 10 12 15 20 50]);
%     ylim(ylims); set(gca,'Yscale','log','Ytick',Ytick_num);
%     xlabel('Freq (kHz)');ylabel('Magnitude (nm)');
%     title(repTitles{g});
%     text(x90(g),y90_1(g),'90','Color','k','FontSize',9)
%     text(x10(g),y10(g),'10','Color','k','FontSize',9)

    %averaged mag 
    nexttile(1+(g-1)*6*2,[2,2]);
    y=dataAve1(g).mag;
    z=(dataAve1(g).magSD./sqrt(dataAve1(g).magN));
    errorbar(x,y,z,'-', 'LineWidth', 1,'Color',TOLC(g,:),'CapSize',capsize, 'LineWidth', 2)
    xlim(xlims); set(gca,'Xscale','log','Xtick',[1 2 4 5 6 7 8 10 12 15 20 50]);
    ylim(ylims); set(gca,'Yscale','log','Ytick',Ytick_num);
    xtickangle(0);
    xlabel('Freq (kHz)');ylabel('Magnitude (nm)');
    G=mouse_groups{g};
    if G=="Alpha9KO"
        G="Alpha9^{-/-}";
    elseif G=="VGLUT3KO"
        G="VGLUT3^{-/-}";
    elseif G=="VGLUT3KOAlpha9KO"
        G="VGLUT3^{-/-}Alpha9^{-/-}";
    end
    title(G);
    txtString=sprintf('n=%i',max(max(dataAve1(g).magN)));
    text(11,200,txtString,'Color','k','FontSize',9)
    text(x90_2(g),y90_1(g),'90','Color','k','FontSize',9)

    %averaged sensitivity 
    nexttile(3+(g-1)*6*2,[2,2]);
    yS=dataAve1(g).sens;
    yS(dataAve1(g).sensN<3)=nan;
    z=(dataAve1(g).sensSD./sqrt(dataAve1(g).sensN));
    errorbar(x,yS,z,'-', 'LineWidth', 1,'Color',TOLC(g,:),'CapSize',capsize, 'LineWidth', 2);
    hold on;
    plot(xlims,Yhigh,'--', 'LineWidth', 1.5,'Color','r');
    hold off;
    xlim(xlims); set(gca,'Xscale','log','Xtick',[1 2 4 5 6 7 8 10 12 15 20 50]);
    ylim(ylimsS); set(gca,'Yscale','lin','Ytick',[0 20 40 60 80 100]);
    xtickangle(0);
    xlabel('Freq (kHz)');ylabel('Sensitivity (dB re:1nm/Pa)');
    G=mouse_groups{g};
    if G=="Alpha9KO"
        G="Alpha9^{-/-}";
    elseif G=="VGLUT3KO"
        G="VGLUT3^{-/-}";
    elseif G=="VGLUT3KOAlpha9KO"
        G="VGLUT3^{-/-}Alpha9^{-/-}";
    end
    title(G);
    text(x90(g),y90_2(g),'90','Color','k','FontSize',9)

    %averaged phase 
    nexttile(5+(g-1)*6*2,[2,2]);
    yP1=dataAve1(g).phi;
    yP1(y<NF)=nan;
    yP1(dataAve1(g).phiN<3)=nan;
    z=(dataAve1(g).phiSD./sqrt(dataAve1(g).phiN));
    errorbar(x,yP1,z,'-', 'LineWidth', 1,'Color',TOLC(g,:),'CapSize',capsize, 'LineWidth', 2);
    xlim(xlims); set(gca,'Xscale','log','Xtick',[1 2 4 5 6 7 8 10 12 15 20 50]);
    xtickangle(0);
    ylim(ylimsP); 
    xlabel('Freq (kHz)');ylabel('Phase (cycles)');   
    G=mouse_groups{g};
    if G=="Alpha9KO"
        G="Alpha9^{-/-}";
    elseif G=="VGLUT3KO"
        G="VGLUT3^{-/-}";
    elseif G=="VGLUT3KOAlpha9KO"
        G="VGLUT3^{-/-}Alpha9^{-/-}";
    end
    title(G);

end

fontsize(f1,scale=1.5);
t.TileSpacing = 'compact';
t.Padding = 'compact';
shg

%%
% Plot level dependence of gain/BF/Q10dB
f2=figure(2);
f2.Position=[1000,70,750,500];
t=tiledlayout(4,3);
nexttile(1,[2,1]);
capsize=3;
x=freqs/1000;
for g = 1 : nGroups
    y=(dataAve1(g).gain1);
    z=(dataAve1(g).gain1SD./sqrt(dataAve1(g).gain1N));
    errorbar(x,y,z,'-', 'LineWidth', 2,'Color',TOLC(g,:),'CapSize',capsize)
    hold on;
end
xlim([7 14]); set(gca,'Xscale','log','Xtick',[7 8 9 10 12 14 15]);
xtickangle(0);
ylim([17 36]); set(gca,'Yscale','linear','Ytick',[0 10 20 30 40 50]);
xlabel('Freq (kHz)');ylabel('Gain (dB)');
text(7.1,34.5,'20-80 dB SPL','Color','k','FontSize',8)

hold off;

nexttile(2,[2,1]);
x=levels;
for g = 1 : nGroups
    y=(dataAve1(g).BF)/1000;
    z=((dataAve1(g).BFSD/1000)./sqrt(dataAve1(g).BFN));
    errorbar(x,y,z,'-','Color',TOLC(g,:),'CapSize',capsize, 'LineWidth', 2)
    hold on;
end
xlim([0 90]); set(gca,'Xscale','linear','Xtick',[0 20 40 60 80 100]);
ylim([5 9.5]); set(gca,'Yscale','linear','Ytick',[5 6 7 8 9 10]);
xlabel('Level (dB SPL)');ylabel('BF (kHz)');
hold off;

nexttile(3,[2,1]);
x=levels;
for g = 1 : nGroups
    y=(dataAve1(g).Q);
    z=(dataAve1(g).QSD./sqrt(dataAve1(g).QN));
    errorbar(x,y,z,'-', 'Color',TOLC(g,:),'CapSize',capsize, 'LineWidth', 2)
    hold on;
    text(15,4.05-(g*0.22),mouse_groups1{g},'Color',TOLC(g,:),'FontSize',7)

end
xlim([0 90]); set(gca,'Xscale','linear','Xtick',[0 20 40 60 80 100]);
ylim([1 4]); set(gca,'Yscale','linear','Ytick',[0 1 2 3 4 5 6]);
xlabel('Level (dB SPL)');ylabel('Q10dB');
hold off;

% box and whisker plots of data at CF
bpcolor=[TOLC(1:4,:)];
pltTable.Group=categorical(pltTable.Group,mouse_groups);
m_g=categorical(mouse_groups,mouse_groups);
nexttile(7,[2,1]);
hold on;
for g=1:nGroups
    test_table=pltTable(pltTable.Group==m_g(g),:);
    boxchart(test_table.Group,test_table.GainMax,'BoxFaceColor',TOLC(g,:),'MarkerStyle','none','BoxFaceAlpha',0.7, 'LineWidth', 2)
end
hold off;
set(gca,'XTickLabel',mouse_groups1);
ylabel('Max Gain (dB)');
ylim([25 70]);

nexttile(8,[2,1]);
% boxchart(pltTable.Group,pltTable.CF/1000);
% ylabel('CF (kHz)');
hold on;
for g=1:nGroups
    test_table=pltTable(pltTable.Group==m_g(g),:);
    boxchart(test_table.Group,test_table.SensLow,'BoxFaceColor',TOLC(g,:),'MarkerStyle','none','BoxFaceAlpha',0.7, 'LineWidth', 2)
end
hold off; 
set(gca,'XTickLabel',mouse_groups1);
ylabel('Sens at CF (dB re:1nm/Pa)');
ylim([60 96]);

nexttile(9,[2,1]);
% boxchart(pltTable.Group,pltTable.Q10dB);
% ylabel('Q10dB');
hold on;
for g=1:nGroups
    test_table=pltTable(pltTable.Group==m_g(g),:);
    boxchart(test_table.Group,test_table.SensHigh,'BoxFaceColor',TOLC(g,:),'MarkerStyle','none','BoxFaceAlpha',0.7, 'LineWidth', 2)
end
hold off; 
set(gca,'XTickLabel',mouse_groups1);
ylabel('Sens at 5kHz (dB re:1nm/Pa)');
ylim([30 60]);

fontsize(f2,scale=1.5);
t.TileSpacing = 'compact';
t.Padding = 'compact';
shg
%% Save the data for R analyses
writetable(pltTable,'Stats_jso_boxplots');
writetable(pltTable1,'Stats_jso_gain');
writetable(pltTable2,'Stats_jso_bf_q10db');

%%%%
% Save figures f1, f2, f3, f4
fprintf('Saving f1\n');
saveas(f1, 'f1.fig');
fprintf('Saving f2\n');
saveas(f2, 'f2.fig');
fprintf('Saving f3\n');
saveas(f3, 'f3.fig');
fprintf('Saving f4\n');
saveas(f4, 'f4.fig');



%% color picker function
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
