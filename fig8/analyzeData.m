% Program to load in Michele Pei's data and create time-stamped variables 
% of the awake mice data
% Pupil diameter, gain, Q10dB, mag, phase, etc. 

% The goal is for these variables to then be analyzed for trends, frequency
% variations, etc.
close all;clear;clc;
warning('off','all')
conditions={'Awake'};
protocolLabels={'*single-tone*'};
protocols={'Singletone'};

genotypes={'WT','Alpha9KO'};
pSize={'Small','Medium','Large'};
percent1=30; %the top and bottom percent of pupil diameters for thresholding
percent2=50; %the top and bottom percent of pupil diameter derivatives for thresholding

%%
data1=loadData1(genotypes);
data1=rmmissing(data1,'DataVariables',{'experiment','time','diameter'});
experiments=unique(data1.experiment);

%%

%%
for e=1:length(experiments)
    data2=data1(strcmp(data1.experiment,experiments(e)),:);
    [C,ia] = unique(data2.time);
    data2=data2(ia,:);
    
    x=data2.time-data2.time(1);
    y=data2.diameter;
    y1=get_derivative(x,y);
    data2=addvars(data2,y1,'NewVariableNames','derivative');

    %now crop out the outliers
    d1=sort(data2.F10000);
    d1=d1(~isnan(d1));
    n1=length(d1);
    n11=int16(n1*(15/100));
    n12=n1-int16(n1*(7/100));
    ind1=find(data2.F10000>d1(n11)&data2.F10000<d1(n12));
    data2=data2(ind1,:);

    %now determine the small,medium, and large pupils
    clear d1 d2 n1 n2 n11 n12 ;   
    d1=sort(data2.diameter);
    d2=sort(data2.derivative);
    n1=length(d1);
    n11=int16(n1*(percent1/100));
    n2=length(d2);
    n21=int16(n2*(percent2/100));
    pupilSize=[d1(n11) d1(n1-n11) d2(n21) d2(n2-n21)];
    data2.size(data2.diameter<=pupilSize(1))=pSize(1);
    data2.size(data2.diameter>pupilSize(1)&data2.diameter<pupilSize(2))=pSize(2);
    data2.size(data2.diameter>=pupilSize(2))=pSize(3);
    if e==1
        data3=data2;      
    else
        data3=[data3;data2];
    end
    
    if e==6
        time1=zeros(1,length(data2.time));
        time1(1)=0;
        for i=2:length(data2.time)
            tDiff=data2.time(i)-data2.time(i-1);
            if tDiff>2
                tDiff=2;
            end
            time1(i)=time1(i-1)+tDiff;
        end
        diam1=data2.diameter;
        mag1=smooth(data2.F10000,4);
        pupilSize1=pupilSize;
    end
end

%%
if length(experiments)>5
    p1=2;
    p2=7;
elseif length(experiments)>1
    p1=5;
    p2=1;
else
    p1=1;
    p2=1;
end

delayTime=1; %time delay between vibrometry measurement and pupil measurement (positive is pupil after vibrometry)
minNum=50; %minimum number of points to include in averaging
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
iTable=1;

% now calculate averaged data for small, medium, and large pupils
% for each experiment
c=1;
p=1;
F=[5000,6000,7000,8000,9000,10000,11000,12000,13000];
L=[50];
for e=1:length(experiments)
    clear mag phi gain1 gain 2 Q BF diameter;
    data4=data3(strcmp(data3.experiment,experiments(e)),:);

    idx{1}=strcmp(data4.size,pSize(1));
    idx{2}=strcmp(data4.size,pSize(2));
    idx{3}=strcmp(data4.size,pSize(3));
    for s=1:3
        data(c,p,e,s).size=pSize(s);
        data(c,p,e,s).F=F;
        data(c,p,e,s).Fmat=repmat(F,[length(L) 1])';
        data(c,p,e,s).L=L;

        mag=table2array(data4(idx{s},5:13));
        data(c,p,e,s).mag=mean(mag,'omitnan');
        data(c,p,e,s).Q=mean(data4.q(idx{s}),'omitnan');
        data(c,p,e,s).BF=mean(data4.bf(idx{s}),'omitnan');
        data(c,p,e,s).diameter=mean(data4.diameter(idx{s}),'omitnan');

        data(c,p,e,s).magSD=std(mag,'omitnan');
        data(c,p,e,s).QSD=std(data4.q(idx{s}),'omitnan');
        data(c,p,e,s).BFSD=std(data4.bf(idx{s}),'omitnan');
        data(c,p,e,s).diameterSD=std(data4.diameter(idx{s}),'omitnan');

        dim=size(mag);
        data(c,p,e,s).magN=dim(1);
        data(c,p,e,s).QN=length(data4.q(idx{s})) - sum(isnan(data4.q(idx{s})));
        data(c,p,e,s).BFN=length(data4.bf(idx{s})) - sum(isnan(data4.bf(idx{s})));
        data(c,p,e,s).diameterN=length(data4.diameter(idx{s})) - sum(isnan(data4.diameter(idx{s})));
        data(c,p,e,s).genotype=data4.genotype(1);
        
        % mag/phase
        temp=reshape(mag',[],1);
        npts=length(temp);
        v1Temp=table('Size',[npts 14], ...
            'VariableTypes',{'int32','string','string','string','string','double', ...
            'double','double','double','double', 'double', ...
            'double','double','double'},'VariableNames', ...
            {'experiment','condition','protocol','genotype','pupilSize', ...
            'freq','level','mag','phi','gain1','gain2','q','bf','diameter'});
        v1Temp.experiment=repmat(e,[npts 1]);
        v1Temp.condition=repmat(conditions(c),[npts 1]);
        v1Temp.protocol=repmat(protocols(p),[npts 1]);
        v1Temp.genotype=repmat(data(c,p,e,s).genotype,[npts 1]);
        v1Temp.pupilSize=repmat(pSize(s),[npts 1]);
        v1Temp.freq=reshape(repmat(F,[length(L) size(mag,1)])',[],1);
        v1Temp.level=reshape(repmat(L,[length(F) size(mag,1)]),[],1);
        v1Temp.mag=temp;
        v1Temp.phi=temp;
        
        % q, bf
        temp=reshape(data4.q(idx{s}),[],1);
        npts=length(temp);
        v3Temp=table('Size',[npts 14], ...
            'VariableTypes',{'int32','string','string','string','string','double', ...
            'double','double','double','double', 'double', ...
            'double','double','double'},'VariableNames', ...
            {'experiment','condition','protocol','genotype','pupilSize', ...
            'freq','level','mag','phi','gain1','gain2','q','bf','diameter'});
        v3Temp.experiment=repmat(e,[npts 1]);
        v3Temp.condition=repmat(conditions(c),[npts 1]);
        v3Temp.protocol=repmat(protocols(p),[npts 1]);
        v3Temp.genotype=repmat(data(c,p,e,s).genotype,[npts 1]);
        v3Temp.pupilSize=repmat(pSize(s),[npts 1]);
        v3Temp.level=reshape(repmat(L,[size(data4.q(idx{s}),1) 1]),[],1);
        v3Temp.q=temp;
        v3Temp.bf=reshape(data4.bf(idx{s}),[],1);

        % diameter
        temp=data4.diameter(idx{s});
        npts=length(temp);
        v4Temp=table('Size',[npts 14], ...
            'VariableTypes',{'int32','string','string','string','string','double', ...
            'double','double','double','double', 'double', ...
            'double','double','double'},'VariableNames', ...
            {'experiment','condition','protocol','genotype','pupilSize', ...
            'freq','level','mag','phi','gain1','gain2','q','bf','diameter'});
        v4Temp.experiment=repmat(e,[npts 1]);
        v4Temp.condition=repmat(conditions(c),[npts 1]);
        v4Temp.protocol=repmat(protocols(p),[npts 1]);
        v4Temp.genotype=repmat(data(c,p,e,s).genotype,[npts 1]);
        v4Temp.pupilSize=repmat(pSize(s),[npts 1]);
        v4Temp.diameter=temp;

        v1=[v1;v1Temp];
        v3=[v3;v3Temp];
        v4=[v4;v4Temp];
    end
end                    

%average data across experiments
for s=1:3
    for g=1:2
        index=[];
        tempMag=[];
        tempQ=[];
        tempBF=[];
        tempDiameter=[];
        for e=1:length(experiments)
            if strcmp(genotypes{g},data(c,p,e,s).genotype)
                index=[index e];
                tempMag=[tempMag;data(c,p,e,s).mag];
                tempQ=[tempQ;data(c,p,e,s).Q];
                tempBF=[tempBF;data(c,p,e,s).BF];
                tempDiameter=[tempDiameter;data(c,p,e,s).diameter];
            end
        end    
        dataAve1(c,p,g,s).mag=mean(tempMag,'omitnan');
        dataAve1(c,p,g,s).magSD=std(tempMag,'omitnan');
        dataAve1(c,p,g,s).magN=size(tempMag,1);
        dataAve1(c,p,g,s).Q=mean(tempQ,'omitnan');
        dataAve1(c,p,g,s).QSD=std(tempQ,'omitnan');
        dataAve1(c,p,g,s).QN=length(tempQ) - sum(isnan(tempQ));
        dataAve1(c,p,g,s).BF=mean(tempBF,'omitnan');
        dataAve1(c,p,g,s).BFSD=std(tempBF,'omitnan');
        dataAve1(c,p,g,s).BFN=length(tempBF) - sum(isnan(tempBF));
        dataAve1(c,p,g,s).diameter=mean(tempDiameter,'omitnan');
        dataAve1(c,p,g,s).diameterSD=std(tempDiameter,'omitnan');
        dataAve1(c,p,g,s).diameterN=length(tempDiameter) - sum(isnan(tempDiameter));
        dataAve1(c,p,g,s).genotype=genotypes{g};
    end
end



plotData=1;  % 0=do not create plots
% now plot summary figure
if plotData==1 & c==1 & (p==1) 
    xlims = [5 13];
    ylims = [1 21];
    capsize=3;
    TOLC_base=viridis(7); % get distinct colors from a palette that works for colorblind
    TOLC=[TOLC_base(2,:);TOLC_base(4,:);TOLC_base(6,:)];
    f1=figure(1);
    f1.Position=[20,70,650,650];

    t=tiledlayout(4,4);
    nexttile([1,2]);
    plot(time1,diam1,'k-','LineWidth',2);
    xlabel('Time (s)');
    ylabel('Pupil (mm)');
    xlim([0 105]);
    xtickangle(0);

    nexttile([2 2]);
    g=1;
    x=data(c,p,e,1).Fmat/1000;
    y=(dataAve1(c,p,g,1).mag);
    z=(dataAve1(c,p,g,1).magSD./sqrt(dataAve1(c,p,g,1).magN));
    errorbar(x,y,z,'-', 'LineWidth', 2,'Color',TOLC(1,:),'CapSize',capsize)
    hold on;
    y=(dataAve1(c,p,g,2).mag);
    z=(dataAve1(c,p,g,2).magSD./sqrt(dataAve1(c,p,g,2).magN));
    errorbar(x,y,z,'-', 'LineWidth', 2,'Color',TOLC(2,:),'CapSize',capsize)
    y=(dataAve1(c,p,g,3).mag);
    z=(dataAve1(c,p,g,3).magSD./sqrt(dataAve1(c,p,g,3).magN));
    errorbar(x,y,z,'-', 'LineWidth', 2,'Color',TOLC(3,:),'CapSize',capsize)
    xlim(xlims); set(gca,'Xscale','log','Xtick',[1 2 5 10 20 50]);
    ylim(ylims); set(gca,'Yscale','log','Ytick',[0.1 0.2 0.5 1 2 5 10 20 50]);
    xlabel('Frequency (kHz)');ylabel('Magnitude (nm)');
    xticks([5,6,7,8,9,10,12]);
    text(11,15,'ns','Color','black')
%     annotation('textarrow', [0.841 0.841], [0.8 0.88],'LineWidth', 2);
    xtickangle(0);
    G=genotypes{g};
    title(G);
    text(5.1,17,'Small','Color',TOLC(1,:))
    text(5.1,13,'Medium','Color',TOLC(2,:))
    text(5.1,10,'Large','Color',TOLC(3,:))

    nexttile([1 2]);
    plot(time1,mag1,'k-','LineWidth',2);
    xlabel('Time (s)');
    ylabel('Magnitude (nm)');
    xlim([0 105]);
    xtickangle(0);
    
    nexttile([2 2])
    plot(diam1,mag1,'ko','MarkerFaceColor','k');
    ylabel('Magnitude (nm)');
    xlabel('Pupil (mm)');
    ylim([6 35]);
    hold on;
    y1=[6 35];
    x1=[pupilSize1(1) pupilSize1(1)];
    x2=[pupilSize1(2) pupilSize1(2)];
    plot(x1,y1,'k-.', 'LineWidth', 2);
    plot(x2,y1,'k-.', 'LineWidth', 2);
    text(0.535,33.5,'Small','Color',TOLC(1,:),'FontSize',9)
    text(0.73,33.5,'Med','Color',TOLC(2,:),'FontSize',9)
    text(0.90,33.5,'Large','Color',TOLC(3,:),'FontSize',9)
    xtickangle(0);
    hold off;
    mdl=fitlm(time1,mag1)

    nexttile([2 2]);
    g=2;
    x=data(c,p,e,1).Fmat/1000;
    y=(dataAve1(c,p,g,1).mag);
    z=(dataAve1(c,p,g,1).magSD./sqrt(dataAve1(c,p,g,1).magN));
    errorbar(x,y,z,'-', 'LineWidth', 2,'Color',TOLC(1,:),'CapSize',capsize)
    hold on;
    y=(dataAve1(c,p,g,2).mag);
    z=(dataAve1(c,p,g,2).magSD./sqrt(dataAve1(c,p,g,2).magN));
    errorbar(x,y,z,'-', 'LineWidth', 2,'Color',TOLC(2,:),'CapSize',capsize)
    y=(dataAve1(c,p,g,3).mag);
    z=(dataAve1(c,p,g,3).magSD./sqrt(dataAve1(c,p,g,3).magN));
    errorbar(x,y,z,'-', 'LineWidth', 2,'Color',TOLC(3,:),'CapSize',capsize)
    xlim(xlims); set(gca,'Xscale','log','Xtick',[1 2 5 10 20 50]);
    ylim(ylims); set(gca,'Yscale','log','Ytick',[0.1 0.2 0.5 1 2 5 10 20 50]);
    xlabel('Frequency (kHz)');ylabel('Magnitude (nm)');
    xticks([5,6,7,8,9,10,12]);
    text(11,15,'ns','Color','black')
%     annotation('textarrow', [0.796 0.796], [0.3 0.38],'LineWidth', 2);
    xtickangle(0);
    G=genotypes{g};
    if G=="Alpha9KO"
        G="Alpha9^{-/-}";
    end

    title(G);

    fontsize(gcf,scale=1.5);
    t.TileSpacing = 'compact';
    t.Padding = 'compact';
end

v1(1,:)=[];
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

%%
function data1=loadData1(genotypes)
    %load in Michele's organized Excel data
    data1=readtable("PeiData1.xlsx");
end

function y1=get_derivative(x,y)
    f = spline(x,y);
    deriv = fnder(f,1);
    y1=ppval(deriv,x);
%     plot(x,y,'b',x,y1,'r');
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