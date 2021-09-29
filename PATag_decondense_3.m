function PATag_decondense_3

%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Author: Sébastien Huet
%   Email: sebastien.huet (at) univ-rennes1.fr
%   Date: 29-07-2021
%   
%   This is the main function used to assess chromatin relaxation at sites of damage
%
%%%%%%%%%%%%%%%%%%%%%%%%%


%%%Parameters
nbclust=2; %%%nb of clusters for the kmeans segmentation
signalclust=2; %index of the cluster containing the signal of interest
Rclose=2; %radius of the closing operator for cleaning the mask
Lengthclose=5; %length of the closing line for cleaning the mask (reassociate broken chromatin lines)
Ropen=1; %radius of the opening operator for cleaning the mask
ratio1=16;
ratio2=32;
ratio3=60;
ratio4=120;


%%%User infos
dial=inputdlg({'First time-point after 405'},'Stack info',1,{'1'});
Firsttime=str2num(dial{1});
resolution=inputdlg({'x,y resolution (um)','time resolution (s)'},'Resolution',1,{'0.108','4'});
PixRes=str2num(resolution{1});
TimeRes=str2num(resolution{2});
ratios=inputdlg({'ratio1','ratio2','ratio3','ratio4'},'Ratios',1,{'16','32','60','120'});
ratio1=round(str2num(ratios{1})/TimeRes);
ratio2=round(str2num(ratios{2})/TimeRes);
ratio3=round(str2num(ratios{3})/TimeRes);
ratio4=round(str2num(ratios{4})/TimeRes);

%%%Stacks selection

[fnstack fpstack]=uigetfile('*.TIF','Select the stacks you want to analyze','MultiSelect','on');
Nbstack=size(fnstack,2);

[fntxt fptxt]=uiputfile('*.txt','Create an text file for the raw output');
[fntxt2 fptxt2]=uiputfile('*.txt','Create an text file for the ratio output');
[fntxt3 fptxt3]=uiputfile('*.txt','Create an text file for the fit output');

fidtext=fopen([fptxt fntxt],'w');
tic
Ratiovect=zeros(4,Nbstack);
%fitres=zeros(6,Nbstack);
fitres=zeros(2,Nbstack);
for istack=1:Nbstack
    disp(['analysed stack : ' fnstack{istack}])
    fprintf(fidtext, '%s\n',fnstack{istack});
    %%%Estimating the number of time points to analyze
    timepoint=Firsttime;
    I=0;
    while I~=1
        I=isempty(double(tiffread([fpstack fnstack{istack}],timepoint)));
        timepoint=timepoint+1;
    end
    Lasttime=timepoint-2;
    
    charactlines=zeros(Lasttime-Firsttime+1,2);
    charactlines(:,1)=Firsttime:Lasttime;
    ntime=1;
    for timepoint=Firsttime:Lasttime
        disp(['time point = ' num2str(timepoint)]);
        %%%Image loading and segmentation by kmeans
        I=double(tiffread([fpstack fnstack{istack}],timepoint));
        [mu,mask]=kmeans_image(I,nbclust);
        mask=(mask==signalclust);
        
        %%%Mask cleaning
        %figure
        %imagesc(mask)
        mask=imclose(mask,strel('disk',Rclose));
        mask=imclose(mask,strel('line',Lengthclose,90));
        mask=imfill(mask,'holes');
        mask=imopen(mask,strel('disk',Ropen));
        [chromlines,volchromlines]=findobject(mask,6); 

        if length(volchromlines)>1 %%%if more than 1 detected objects only keep the biggest
            volsort=sort(volchromlines);
            line1=find(volchromlines==volsort(length(volsort)),1);
            for j=1:length(volchromlines)
                if all(j~=line1)
                    mask(chromlines{j})=0;
                end
            end
        end
        %figure
        %imagesc(mask)
        %imwrite(mask,[fp fn 't' sprintf('%d',timepoint) '-segment.tif'],'tif','Compression','none');
        [chromlines,numchromlines]=bwlabeln(mask);
        STATS=regionprops(chromlines,'Centroid','MajorAxisLength','MinorAxisLength','Orientation');
        charactlines(ntime,2)=[STATS.MinorAxisLength]*PixRes;
        ntime=ntime+1;
    end
    Ratiovect(1,istack)=charactlines(ratio1+1,2)/charactlines(1,2);
    Ratiovect(2,istack)=charactlines(ratio2+1,2)/charactlines(1,2);
    Ratiovect(3,istack)=charactlines(ratio3+1,2)/charactlines(1,2);
    Ratiovect(4,istack)=charactlines(ratio4+1,2)/charactlines(1,2);
    
    normline=charactlines(:,2)/charactlines(1,2);
    timevect=(charactlines(:,1)-Firsttime)*TimeRes;
    %%% Single exponential fit
    initparam=zeros(2,1);
    initparam(1)=normline(length(normline))-1;
    initparam(2)=timevect(find(abs(normline-1-initparam(1)/2)==min(abs(normline-1-initparam(1)/2)),1));
    Fittedparam1 = lsqcurvefit(@fitfun1,initparam,timevect,normline);
    fitres(1,istack)=Fittedparam1(1);
    fitres(2,istack)=Fittedparam1(2);
    %%% Double exponential fit
    %initparam=zeros(4,1);
    %initparam(1)=normline(length(normline))-1;
    %initparam(2)=timevect(find(abs(normline-1-initparam(1)/4)==min(abs(normline-1-initparam(1)/4)),1));
    %initparam(3)=0.5;
    %initparam(4)=timevect(find(abs(normline-1-3*initparam(1)/4)==min(abs(normline-1-3*initparam(1)/4)),1));
    %Fittedparam2 = lsqcurvefit(@fitfun2,initparam,timevect,normline);
    %fitres(3,istack)=Fittedparam2(1);
    %fitres(4,istack)=Fittedparam2(2);
    %fitres(5,istack)=Fittedparam2(3);
    %fitres(6,istack)=Fittedparam2(4);

    figure
    set(gcf,'Name',fnstack{istack})
    fitvect1=1+Fittedparam1(1)*(1-exp(-timevect/Fittedparam1(2)));
    %fitvect2=1+Fittedparam2(1)*(1-Fittedparam2(3)*exp(-timevect/Fittedparam2(2))-(1-Fittedparam2(3))*exp(-timevect/Fittedparam2(4)));
    %plot(timevect,normline,timevect,fitvect1,timevect,fitvect2)
    plot(timevect,normline,timevect,fitvect1)
    %plot(timevect,normline)
    title('average thickness of the chromatin lines')
    fprintf(fidtext, '%f %f\n',charactlines');    
end
toc
fclose(fidtext);
csvwrite([fptxt2 fntxt2],Ratiovect);
csvwrite([fptxt3 fntxt3],fitres);


%%%Fitting functions
function yfit = fitfun1(params,T)
A1=params(1);
Tau1=params(2);
yfit=1+A1*(1-exp(-T/Tau1)); 

function yfit = fitfun2(params,T)
A1=params(1);
Tau1=params(2);
A2=params(3);
Tau2=params(4);
yfit=1+A1*(1-A2*exp(-T/Tau1)-(1-A2)*exp(-T/Tau2)); 
