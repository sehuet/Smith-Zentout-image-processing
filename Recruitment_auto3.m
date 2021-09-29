function Recruitment_auto3

%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Author: Sébastien Huet
%   Email: sebastien.huet (at) univ-rennes1.fr
%   Date: 11-05-2021
%   
%   This is the main function used to assess recruitment kinectics of PARP1 at sites of damage
%
%%%%%%%%%%%%%%%%%%%%%%%%%


%%% Segmentation parameters for the cell finding
maxsearch=4; % search for cells till number of planes N/maxsearch
nbclust_cell=2;      % number of clusters for the k-means segmentation
signalclust_cell=2; %index of the cluster containing the signal of interest
Lengthopen=3; %line opening to keep only the line at the time of UV irrad
minlinesize=1000; %min size of the object detected as an irrad line
nbstd=4; %threshold for nuc detection is set at the mean in the irrad ROI - nbstd * std in the irrad ROI
RnucREC=5; % radius of the opening operator for cleaning the mask
minnucsize=1000; %min size of the object detected as nucleus
ROIborder=25; %nb of pix added around the ROI defined for the analyzed nuc

%%% Segmentation parameters for the H2B photo-activated line
nbclust=2;      % number of clusters for the first k-means segmentation
minsize=1000;    % minimum size of the object kept by bwareapoen for cleaning the mask
Ropen=4;        % radius of the opening operator for cleaning the mask
Lengthclose=9;  % length of the closing line for cleaning the mask (reassociate broken lines)

%%% Segmentation parameters for the entire nucleus (using the GFP channel)
Rnuc=2;         % radius of the closing and opening operators for cleaning the mask
BGREClimit=150; %max value of the BG in the recruitment channel over which the manual value is used

%%%Fitting of the recruitment curves
maxiter=1000; %max nb of iteration for the fit
nbparamfit=4;
Nbpeak=10; %%%nb of max intensities values used to estimate the plateau time tau_0

%%% User information
imageinfos=inputdlg({'x,y resolution (µm)','time resolution (s)','Tag for the recruitment images','Manual BG H2B','Manual BG REC','Max nb of cells in fov','Max nb of planes in movie'},'Image info',1,{'0.1075','4','GFP','110','110','4','40'});
PixRes=str2num(imageinfos{1});      
TimeRes=str2num(imageinfos{2});
RECTag=imageinfos{3};
BGH2Bmanual=str2num(imageinfos{4});
BGRECmanual=str2num(imageinfos{5});
maxnbcells=str2num(imageinfos{6});
maxnbplanes=str2num(imageinfos{7});
choiceseg=questdlg('Segmentation method for the REC channel','User info','auto','manu','auto');

%%% Stacks selection
[fnstack fpstack]=uigetfile('*.TIF','Select stacks to analyze (chromatin marker channel)','MultiSelect','on');
Nbstack=size(fnstack,2)
[fntxt fptxt]=uiputfile('*.txt','Create a text file for the raw output');
fncommontxt=fntxt(1:(find(fntxt=='.',1,'last')-1));
fidtext=fopen([fptxt fncommontxt '-info.txt'],'w'); 

cellindex=1;
Results_int=nan(maxnbplanes,8*Nbstack*maxnbcells);
Results_mean=nan(maxnbplanes,8*Nbstack*maxnbcells);
Results_fit=nan(Nbstack*maxnbcells,nbparamfit);
for istack=1:Nbstack
    disp(['analysed stack : ' fnstack{istack}])
    [Lim,Cim]=size(double(tiffread([fpstack fnstack{istack}],1)));
    %%% Estimating the number of planes in the stack
    Nbplanes=1;
    I=0;
    while I~=1
        I=isempty(double(tiffread([fpstack fnstack{istack}],Nbplanes)));
        Nbplanes=Nbplanes+1;
    end
    Nbplanes=Nbplanes-2;
    
    %%% Loading of the name of the REC images
    fnstackH2B=fnstack{istack};
    fncommoni=fnstackH2B(1:find(fnstackH2B==' ',1,'last'));
    fnstackREC=[fncommoni RECTag '.tif'];
    fprintf(fidtext, '%s\n ',regexprep(fncommoni,' ','_')); %%%write the name of the movie in the info text file
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%identifying the indididual cells and the time of irrad along the whole
    %%%stack
    if all(choiceseg=='manu')
        maskRECmanu=zeros(Lim,Cim,maxnbcells);
    end
    timepoint=1;
    nbcells=1;
    timeirrad=zeros(maxnbcells,1);
    ROInuc=zeros(maxnbcells,4);
    while (timepoint<=Nbplanes/maxsearch & nbcells<=maxnbcells)
        timepoint
        IH2Bdiff=double(tiffread([fpstack fnstackH2B],timepoint+1))-double(tiffread([fpstack fnstackH2B],timepoint));
        IH2Bdiff=IH2Bdiff.*(IH2Bdiff>0); % only keep positive values (we look for an increase)
        
        [mu,maskH2B]=kmeans_image(IH2Bdiff,nbclust_cell);
        maskH2B=(maskH2B==signalclust_cell);
        %maskH2B=imopen(maskH2B,strel('line',Lengthopen,90));
       
        maskH2B=bwareaopen(maskH2B,minlinesize);
        %figure
        %imagesc(maskH2B)
        
        %%%ROI for the individual cells and times of irrad
        if any(maskH2B(:)>0)
            disp(['nb of detected cells = ' num2str(nbcells)]);
            timeirrad(nbcells)=timepoint+1;
            
            if all(choiceseg=='auto')
                IREC=double(tiffread([fpstack fnstackREC],timepoint));
                IREC_inmask=IREC.*maskH2B;
                IREC_thresh_low=mean(IREC_inmask(IREC_inmask(:)>0))-nbstd*std(IREC_inmask(IREC_inmask(:)>0));
                IREC_thresh_high=mean(IREC_inmask(IREC_inmask(:)>0))+nbstd*std(IREC_inmask(IREC_inmask(:)>0));
                maskREC=(IREC>IREC_thresh_low & IREC<IREC_thresh_high);
                maskREC=imfill(maskREC,'holes');
                maskREC=imopen(maskREC,strel('disk',RnucREC));      
                maskREC=bwareaopen(maskREC,minnucsize);
                [maskREC,numnuc] = bwlabel(maskREC);
                overlaparea=zeros(numnuc,1);
                for nuc=1:numnuc
                    maskRECnuc=(maskREC==nuc);
                    overlaparea(nuc)=sum(sum(maskRECnuc.*maskH2B));
                end
                nuc_select=find(overlaparea==max(overlaparea));
                maskREC=(maskREC==nuc_select);
                
            else
                IREC=zeros(Lim,Cim,3);
                IREC(:,:,1)=maskH2B;
                temp=double(tiffread([fpstack fnstackREC],timepoint));
                meanI=sum(sum(temp.*maskH2B))/sum(maskH2B(:));
                IREC(:,:,2)=temp/meanI*0.75;
                h=figure;
                maskRECmanu(:,:,nbcells)=roipoly(IREC);
                close(h)
                maskREC=maskRECmanu(:,:,nbcells);
            end
            STATS=regionprops(double(maskREC),'BoundingBox','Centroid');
            ROInuc(nbcells,:)=(round([STATS.BoundingBox]));
            ROInuc(nbcells,:)=[max([ROInuc(nbcells,2)-ROIborder,1]),min([ROInuc(nbcells,2)+ROInuc(nbcells,4)-1+ROIborder,size(maskREC,1)]),max([ROInuc(nbcells,1)-ROIborder,1]),min([ROInuc(nbcells,1)+ROInuc(nbcells,3)-1+ROIborder,size(maskREC,2)])];
            
            %%info written in the output info text file to be able to back track the cell
            PosCell=round([STATS.Centroid]); 
            fprintf(fidtext, '%f %f %f \n',[timeirrad(nbcells) PosCell(1) PosCell(2)]);
            masksum=maskH2B+maskREC;
            %figure
            %imagesc(masksum(ROInuc(nbcells,1):ROInuc(nbcells,2),ROInuc(nbcells,3):ROInuc(nbcells,4)))
            nbcells=nbcells+1;
        end
        timepoint=timepoint+1;
    end
    nbcellstot=nbcells-1;
    timeirrad
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%Creation of the loop to analyze all timepoints after photo-perturbation
    for nbcells=1:nbcellstot
        
        Firsttime=timeirrad(nbcells);        
        ntime=1;
        
        for timepoint=Firsttime:Nbplanes
            
            disp(['time point = ' num2str(timepoint)]);
            
            %%%%%%% Segmentation of the chromatin line (H2B channel) %%%%%%%
            %%% Image loading and segmentation of the photo-irradiated ROI
            IH2B=double(tiffread([fpstack fnstack{istack}],timepoint));
            IH2B=IH2B(ROInuc(nbcells,1):ROInuc(nbcells,2),ROInuc(nbcells,3):ROInuc(nbcells,4));
            %IH2B=medfilt2(IH2B);
            [mu,maskH2B]=kmeans_image(IH2B,nbclust);
            maskH2B=(medfilt2(IH2B)>=(mean(mu)));
            
            %%%Mask cleaning
            %%maskH2B=imclose(maskH2B,strel('disk',Rclose));
            maskH2B=bwareaopen(maskH2B,minlinesize);
            maskH2B=imfill(maskH2B,'holes');
            maskH2B=imopen(maskH2B,strel('disk',Ropen));
            maskH2B=imclose(maskH2B,strel('line',Lengthclose,90));
            %figure
            %imagesc(maskH2B)   
            [chromlines,volchromlines]=findobject(maskH2B,6);
            
            %%% Keeping only the biggest segmented object
            if length(volchromlines)>1
                volsort=sort(volchromlines);
                line1=find(volchromlines==volsort(length(volsort)),1);
                for j=1:length(volchromlines)
                    if all(j~=line1)
                        maskH2B(chromlines{j})=0;
                    end
                end
            end
            
            %%% Fitting the ellipse on the mask of the photo-converted line and extracting the thickness
            if sum(maskH2B(:))~=0
                STATS=regionprops(double(maskH2B),'Centroid','MajorAxisLength','MinorAxisLength','Orientation');
                thickness=[STATS.MinorAxisLength];
                PosCell=[STATS.Centroid]; 
            end
            
            %%%%%%% Segmentation of the whole nucleus (REC channel) %%%%%%%
            %%% The threshold is estimated based on the images before irradiation
            if timepoint==Firsttime
                PosCell_1=PosCell;
                maskH2B_1=maskH2B;
                if all(choiceseg=='auto')
                    thresh_low=zeros(Firsttime-1,1);
                    thresh_high=zeros(Firsttime-1,1);
                    for i=1:Firsttime-1
                        IREC=double(tiffread([fpstack fnstackREC],i));
                        IREC=IREC(ROInuc(nbcells,1):ROInuc(nbcells,2),ROInuc(nbcells,3):ROInuc(nbcells,4));
                        IREC_inmask=IREC.*maskH2B;
                        thresh_low(i)=mean(IREC_inmask(IREC_inmask(:)>0))-nbstd*std(IREC_inmask(IREC_inmask(:)>0));
                        thresh_high(i)=mean(IREC_inmask(IREC_inmask(:)>0))+nbstd*std(IREC_inmask(IREC_inmask(:)>0));
                    end
                    thresh_low=sum(thresh_low)/(Firsttime-1);
                    thresh_high=sum(thresh_high)/(Firsttime-1);
                    IREC=double(tiffread([fpstack fnstackREC],Firsttime-1));
                    IREC=IREC(ROInuc(nbcells,1):ROInuc(nbcells,2),ROInuc(nbcells,3):ROInuc(nbcells,4));
                    masknuc=(IREC>thresh_low & IREC<thresh_high);
                
                    %%% Cleaning of the mask of the nucleus
                    masknuc=imfill(masknuc,'holes');
                    masknuc=imopen(masknuc,strel('disk',Rnuc));   
                    masknuc=imclose(masknuc,strel('disk',Rnuc));             
                    [nucreg,volnucreg]=findobject(masknuc,6);
                
                    %%% Keeping only the biggest segmented object
                    if length(volnucreg)>1
                        volsort=sort(volnucreg);
                        obj1=find(volnucreg==volsort(length(volsort)),1);
                        for j=1:length(volnucreg)
                            if all(j~=obj1)
                                masknuc(nucreg{j})=0;
                            end
                        end
                    end
                else
                    masknuc=maskRECmanu(ROInuc(nbcells,1):ROInuc(nbcells,2),ROInuc(nbcells,3):ROInuc(nbcells,4),nbcells);
                end
                masknuc_1=masknuc;
                
                %%% Measurements of the backgrounds in the 2 channels and calculation of pre-photo-perturbation intensities for normalization
                BGH2B=zeros(Firsttime-1,1);
                BGREC=zeros(Firsttime-1,1);
                RECnorm=zeros(Firsttime-1,1);
                for i=1:Firsttime-1
                    IH2B=double(tiffread([fpstack fnstackH2B],i));
                    IH2B=IH2B(ROInuc(nbcells,1):ROInuc(nbcells,2),ROInuc(nbcells,3):ROInuc(nbcells,4));
                    IREC=double(tiffread([fpstack fnstackREC],i));
                    IREC=IREC(ROInuc(nbcells,1):ROInuc(nbcells,2),ROInuc(nbcells,3):ROInuc(nbcells,4));
                    BGH2B(i)=sum(sum(IH2B.*masknuc))/sum(sum(masknuc));
                    BGREC(i)=sum(sum(IREC.*~masknuc))/sum(sum(~masknuc));
                    RECnorm(i)=sum(sum((IREC-BGREC(i)).*maskH2B))/sum(sum((IREC-BGREC(i)).*masknuc));
                end
                BGH2B=sum(BGH2B)/(Firsttime-1);
                BGREC=sum(BGREC)/(Firsttime-1);
                
                
                %%%if nucleus segmentation is crap, use manually defined
                %%%values for the BG intensities
                if BGREC>BGREClimit
                    disp('Background determination is crap. Fixed value is used')
                    BGH2B=BGH2Bmanual;
                    BGREC=BGRECmanual;
                    for i=1:Firsttime-1            
                        IREC=double(tiffread([fpstack fnstackREC],i));
                        IREC=IREC(ROInuc(nbcells,1):ROInuc(nbcells,2),ROInuc(nbcells,3):ROInuc(nbcells,4));
                        RECnorm(i)=sum(sum((IREC-BGREC).*maskH2B))/sum(sum((IREC-BGREC).*masknuc));
                    end
                end
                RECnorm=sum(RECnorm)/(Firsttime-1);
                RECnorm_mean=RECnorm*sum(sum(masknuc))/sum(sum(maskH2B));
                
                %%% Normalized intensity pre-damage considered equal to 1
                
                %%% Integrated intensities
                Results_int(ntime,8*(cellindex-1)+1)=Firsttime-1;
                Results_int(ntime,8*(cellindex-1)+7)=1;
                %%% Mean intensities
                Results_mean(ntime,8*(cellindex-1)+1)=Firsttime-1;
                Results_mean(ntime,8*(cellindex-1)+7)=1;
                ntime=ntime+1;
                          
            end
            
            %%%Nuc registration (shifts the mask of the nucleus based
            %%%on the movement of the irradiated area (segmented at
            %%%each frame)
            ShiftNuc=round(PosCell-PosCell_1);
            masknuc=circshift(masknuc_1,ShiftNuc(2),1);
            masknuc=circshift(masknuc,ShiftNuc(1),2);
            
            
            
            %%% Segmentation Results
            if timepoint==Nbplanes

                IH2B_1=double(tiffread([fpstack fnstackH2B],Firsttime-1));
                IH2B_1=IH2B_1(ROInuc(nbcells,1):ROInuc(nbcells,2),ROInuc(nbcells,3):ROInuc(nbcells,4));
                IREC_1=double(tiffread([fpstack fnstackREC],Firsttime-1));
                IREC_1=IREC_1(ROInuc(nbcells,1):ROInuc(nbcells,2),ROInuc(nbcells,3):ROInuc(nbcells,4));
                IH2B_2=double(tiffread([fpstack fnstackH2B],Firsttime));
                IH2B_2=IH2B_2(ROInuc(nbcells,1):ROInuc(nbcells,2),ROInuc(nbcells,3):ROInuc(nbcells,4));
                IREC_2=double(tiffread([fpstack fnstackREC],Firsttime));
                IREC_2=IREC_2(ROInuc(nbcells,1):ROInuc(nbcells,2),ROInuc(nbcells,3):ROInuc(nbcells,4));
                IH2B_3=double(tiffread([fpstack fnstackH2B],Nbplanes));
                IH2B_3=IH2B_3(ROInuc(nbcells,1):ROInuc(nbcells,2),ROInuc(nbcells,3):ROInuc(nbcells,4));
                IREC_3=double(tiffread([fpstack fnstackREC],Nbplanes));
                IREC_3=IREC_3(ROInuc(nbcells,1):ROInuc(nbcells,2),ROInuc(nbcells,3):ROInuc(nbcells,4));
                
                [Lim,Cim]=size(IH2B_1);
                overlay1=zeros(Lim,Cim,3);
                overlay1_2=zeros(Lim,Cim,3);
                overlay2=zeros(Lim,Cim,3);
                overlay2_2=zeros(Lim,Cim,3);
                overlay3=zeros(Lim,Cim,3);
                overlay3_2=zeros(Lim,Cim,3);
                
                overlay1(:,:,1)=bwperim(maskH2B_1);
                overlay1(:,:,2)=(IH2B_1-min(IH2B_1(:)))/(max(IH2B_1(:))-min(IH2B_1(:)));
                overlay1(:,:,3)=bwperim(masknuc_1);
                overlay1_2(:,:,1)=bwperim(maskH2B_1);
                overlay1_2(:,:,2)=(IREC_1-min(IREC_1(:)))/(max(IREC_1(:))-min(IREC_1(:)));
                overlay1_2(:,:,3)=bwperim(masknuc_1);
                overlay2(:,:,1)=bwperim(maskH2B_1);
                overlay2(:,:,2)=(IH2B_2-min(IH2B_2(:)))/(max(IH2B_2(:))-min(IH2B_2(:)));
                overlay2(:,:,3)=bwperim(masknuc_1);
                overlay2_2(:,:,1)=bwperim(maskH2B_1);
                overlay2_2(:,:,2)=(IREC_2-min(IREC_2(:)))/(max(IREC_2(:))-min(IREC_2(:)));
                overlay2_2(:,:,3)=bwperim(masknuc_1);
                overlay3(:,:,1)=bwperim(maskH2B);
                overlay3(:,:,2)=(IH2B_3-min(IH2B_3(:)))/(max(IH2B_3(:))-min(IH2B_3(:)));
                overlay3(:,:,3)=bwperim(masknuc);
                overlay3_2(:,:,1)=bwperim(maskH2B);
                overlay3_2(:,:,2)=(IREC_3-min(IREC_3(:)))/(max(IREC_3(:))-min(IREC_3(:)));
                overlay3_2(:,:,3)=bwperim(masknuc);
                
                figure;
                set(gcf,'Position',[50 50 1000 600])
                subplot(2,3,1)
                imagesc(overlay1)
                title('Last frame before photo-perturbation')
                ylabel('H2B channel')
                axis equal
                subplot(2,3,2)
                imagesc(overlay2)
                title('First frame after photo-perturbation')
                axis equal
                subplot(2,3,3)
                imagesc(overlay3)
                title('Last frame of the movie')
                axis equal
                subplot(2,3,4)
                imagesc(overlay1_2)
                ylabel('REC channel')
                axis equal
                subplot(2,3,5)
                imagesc(overlay2_2)
                title(fncommoni,'Interpreter','none')
                axis equal
                subplot(2,3,6)
                imagesc(overlay3_2)
                axis equal
            end
            
            %%% Intensity Measurements
            %%%Integrated intensities
            IREC=double(tiffread([fpstack fnstackREC],timepoint));
            IREC=IREC(ROInuc(nbcells,1):ROInuc(nbcells,2),ROInuc(nbcells,3):ROInuc(nbcells,4));
            IH2B=double(tiffread([fpstack fnstackH2B],timepoint));
            IH2B=IH2B(ROInuc(nbcells,1):ROInuc(nbcells,2),ROInuc(nbcells,3):ROInuc(nbcells,4));
            IbreakREC=sum(sum((IREC-BGREC).*maskH2B));
            IbreakH2B=sum(sum((IH2B-BGH2B).*maskH2B));
            InucREC=sum(sum((IREC-BGREC).*masknuc));
            InucH2B=sum(sum((IH2B-BGH2B).*masknuc));
            %%%Mean intensities
            IbreakREC_mean=IbreakREC/sum(sum(maskH2B));
            InucREC_mean=InucREC/sum(sum(masknuc));
            IbreakH2B_mean=IbreakH2B/sum(sum(maskH2B));
            InucH2B_mean=InucH2B/sum(sum(masknuc));
            
            %%% Data collection in the 'Results' matrices
            %%% Integrated intensities
            Results_int(ntime,8*(cellindex-1)+1)=timepoint;
            Results_int(ntime,8*(cellindex-1)+2)=IbreakREC;
            Results_int(ntime,8*(cellindex-1)+3)=InucREC;
            Results_int(ntime,8*(cellindex-1)+4)=IbreakH2B;
            Results_int(ntime,8*(cellindex-1)+5)=InucH2B;
            Results_int(ntime,8*(cellindex-1)+6)=thickness*PixRes;
            Results_int(ntime,8*(cellindex-1)+7)=(IbreakREC/InucREC)/RECnorm;
            Results_int(ntime,8*(cellindex-1)+8)=IbreakH2B/InucH2B;
            %%% Mean intensities
            Results_mean(ntime,8*(cellindex-1)+1)=timepoint;
            Results_mean(ntime,8*(cellindex-1)+2)=IbreakREC_mean;
            Results_mean(ntime,8*(cellindex-1)+3)=InucREC_mean;
            Results_mean(ntime,8*(cellindex-1)+4)=IbreakH2B_mean;
            Results_mean(ntime,8*(cellindex-1)+5)=InucH2B_mean;
            Results_mean(ntime,8*(cellindex-1)+6)=thickness*PixRes;
            Results_mean(ntime,8*(cellindex-1)+7)=(IbreakREC_mean/InucREC_mean)/RECnorm_mean;
            Results_mean(ntime,8*(cellindex-1)+8)=IbreakH2B_mean/InucH2B_mean;
            
            ntime=ntime+1;
        end
        
        
        %%%%Interpolate with spline prior to recruit peak to balance the
        %%%%number of points before and after the recruit peak
        fitpoints=find(~isnan(Results_mean(:,8*(cellindex-1)+1)));   
        time=(Results_mean(fitpoints,8*(cellindex-1)+1)-Firsttime+1)*TimeRes;
        I_exp=Results_mean(fitpoints,8*(cellindex-1)+7);
        
        h=figure;
        plot(time,I_exp,'b')
        
        
        options.Interpreter = 'tex';
        options.Default = 'Yes';
        choice = questdlg('Would you like to fit this recruitment curve?','Continue', 'Yes','No',options);
        close(h)
        switch choice
        case 'Yes'
            nbpts=length(I_exp);
            idxmax=find(I_exp==max(I_exp));

            trec1=time(1:idxmax);
            trec2=(time(1):(TimeRes*(nbpts/idxmax).^(-1)):time(idxmax))';
            I_exp_interp=spline(trec1,I_exp(1:idxmax),trec2);
            I_exp_interp=[I_exp_interp ; I_exp(idxmax+1:nbpts)];
            time_interp=[trec2 ; time(idxmax+1:nbpts)];

            %%%Fit of the interpolated recruitment curves with the CRC model%%%

            %%%setting the initial guesses
            idxmax=find(I_exp_interp==max(I_exp_interp));        
            A_0=I_exp_interp(idxmax)-1;
            halfrec=abs(I_exp_interp-(A_0/2+1));
            idx1=find(halfrec(1:idxmax)==min(halfrec(1:idxmax)));
            idx2=idxmax+find(halfrec((idxmax+1):length(halfrec))==min(halfrec((idxmax+1):length(halfrec)))); 
            k1_0=1/time_interp(idx1);
            k2_0=1/(time_interp(idx2)-time_interp(idxmax));

            temp=sortrows([I_exp_interp time_interp],-1);
            temp=temp(1:Nbpeak,2);
            tau_0=max(temp)-min(temp);

            %%%fitting the recruitment curve
            fit_fun= @(PAR, x) CRC1(PAR,x);
            COEF_0=[A_0,tau_0,k1_0,k2_0]
            [COEF,res,J] = nlinfit(time_interp,I_exp_interp,fit_fun, COEF_0);
            Results_fit(cellindex,:)=[COEF(1),COEF(2),1/COEF(3),1/COEF(4)];
            %%%plotting the fit result
            recfit_interp = CRC1(COEF,time_interp);
            figure
            plot(time,I_exp,'b',time_interp,I_exp_interp,'b.',time_interp,recfit_interp,'r')
            legend('Exp','Interp','Fit')
            title(fncommoni,'Interpreter','none')
            axis tight
        case 'No'
            Results_fit(cellindex,:)=[NaN,NaN,NaN,NaN];
        end
        
        cellindex=cellindex+1;
    end
end
    
fclose(fidtext);
csvwrite([fptxt fncommontxt '-rawint.txt'],Results_int);
csvwrite([fptxt fncommontxt '-rawmean.txt'],Results_mean);
csvwrite([fptxt fncommontxt '-thick.txt'],Results_int(:,6:8:size(Results_int,2)));
csvwrite([fptxt fncommontxt '-normint.txt'],Results_int(:,7:8:size(Results_int,2)));
csvwrite([fptxt fncommontxt '-fitres.txt'],Results_fit);