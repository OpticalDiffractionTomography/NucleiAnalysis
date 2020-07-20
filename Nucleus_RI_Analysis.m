clear all
close all
%% Parameters
% input variables:
% Reconimg    : RI tomogram (x,y,z)
% res3        : Lateral resolution [um]
% res4        : Axial resolution [um]
% lambda      : Illumination wavelength [um]
% fluoImg     : colocalized epi-fluorescence image (Green and Red channels
%               for FUCCI marker (fluoImg(:,:,1) is Green, fluoImg(:,:,2)
%               is red channel fluorescence.
RII=0.19;       % RI increment [mL/g] = [fL/pg]
n_m=1.337;      % RI of medium
n_s=n_m+0.04;   % maximum RI value for images (only visualization pulpose)
%% Matching resolution of fluorescence and RI tomogram
% measured epi-fluorescence image is oversampled.
% By cropping fluorescence images in the Fourier space, the spatial
% resolution of fluorescence images and RI tomogram are matched.
combinedAnalysis=zeros(size(Reconimg,1),size(Reconimg,2),3);
for mm=1:2
    temp=fluoImg(:,:,mm);
    temp=temp-mean2(temp(end-15:end-5,5:15));
    temp=conv2(temp,fspecial('disk',2),'same');
    pmap=imgaussfilt(temp,300);                 % Flattening images
    temp=temp-pmap;
    temp = wiener2(temp,[5 5]);                 % Reducing noise
    ii=length(temp);
    Ftemp=fftshift(fft2(temp./max(max(temp))*255))/(ii^2);  % 2D Fourier cropping
    Ftemp=Ftemp(end/2-size(Reconimg,1)/2+1:end/2+size(Reconimg,1)/2,end/2-size(Reconimg,1)/2+1:end/2+size(Reconimg,1)/2);
    sizeFimg=length(Ftemp);
    Ftemp=ifft2(ifftshift(Ftemp))*(sizeFimg^2);
    combinedAnalysis(:,:,3-mm)=abs(Ftemp);
end
% Reducing noise in the tomogram
for zz=1:size(Reconimg,3)
    temp=conv2(squeeze(Reconimg(:,:,zz)),fspecial('disk',0.7),'same');
    temp=temp(3:end-2,3:end-2);
    temp=padarray(temp,[2 2],n_m);
    Reconimg(:,:,zz)=temp;
end
%% Cell Segmentation #1
% Applying the Otsu's method to segment cell
% Output: 3D masks of individual cell
level = multithresh(Reconimg,1);
cellMap_temp=((Reconimg)>level);
% connect neighboring fractions
se=strel('sphere',4);
cellMap_temp=imfill(cellMap_temp,'holes');
cellMap_temp=imdilate(cellMap_temp,se);
cellMap_temp=imfill(cellMap_temp,'holes');
se=strel('sphere',4);
cellMap_temp=imerode(cellMap_temp,se);
cellMap_temp=bwareaopen(cellMap_temp,500);
% Labeling individual cell by the watershed algorithm
D = bwdist(~cellMap_temp);
D = -D;
D(~cellMap_temp) = -Inf;
cellMap_seg1 = watershed(imhmin(D,3));
cellMap_seg1(~cellMap_temp) = 0;
se=strel('sphere',1);
cellMap_seg1=single(imdilate(cellMap_seg1,se));
%% Nuclei Segmentation #1
%   Nuclei are segmented from epi-fluorescence images of FUCCI-stained cells
%   Output: 2D mask of nuclei in the FOV
nucMap_seg1=zeros(size(combinedAnalysis,1),size(combinedAnalysis,2));
for cellNum=1:max(max(max(cellMap_seg1)))
    cellMap=(cellMap_seg1==(cellNum));
    cellMap=imfill(cellMap,'holes');
    try
        %%
        z=-1;           % The axial position at the epi-fluorescence is lying at the center of the RI tomogram.
        nucMap2=ones(size(combinedAnalysis,1),size(combinedAnalysis,2));
        tempVal=max(combinedAnalysis,[],3).*cellMap(:,:,end/2+z);
        level = multithresh(tempVal(tempVal>0),1);
        nucMap_temp=(max(combinedAnalysis,[],3).*cellMap(:,:,end/2+z))>level;
        % segment nuclei in both channels
        if level>median(median(max(combinedAnalysis,[],3)))+10
            nucMap_temp=imfill(nucMap_temp,'holes');
            nucMap_temp=bwareafilt(nucMap_temp,[500 100000]);
            nucMap2=nucMap2.*~nucMap_temp;
        end
        nucMap2=1-nucMap2;
        % connect neighboring fractions
        se=strel('disk',2);
        nucMap_temp=imfill(nucMap2,'holes');
        nucMap_temp=imdilate(nucMap_temp,se);
        nucMap_temp=imfill(nucMap_temp,'holes');
        nucMap_temp=imerode(nucMap_temp,se);
        nucMap_seg1=nucMap_seg1+nucMap_temp;
    catch
    end
end
% Labeling individual nuclei by the watershed algorithm
D = bwdist(~nucMap_seg1);
D = -D;
D(~nucMap_seg1) = -Inf;
nucMap_seg2 = watershed(imhmin(D,3));
nucMap_seg2(~nucMap_seg1) = 0;
se=strel('sphere',1);
nucMap_seg2=single(imdilate(nucMap_seg2,se));
%% Cell Segmentation #2
% weighting nuclear regions labeled in the previous step,
% Nuclei Segmentation #1, to the 3D masks of individual cell
cellMap=cellMap_seg1>0;
cellMap=bwareaopen(cellMap,500);
nucMap=repmat(nucMap_seg2,[1 1 size(cellMap,3)]);
cellMap=cellMap+nucMap;
D = bwdist(~cellMap);
D = -D;
D(~cellMap) = -Inf;
cellMap4 = watershed(imhmin(D,3));
cellMap4(~cellMap) = 0;
se=strel('sphere',1);
cellMap4=single(imdilate(cellMap4,se));
cellMap4(cellMap_seg1==0)=0;
cellMap_seg2=cellMap4;
%% Quantitative Characterization
qData=[];
tic
for cellNum=1:max(max(max(cellMap_seg2)))
    cellMap_selected=(cellMap_seg2==(cellNum));          % Choose labeled cells
    cellMap_selected=imfill(cellMap_selected,'holes');
    [tempX tempY tempZ]=ind2sub(size(cellMap_selected),find(cellMap_selected));
    qDataTemp=zeros(1,20);
    try
        %%
        % Erosion of cell mask along the axial direction to take into
        % account the lower spatial resolution along the axial direction
        kk=1;
        [xxx yyy zzz]=meshgrid(-15:15,-15:15,-15:15);
        se=real(single((xxx.^2/kk.^2+yyy.^2/kk.^2+zzz.^2/(kk/(res3/res4)+1).^2)<1^2));
        cellMap_selected=imerode(cellMap_selected,se);
        Reconimg_selected=Reconimg.*cellMap_selected;
        %% Surface Area Calculation
        % Summation of areas of triangles in the mesh of isosurface of cells
        p2 = isosurface(cellMap_selected,0.8);
        v = p2.vertices;
        v(:,1:2)=v(:,1:2)*res3;
        v(:,3)=v(:,3)*res4;
        f = p2.faces;
        a = v(f(:,2),:) - v(f(:,1),:);
        b = v(f(:,3),:) - v(f(:,1),:);
        c = cross(a,b,2);
        area = 0.5*sum(sqrt(sum(c.^2,2)));
        %% Nucleus Segmentation #2
        % segment nucleus in the cell selected in 3D
        nucleusMap=ones(size(combinedAnalysis,1),size(combinedAnalysis,2));
        for mm=1:2
            tempVal=combinedAnalysis(:,:,3-mm).*cellMap_selected(:,:,end/2+z);
            level = multithresh(tempVal(tempVal>0),1);
            nucMap_temp=((combinedAnalysis(:,:,3-mm).*cellMap_selected(:,:,end/2+z))>level(1));
            
            if level>median(median(combinedAnalysis(:,:,3-mm)))+10
                nucMap_temp=imfill(nucMap_temp,'holes');
                nucMap_temp=bwareafilt(nucMap_temp,[500 100000]);
                nucleusMap=nucleusMap.*~nucMap_temp;
            end
        end
        nucleusMap=1-nucleusMap;
        % connect neighboring fractions
        se=strel('disk',2);
        nucMap_temp=imfill(nucleusMap,'holes');
        nucMap_temp=imdilate(nucMap_temp,se);
        nucMap_temp=imfill(nucMap_temp,'holes');
        se=strel('disk',3);
        nucMap_temp=imerode(nucMap_temp,se);
        se=strel('disk',1);
        nucMap_temp=imdilate(nucMap_temp,se);
        Reconimg3=Reconimg_selected.*repmat(nucMap_temp,[1 1 size(Reconimg_selected,3)]);
        Reconimg3=Reconimg3.*imerode(cellMap_selected,strel('sphere',3));
        %
        nucleusMap=((Reconimg3)>0);
        nucleusMap=imfill(nucleusMap,'holes');
        nucleusMap=bwareaopen(nucleusMap,500);
        %% Nucleoli Segmentation
        % segment nucleoli inside nucleus in RI tomograms, and segment nucleoplasm
        level = multithresh(Reconimg3(Reconimg3>0),2);
        nucleoliMap=((Reconimg3)>level(2));
        nucleoliMap=nucleoliMap.*cellMap_selected;
        % connect neighboring fractions
        se=strel('sphere',1);
        nucleoliMap=imfill(nucleoliMap,'holes');
        nucleoliMap=imdilate(nucleoliMap,se);
        nucleoliMap=imfill(nucleoliMap,'holes');
        nucleoliMap=imerode(nucleoliMap,se);
        nucleoliMap=bwareaopen(nucleoliMap,500).*imerode(cellMap_selected,strel('sphere',4));
        %
        [tempXR tempYR tempZR]=ind2sub(size(cellMap),find(nucleoliMap));
        nucleoliRI=(Reconimg_selected).*nucleoliMap;
        nucleoplasmRI=(Reconimg_selected).*nucleusMap.*~imdilate(nucleoliMap,strel('sphere',1));
        %% Peripheral Region Mask
        se=strel('sphere',8);
        periMap=imdilate((nucleusMap),se);
        se=strel('sphere',2);
        periMap=periMap-((nucleusMap+imdilate((nucleoliMap),se))>0);
        periMap=(periMap.*cellMap_selected)>0;
        %%
        nucleusRI=(Reconimg).*((nucleusMap+((nucleoliMap)))>0);
        periRI=(Reconimg).*periMap;
        periRI(:,:,max(tempZ)-3:end)=0;
        cytoplasmRI=(Reconimg_selected).*(~nucleusMap.*~imdilate((nucleoliMap),se));
        cytoplasmRI(:,:,max(tempZ)-3:end)=0; % discard RI under coverslip
        %% Nuclear Surface Area Calculation
        % Summation of areas of triangles in the mesh of isosurface of  nuclei
        p2 = isosurface(((nucleusMap+((nucleoliMap)))>0),0.8);
        v = p2.vertices;
        v(:,1:2)=v(:,1:2)*res3;
        v(:,3)=v(:,3)*res4;
        f = p2.faces;
        a = v(f(:,2),:) - v(f(:,1),:);
        b = v(f(:,3),:) - v(f(:,1),:);
        c = cross(a,b,2);
        NucArea = 0.5*sum(sqrt(sum(c.^2,2)));
        qDataTemp(20)=NucArea;
        %% Quantification #1: Entire cell
        %1: 3D Volume
        %2: 3D Dry mass
        %3: 3D RI
        %4: Surface Area
        qDataTemp(1)=sum(sum(sum(cellMap_selected))).*(res3.^2*res4);%volume
        qDataTemp(2)=sum(sum(sum((Reconimg-1.337).*cellMap_selected))).*(res3.^2*res4)/RII;%mass
        qDataTemp(3)=mean(Reconimg_selected(cellMap_selected));
        qDataTemp(4)=area;
        %% Quantification #2: RI of compartments
        %5: 3D Cytoplasm RI
        %6: 3D Peripheral RI
        %7: 3D Nucleus RI
        %8: 3D Nucleolus RI
        %9: 3D Nucleoplasm RI
        qDataTemp(5)=sum(sum(sum(cytoplasmRI)))./sum(sum(sum(cytoplasmRI>0.1)));
        qDataTemp(6)=sum(sum(sum(periRI(:,:,min(tempZR):max(tempZR)))))./sum(sum(sum(periRI(:,:,min(tempZR):max(tempZR))>0.1)));
        qDataTemp(7)=sum(sum(sum(nucleusRI)))./sum(sum(sum(nucleusRI>0.1)));
        qDataTemp(8)=sum(sum(sum(nucleoliRI)))./sum(sum(sum(nucleoliRI>0.1)));
        qDataTemp(9)=sum(sum(sum(nucleoplasmRI)))./sum(sum(sum(nucleoplasmRI>0.1)));
        %% Quantification #3: Volume and mass of compartments
        %10: 3D Cytoplasm volume
        %11: 3D Cytoplasm drymass
        %12: 3D Nucleus volume
        %13: 3D Nucleus drymass
        %14: 3D Nucleolus volume
        %15: 3D Nucleolus drymass
        %16: 3D Nucleoplasm volume
        %17: 3D Nucleoplasm drymass
        qDataTemp(10)=sum(sum(sum((cytoplasmRI>0)))).*(res3.^2*res4);%volume
        qDataTemp(11)=sum(sum(sum((Reconimg-1.337).*(cytoplasmRI>0)))).*(res3.^2*res4)/RII;%mass
        qDataTemp(12)=sum(sum(sum((nucleusRI>0)))).*(res3.^2*res4);%volume
        qDataTemp(13)=sum(sum(sum((Reconimg-1.337).*(nucleusRI>0)))).*(res3.^2*res4)/RII;%mass
        qDataTemp(14)=sum(sum(sum((nucleoliRI>0)))).*(res3.^2*res4);%volume
        qDataTemp(15)=sum(sum(sum((Reconimg-1.337).*(nucleoliRI>0)))).*(res3.^2*res4)/RII;%mass
        qDataTemp(16)=sum(sum(sum((nucleoplasmRI>0)))).*(res3.^2*res4);%volume
        qDataTemp(17)=sum(sum(sum((Reconimg-1.337).*(nucleoplasmRI>0)))).*(res3.^2*res4)/RII;%mass
        %% Quantification #4: Cell Cycle determination
        %18: Intensity in Green channel
        %19: Intensity in Red channel
        if max(max(squeeze(nucleusMap(:,:,end/2+z))))>0
            qDataTemp(18)=sum(sum(squeeze(combinedAnalysis(:,:,2).*squeeze(nucleusMap(:,:,end/2+z)))))./sum(sum(squeeze(nucleusMap(:,:,end/2+z))));
            qDataTemp(19)=sum(sum(squeeze(combinedAnalysis(:,:,1).*squeeze(nucleusMap(:,:,end/2+z)))))./sum(sum(squeeze(nucleusMap(:,:,end/2+z))));
        else
            qDataTemp(18)=sum(sum(squeeze(combinedAnalysis(:,:,2).*(squeeze(cytoplasmRI(:,:,end/2+z))>0.1))))./sum(sum((squeeze(cytoplasmRI(:,:,end/2+z))>0.1)));
            qDataTemp(19)=sum(sum(squeeze(combinedAnalysis(:,:,1).*(squeeze(cytoplasmRI(:,:,end/2+z))>0.1))))./sum(sum((squeeze(cytoplasmRI(:,:,end/2+z))>0.1)));
        end
        qData=[qData;qDataTemp];
    catch
    end
end