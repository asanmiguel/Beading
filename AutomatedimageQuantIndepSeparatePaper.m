clc
clear all

%%
dirname='D:\PVD\Agingtry2V2\Agingtry2V2\Day 2\';
jjj=0;
d=dir(fullfile(dirname));
%%
for i=3:length(d)
    number=extractAfter(d(i).name,length(d(i).name)-1);
    number=str2double(number);
    if rem(number,2)==0
    
    dd1=strcat(dirname,d(i).name,'\','pred_masks\');
    dd=dir(fullfile(dd1,'*.png'));
    pp=strcat(dirname,d(i).name,'\','images\',d(i).name,'.tif');
    
    image2=imread(pp);
yy=0;
for l=1:size(dd)
dd(l).name;
p=extractBetween(ans,'score_','.png');
ss=str2num(p{1,1});
if ss>0.7
yy=yy+1;
zz(yy,1)=dd(l);
end
end
length(zz)
image=zeros(2048);
x=0;
xx=0;
for j=1:length(zz)
    if length(zz)
        jjj=jjj+1;
    end
    objectfile=strcat(dd1,zz(j).name);
    object=double(imread(objectfile));
    object=logical(object);
    C=regionprops(object,'all');
    area(j,i)=C(1).Area;
    area3(jjj,1)=C(1).Area;
        intensity(j,i)=mean(image2(C(1).PixelIdxList));
        cntr(j,1:2)=C(1).Centroid;
        dis=pdist(cntr);
        if area(j,i)>100
            xx=xx+1;
            up100(xx,i)=area(j,i);
        else
            x=x+1;
            less100(x,i)=area(j,i);
        end
    image=image+object;
end
data(i,1)=length(zz); %% Number of beads
image=logical(image);
% imshow(image)
   area_sort=area(1:j,i);
   area_sort=sort(area_sort);
  intensity_sort=intensity(1:j,i);
  intensity_sort=sort(intensity_sort);
    less100(less100==0)=NaN;
    up100(up100==0)=NaN;
     area(area==0)=NaN;
     intensity(intensity==0)=NaN;
    data(i,2)=sum(area(1:j,i)); % Sum of area occupied by puncta 
    data(i,3)=mean(area(1:j,i)); % Average of puncta size 
    data(i,4)=std(area(1:j,i)); % Standard deviation of bead size 
    data(i,5)=std(area(1:j,i))/mean(area(1:j,i)); % Standard deviation of mean of bead size 
    data(i,6)=std(area(1:j,i))/(length(zz)^0.5); % Standard Error of Mean (SEM)
    data(i,7)=prctile(area(1:j,i),90); % 90th percentile of bead size 
    data(i,8)=prctile(area(1:j,i),75); % 75th percentile of bead size 
    data(i,9)=prctile(area(1:j,i),50); % 50th percentile of bead size 
    data(i,10)=prctile(area(1:j,i),25); % 25th percentile of bead size 
    if xx~=0
    data(i,11)=mean(up100(1:xx,i)); % Average bead size for beads larger than 100 pixels 
    end
    if x~=0
    data(i,12)=mean(less100(1:x,i)); % Average bead size for beads smaller than 100 pixels 
    end
    data(i,13)=xx; % Number of beads with area larger than 100 pixels 
    data(i,14)=x;  % Number of beads with area smaller than 100 pixels
    data(i,15)=xx/(xx+x); % Percentage of beads with area larger than 100 Pixels 
    data(i,16)=x/(xx+x); % Percentage of beads with area smaller than 100 Pixels 
    data(i,17)=sum(dis)/length(dis); % Average distance betweeen beads
    data(i,18)=std(dis); % Standard Deviation of interbead distance 
    data(i,19)=std(dis)/mean(dis); % Standard deviation divided by mean of interbead distance 
    data(i,20)=std(dis)/(length(zz)^0.5); % Standard Error of Mean for distance 
    
    s=squareform(dis);
    y=zeros(length(s),8);
for iii=1:length(s)
for k=1:length(s)
if s(iii,k)<150
y(iii,1)=y(iii,1)+1;
elseif s(iii,k)<300&&s(iii,k)>150
y(iii,2)=y(iii,2)+1;
elseif s(iii,k)<450&&s(iii,k)>300
y(iii,3)=y(iii,3)+1;
else 
y(iii,4)=y(iii,4)+1;
end
end
y(iii,5)=y(iii,1)/(y(iii,1)+y(iii,2)+y(iii,3)+y(iii,4));
y(iii,6)=y(iii,2)/(y(iii,1)+y(iii,2)+y(iii,3)+y(iii,4));
y(iii,7)=y(iii,3)/(y(iii,1)+y(iii,2)+y(iii,3)+y(iii,4));
y(iii,8)=y(iii,4)/(y(iii,1)+y(iii,2)+y(iii,3)+y(iii,4));
end
data(i,21)=mean(y(:,5)); % Average number of beads with 150 pixels distance 
data(i,22)=mean(y(:,6)); % Average number of beads with 150 to 300 pixels distance 
data(i,23)=mean(y(:,7)); % Average number of beads with 300 to 450 pixels distance 
data(i,24)=mean(y(:,8)); % Average number of beads with more than 450 pixels distance 
data(i,25)=median(area(1:j,i)); % Median of area 
data(i,26)=max(area(1:j,i)); % Max size of the bead
data(i,27)=data(i,11)/data(i,12); % Mean area of bead larger than 100 pixels/Mean area of bead smaller than 100 pixels 
data(i,28)=mean(area_sort(1:floor(length(zz)/2),1)); % Mean area of smallest beads 
data(i,29)=mean(area_sort(floor(length(zz)/2)+1:length(zz),1)); % Mean area of largest beads 
data(i,30)=data(i,29)/data(i,28); % Mean area of largest beads/ Mean area of smallets beads 
data(i,31)=mean(intensity(1:j,i)); % Average of mean bead intensity 
data(i,32)=median(intensity(1:j,i)); % Median of mean bead intensity 
data(i,33)=max(intensity(1:j,i)); % Max of mean bead intensity 
data(i,34)=std(intensity(1:j,i)); % Standard deviation of mean bead internsity 
data(i,35)=std(intensity(1:j,i))/mean(intensity(1:j,i)); % Standard deviation of mean bead intensity/Average of mean bead intensity
data(i,36)=std(intensity(1:j,i))/(length(zz)^0.5); % Standard error of mean for mean bead intensity 
data(i,37)=prctile(intensity(1:j,i),90); % 90th percentile of mean bead intensity 
data(i,38)=prctile(intensity(1:j,i),75); % 75th percentile of mean bead intensity 
data(i,39)=prctile(intensity(1:j,i),50); % 50th percentile of mean bead intensity
data(i,40)=prctile(intensity(1:j,i),25); % 25th percentile of mean bead intensity 
data(i,41)=mean(intensity_sort(1:floor(length(zz)/2),1)); % Mean intensity of smallest beads 
data(i,42)=mean(intensity_sort(floor(length(zz)/2)+1:length(zz),1)); % Mean intensity of largest beads
data(i,43)=prctile(dis,90);
data(i,44)=prctile(dis,75);
data(i,45)=prctile(dis,50);
data(i,46)=prctile(dis,25);


image3=imadjust(image2);
BW=imbinarize(image3,0.12);
BW2=bwareaopen(BW,500);
BW3=imfill(BW2,'holes');
cc=regionprops(BW3,'all');
for iii=1:length(ans)
    
    area2(iii,i)=cc(iii).Area;
    
end
data(i,47)=max(area2(:,i));
clear dd zz dd1 pp dis cntr 
    else 
        
        dd1=strcat(dirname,d(i).name,'\','pred_masks\');
    dd=dir(fullfile(dd1,'*.png'));
    pp=strcat(dirname,d(i).name,'\','images\',d(i).name,'.tif');
    
    image2=imread(pp);
yy=0;
for l=1:size(dd)
dd(l).name;
p=extractBetween(ans,'score_','.png');
ss=str2num(p{1,1});
if ss>0.7
yy=yy+1;
zz(yy,1)=dd(l);
end
end
length(zz)
image=zeros(2048);
x=0;
xx=0;
for j=1:length(zz)
    
    if length(zz)
        jjj=jjj+1;
    end
    objectfile=strcat(dd1,zz(j).name);
    object=double(imread(objectfile));
    object=logical(object);
    C=regionprops(object,'all');
    area(j,i)=C(1).Area;
    area3(jjj,1)=C(1).Area;
        intensity(j,i)=mean(image2(C(1).PixelIdxList));
        cntr(j,1:2)=C(1).Centroid;
        dis=pdist(cntr);
        if area(j,i)>100
            xx=xx+1;
            up100(xx,i)=area(j,i);
        else
            x=x+1;
            less100(x,i)=area(j,i);
        end
    image=image+object;
end
data2(i,1)=length(zz);
image=logical(image);
   area_sort=area(1:j,i);
   area_sort=sort(area_sort);
%    length(area_sort)
  intensity_sort=intensity(1:j,i);
  intensity_sort=sort(intensity_sort);
    less100(less100==0)=NaN;
    up100(up100==0)=NaN;
     area(area==0)=NaN;
     intensity(intensity==0)=NaN;
    data2(i,2)=sum(area(1:j,i)); % Sum of area occupied by puncta 
    data2(i,3)=mean(area(1:j,i)); % Average of puncta size 
    data2(i,4)=std(area(1:j,i)); % Standard deviation of bead size 
    data2(i,5)=std(area(1:j,i))/mean(area(1:j,i)); % Standard deviation of mean of bead size 
    data2(i,6)=std(area(1:j,i))/(length(zz)^0.5); % Standard Error of Mean (SEM)
    data2(i,7)=prctile(area(1:j,i),90); % 90th percentile of bead size 
    data2(i,8)=prctile(area(1:j,i),75); % 75th percentile of bead size 
    data2(i,9)=prctile(area(1:j,i),50); % 50th percentile of bead size 
    data2(i,10)=prctile(area(1:j,i),25); % 25th percentile of bead size 
    if xx~=0
    data2(i,11)=mean(up100(1:xx,i)); % Average bead size for beads larger than 100 pixels 
    end
    if x~=0
    data2(i,12)=mean(less100(1:x,i)); % Average bead size for beads smaller than 100 pixels 
    end
    data2(i,13)=xx; % Number of beads with area larger than 100 pixels 
    data2(i,14)=x;  % Number of beads with area smaller than 100 pixels
    data2(i,15)=xx/(xx+x); % Percentage of beads with area larger than 100 Pixels 
    data2(i,16)=x/(xx+x); % Percentage of beads with area smaller than 100 Pixels 
    data2(i,17)=sum(dis)/length(dis); % Average distance betweeen beads
    data2(i,18)=std(dis); % Standard Deviation of interbead distance 
    data2(i,19)=std(dis)/mean(dis); % Standard deviation divided by mean of interbead distance 
    data2(i,20)=std(dis)/(length(zz)^0.5); % Standard Error of Mean for distance 
    
    s=squareform(dis);
    y=zeros(length(s),8);
for iii=1:length(s)
for k=1:length(s)
if s(iii,k)<150
y(iii,1)=y(iii,1)+1;
elseif s(iii,k)<300&&s(iii,k)>150
y(iii,2)=y(iii,2)+1;
elseif s(iii,k)<450&&s(iii,k)>300
y(iii,3)=y(iii,3)+1;
else 
y(iii,4)=y(iii,4)+1;
end
end
y(iii,5)=y(iii,1)/(y(iii,1)+y(iii,2)+y(iii,3)+y(iii,4));
y(iii,6)=y(iii,2)/(y(iii,1)+y(iii,2)+y(iii,3)+y(iii,4));
y(iii,7)=y(iii,3)/(y(iii,1)+y(iii,2)+y(iii,3)+y(iii,4));
y(iii,8)=y(iii,4)/(y(iii,1)+y(iii,2)+y(iii,3)+y(iii,4));
end
% length(y)
data2(i,21)=mean(y(:,5)); % Average number of beads with 150 pixels distance 
data2(i,22)=mean(y(:,6)); % Average number of beads with 150 to 300 pixels distance 
data2(i,23)=mean(y(:,7)); % Average number of beads with 300 to 450 pixels distance 
data2(i,24)=mean(y(:,8)); % Average number of beads with more than 450 pixels distance 
data2(i,25)=median(area(1:j,i)); % Median of area 
data2(i,26)=max(area(1:j,i)); % Max size of the bead
data2(i,27)=data2(i,11)/data2(i,12); % Mean area of bead larger than 100 pixels/Mean area of bead smaller than 100 pixels 
data2(i,28)=mean(area_sort(1:floor(length(zz)/2),1)); % Mean area of smallest beads 
data2(i,29)=mean(area_sort(floor(length(zz)/2)+1:length(zz),1)); % Mean area of largest beads 
data2(i,30)=data2(i,29)/data2(i,28); % Mean area of largest beads/ Mean area of smallets beads 
data2(i,31)=mean(intensity(1:j,i)); % Average of mean bead intensity 
data2(i,32)=median(intensity(1:j,i)); % Median of mean bead intensity 
data2(i,33)=max(intensity(1:j,i)); % Max of mean bead intensity 
data2(i,34)=std(intensity(1:j,i)); % Standard deviation of mean bead internsity 
data2(i,35)=std(intensity(1:j,i))/mean(intensity(1:j,i)); % Standard deviation of mean bead intensity/Average of mean bead intensity
data2(i,36)=std(intensity(1:j,i))/(length(zz)^0.5); % Standard error of mean for mean bead intensity 
data2(i,37)=prctile(intensity(1:j,i),90); % 90th percentile of mean bead intensity 
data2(i,38)=prctile(intensity(1:j,i),75); % 75th percentile of mean bead intensity 
data2(i,39)=prctile(intensity(1:j,i),50); % 50th percentile of mean bead intensity
data2(i,40)=prctile(intensity(1:j,i),25); % 25th percentile of mean bead intensity 
data2(i,41)=mean(intensity_sort(1:floor(length(zz)/2),1)); % Mean intensity of smallest beads 
data2(i,42)=mean(intensity_sort(floor(length(zz)/2)+1:length(zz),1)); % Mean intensity of largest beads
data2(i,43)=prctile(dis,90);
data2(i,44)=prctile(dis,75);
data2(i,45)=prctile(dis,50);
data2(i,46)=prctile(dis,25);
% sahand=dis;
% sahand2=cntr;
image3=imadjust(image2);
BW=imbinarize(image3,0.12);
BW2=bwareaopen(BW,500);
BW3=imfill(BW2,'holes');
cc=regionprops(BW3,'all');
for iii=1:length(ans)
    
    area2(iii,i)=cc(iii).Area;
    
end
data2(i,47)=max(area2(:,i));
clear dd zz dd1 pp dis cntr
    end
end 

