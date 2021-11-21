clear
c=dir('*.jpg');
for i=1:numel(c)
    mkdir(c(i).name(1:end-4));
    imFinal=imread(c(i).name);
%     imStart=randn(size(imFinal));
    imStart=imread('1.jpg');
    temp=imFinal;
    imFinal=imStart;
    imStart=temp;
    
    
    imSeq=minPhaseInterp(imStart,imFinal,[linspace(0.6,1,50)]);
%     figure;
    colormap gray;
    
    for iSeq=1:50,
        pic=((imSeq(:,:,iSeq)-min(min(imSeq(:,:,iSeq))))./(max(max(imSeq(:,:,iSeq))-min(min(imSeq(:,:,iSeq))))));
          PSF = fspecial('gaussian',60,10);
pic = edgetaper(pic,PSF);
        for jj=1:size(pic,2)
        pic(1:10,jj)=0.5:(pic(10,jj)-0.5)./9:pic(10,jj);
pic(341:350,jj)=pic(341,jj):(0.5-(pic(341,jj)))./9:0.5;
        end
        pic=permute(pic,[2 1]);
          for jj=1:size(pic,2)
        pic(1:10,jj)=0.5:(pic(10,jj)-0.5+0.0001)./9:(pic(10,jj))+0.0001;
pic(291:300,jj)=(pic(291,jj)-0.0001):(0.5-(pic(291,jj))+0.0001)./9:0.5;
          end
        pic=permute(pic,[2 1]);
        
        imwrite(pic,[c(i).name(1:end-4),'/',num2str(iSeq),'.jpg']);        
    end
end