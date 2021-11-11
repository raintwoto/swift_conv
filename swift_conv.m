dirs=dir('images_crop/*.jpg');
refMean=128;
for i=1:size(dirs,1)
    img=imread(['images_crop/',dirs(i).name]);
img= histeq(imresize(rgb2gray(img),[480 640]));

     [sequence] = swift(1,200,img);
% 
mkdir('images_swift',strrep(dirs(i).name,'.jpg',''))
 for j=1:50
     imwrite(uint8(squeeze(sequence(j,:,:))),['images_swift/',strrep(dirs(i).name,'.jpg',['/',num2str(j),'.jpg'])]);     

 end
end