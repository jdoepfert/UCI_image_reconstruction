
function handle=ShowImages(imagestack, Nrows, NCols,DisplayRange)

if ~exist('Nrows')|| isempty(Nrows)
    Nrows=NaN;
end

if ~exist('NCols')|| isempty(NCols)
    NCols=NaN;
end

if ~exist('DisplayRange') || isempty(DisplayRange)
    DisplayRange=[min(min(min(imagestack))) max(max(max(imagestack)))+0.00000000001];
end

[sy,sx]=size(imagestack(:,:,1)); 
imagestack=reshape(imagestack,sy,sx,1,size(imagestack,3));
handle=montage(imagestack,'Size', [Nrows NCols],'DisplayRange',DisplayRange);

end
