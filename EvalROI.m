function [Mean,Std] = EvalROI(Image,RoiMap)

%   evaluate mean and stddeviation in a roi given by matrix RoiMap

    Mean=mean(Image(RoiMap>0));
    Std=std(double(Image(RoiMap>0)));

end


