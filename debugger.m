figure;
scatter(obj.points(:,1),obj.points(:,2),50)
hold on
scatter(pEnd(:,1),pEnd(:,2),50,'filled')
scatter(pStart(:,1),pStart(:,2),50,'filled')
scatter(points4Link(:,1),points4Link(:,2),10,'filled')
vecnorm(pEnd-pStart,2)