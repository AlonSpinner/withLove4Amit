%% Points1
load(fullfile('data','Points1.mat'));
obj = Points2Joints(points);
obj.LinkLength = 15;
obj.maxDistance = 4;
obj.findJoints();
obj.plot();
disp(obj.pJoints);
%% Points2
load(fullfile('data','Points2.mat'));
obj = Points2Joints(points);
obj.bezOrder = 4;
obj.weightRange = [1,1.5];
obj.LinkLength = 30;
obj.findJoints();
obj.plot();
disp(obj.pJoints);

%% Points3
load(fullfile('data','Points3.mat'));
obj = Points2Joints(points);
obj.LinkLength = 20;
obj.findJoints;
obj.plot();
disp(obj.pJoints);