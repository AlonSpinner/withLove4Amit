classdef Points2Joints < handle
    properties
        LinkLength = 35;
        dTheta = 0.01; %rad
        p0 = [282.9951  228.9515];
        maxDistance = 3;
        points = [];
        pJoints = [];
    end
    methods (Access = public)
        function obj = Points2Joints(points)
            obj.points = points;
        end
        function findJoints(obj)
            mpoints = obj.points;
            mpJoints = obj.p0;
            thetaVec = 0:obj.dTheta:2*pi;
            while size(mpoints,1) > 0
                pStart = mpJoints(end,:);
                [points4Link,indcies] = obj.FilterPointsByRadius(pStart,obj.LinkLength,mpoints);
                if isempty(points4Link), disp('no points found around last joint'); break; end
                
                maxInliersAmnt = 0;
                for theta = thetaVec
                    pEndTest = pStart + obj.LinkLength*[cos(theta),sin(theta)];
                    [D,d_t] = obj.DistanceTo2PointsLine(pStart,pEndTest,points4Link);
                    iterationInliersAmnt = sum(D < obj.maxDistance & d_t > 0);
                    if iterationInliersAmnt > maxInliersAmnt 
                        maxInliersAmnt = iterationInliersAmnt;
                        pEnd = pEndTest;
                    end
                end
                mpJoints(end+1,:) = pEnd;
                mpoints(indcies,:) = []; %destroy points that we already went through
                if isempty(mpoints), break, end
            end
        
              %delete last point if there is nothing really near it
%             lastJoint = mpJoints(end,:);
%             pNext2Last = obj.FilterPointsByRadius(lastJoint,obj.LinkLength/2,obj.points);
%             if size(pNext2Last,1) < obj.maxDistance*2
%                 mpJoints(end,:)= [];
%             end

            obj.pJoints = mpJoints;
        end

        function plot(obj)
            Fig=figure('color',[0,0,0]);
            Ax=axes(Fig,'color',[0,0,0],'XColor',[1,1,1],'YColor',[1,1,1]);
            xlabel(Ax,'x'); ylabel(Ax,'y'); zlabel(Ax,'z');
            axis(Ax,'equal'); grid(Ax,'on'); hold(Ax,'on');

            scatter(Ax,obj.points(:,1),obj.points(:,2),10,'filled');
            scatter(Ax,obj.p0(:,1),obj.p0(:,2),500,'yellow');
            plot(Ax,obj.pJoints(:,1),obj.pJoints(:,2),'linewidth',4,...
                'Marker','o','MarkerFaceColor',[1,1,1],'MarkerSize',10);
            
        end
    end
    methods (Access = private, Static)
        function [Pnew,Idx] = FilterPointsByRadius(p0,r,P)
            Idx = vecnorm(P-p0,2,2) < r;
            Pnew = P(Idx,:);
        end
        function [d,d_t] = DistanceTo2PointsLine(p1,p2,q)
            t = (p2-p1)/vecnorm(p2-p1);
            d_t = (q-p1)*t';
            v_t = d_t*t;
            v_n = (q-p1) - v_t;
            d = vecnorm(v_n,2,2);
        end

    end
end