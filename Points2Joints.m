classdef Points2Joints < handle
    properties
        LinkLength = 35;
        maxDistance = 3;
        dTheta = 0.01; %rad
        p0 = [282.9951  228.9515];
        points = [];
        pJoints = [];
        bezOrder = 4;
        weightRange = [1,2]; %[weight on close points, weight on far points]
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

                MaxwInliersAmount = 0;
                sumD = inf;
                for theta = thetaVec
                    pEndTest = pStart + obj.LinkLength*[cos(theta),sin(theta)];
                    [D2line,d_t] = obj.DistanceTo2PointsLine(pStart,pEndTest,points4Link);
                    inlierIndcies = D2line < obj.maxDistance & d_t > 0;
                    D2Start = vecnorm(points4Link-pStart,2,2);
                    w = obj.mapfun(D2Start,min(D2Start),max(D2Start),obj.weightRange(1),obj.weightRange(2));
                    wInliersAmount=sum(w.*inlierIndcies);
                    iterartionSumD = sum(D2line(inlierIndcies));
                    if wInliersAmount > MaxwInliersAmount ||...
                        wInliersAmount == MaxwInliersAmount && iterartionSumD < sumD
                        sumD = iterartionSumD;
                        MaxwInliersAmount = wInliersAmount;
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

            CP = obj.BezFit(obj.pJoints,obj.bezOrder);
            u=linspace(0,1,1000);
            bezCurve = obj.EvalBezCrv_DeCasteljau(u,CP);
            plot(Ax,bezCurve(:,1),bezCurve(:,2),...
                'linestyle','-','linewidth',2,'color',[0.6,0.3,0.6]);
        end
    end
    methods (Access = private, Static) %building joints functions
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
    methods (Access = private, Static) %bezier functions
        function R=EvalBezCrv_DeCasteljau(q,Q)
            %{
            Evaluates Bezier Curve by given nodes and parameter
            value
            
            Q - nodes in format [x,y] dimension matrix. top row is q=0.
            q - running parameter of bezier curve. 0=<q<=1
            R - [x,y] or [x,y,z] format. point on bezier curve
            https://pages.mtu.edu/~shene/COURSES/cs3621/NOTES/spline/Bezier/de-casteljau.html
            %}

            %currently not used in main
            q = q(:);
            R = zeros(length(q),size(Q,2));
            n=size(Q,1)-1; %degree of bezier polynomial
            for ii = 1:length(q)
                qii = q(ii);
                Qii = Q;
                for kk=1:(n+1)
                    for jj=1:(n+1-kk) %doesnt enter in first iteration
                        Qii(jj,:)=(1-qii)*Qii(jj,:)+qii*Qii(jj+1,:);
                    end
                end
                R(ii,:)=Qii(1,:);
            end
        end
        function CP = BezFit(xy,order,xy0,xyn)
            %xy - [x,y] size mx2
            %order - scalar, order of bezier curve

            N = size(xy,1);
            U=deal(zeros(N,order+1));
            q=linspace(0,1,N); %assume points should be evenly spaced on curve
            for i=1:N
                for j=1:(order+1)
                    U(i,j)=q(i)^(j-1); %[q^0, q^1, q^2, q^3]
                end
            end
            B = Points2Joints.calcBezMatrix(order);
            pinv=@(A) (A'*A)\A'; %psuedo inverse matrix (moore-pensrose)

            if nargin>2 %enforce starting and ending points for nonordered optimization
                xy(1,:) = xy0;
                xy(end,:) = xyn;
            end

            %XY=U*B*CP
            CP=zeros(order+1,size(xy,2)); %initalize
            for dim=1:size(xy,2) %R2 dimensions (x,y), they are indepednent of each other
                CP(:,dim)=B\pinv(U)*xy(:,dim);
            end
        end
        function B=calcBezMatrix(order,symbolic)
            %{
            Note: heavily relies on symbolic toolbox if not classic cases
            
            
            Input:
            order - bezier order
            
            note: bernsteinMatrix is in fact the following:
            n=BezO;
            for k=0:n
               b(k+1)= nchoosek(n,k)*t^k*(1-t)^(n-k);
            end
            %}

            if nargin<2
                symbolic = false;
            end

            if symbolic
                B = calcBsymbolic(order);
                return
            end

            %numeric for simple classic cases
            switch order
                case 1
                    B = [-1 0;
                        -1 1];
                case 2
                    B = [1 0 0;
                        -2 2 0;
                        1 -2 1];
                case 3
                    B = [1 0 0 0;
                        -3 3 0 0;
                        3 -6 3 0;
                        -1 3 -3 1];
                case 4
                    B = [1 0 0 0 0;
                        -4 4 0 0 0;
                        6 -12 6 0 0
                        -4 12 -12 4 0;
                        1 -4 6 -4 1];
                otherwise
                    B = Points2Joints.calcBsymbolic(order);
            end
        end

        function B = calcBsymbolic(order)
            syms t
            B=zeros(order+1,order+1);
            b=fliplr(bernsteinMatrix(order,t));
            for j=1:length(b)
                B(j,:)=coeffs(b(j),t,'all');
            end
        end

        function output = mapfun(value,fromLow,fromHigh,toLow,toHigh)
            % from https://viewer.mathworks.com/?viewer=plain_code&url=https%3A%2F%2Fwww.mathworks.com%2Fmatlabcentral%2Fmlc-downloads%2Fdownloads%2Fsubmissions%2F61213%2Fversions%2F1%2Fcontents%2Fmapfun.m&embed=web
            % Created by David Coventry on 1/19/2017
            output = (value - fromLow) .* (toHigh - toLow) ./ (fromHigh - fromLow) + toLow;
        end
    end
end