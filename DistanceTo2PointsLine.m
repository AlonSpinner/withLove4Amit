function d = DistanceTo2PointsLine(p1,p2,q)
t = (p1-p2)/vecnorm(p1-p2);
v_t = ((q-p1)*t')*t;
v_n = (q-p1) - v_t;
d = vecnorm(v_n,2,2);
end