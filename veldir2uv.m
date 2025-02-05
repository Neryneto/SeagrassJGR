function [u,v] = veldir2uv(vel,dir)
u = vel .* cosd(dir);
v = vel .* sind(dir);
end