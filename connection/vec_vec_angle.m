function [angle] = vec_vec_angle(a, b)
%VEC_VEC_ANGLE Return the angle between two vectors.

[n, m] = size(a);

if n == 1 || m == 1
    angle = atan2(norm(cross(a,b)),dot(a,b));
else
    angle = atan2(row_norm(cross(a,b,2)),dot(a,b,2));
end

end

