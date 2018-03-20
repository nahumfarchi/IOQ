function [theta_in_radians] = v_v_angle(u, v)
%V_V_ANGLE Return the angle between two vectors.

theta_in_radians = atan2(norm(cross(u,v)),dot(u,v));

end

