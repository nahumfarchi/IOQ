function frame_diffs = compute_frame_diffs(m, debug)
    %COMPUTE_FRAME_DIFFS Compute the angles between local frames of adjacent
    %faces in the given mesh.

    if nargin < 2
        debug = false;
    end

    EV = m.EVAdj;
    EF = m.EFAdj;
    face_normals = m.FNormals;
    V = m.V;
    F = m.F;
    frame_diffs = zeros(m.nE, 1);
    assert(m.nE == size(EF, 1))
    assert(m.nE == size(EV, 1))

    for eid = 1:m.nE
        if m.isBoundaryEdge(eid)
            continue;
        end

        fid0 = EF(eid, 1);
        fid1 = EF(eid, 2);

        N0 = face_normals(fid0, :);
        %N1 = face_normals(fid1, :);

        v1 = EV(eid, 1);
        v2 = EV(eid, 2);

        %common_edge = normalize_rows(V(EV(eid, 2),:) - V(EV(eid, 1),:));
        common_edge = normalize_rows(V(v2, :) - ...
                                     V(v1, :));

        % Transform the first triangle so that:
        %   x axis is the common edge
        %   z axis is N0
        % the triangle now sits in the XY plane.
        o = V(EV(eid, 1),:);
        tmp = -cross(N0, common_edge, 2);
        P = [common_edge; tmp; N0];

        V0 = [V(F(fid0,1),:) - o; ...
              V(F(fid0,2),:) - o; ...
              V(F(fid0,3),:) - o];

        % Apply the same transformation to the second triangle.
        V1 = [V(F(fid1,1),:) - o; ...
              V(F(fid1,2),:) - o; ...
              V(F(fid1,3),:) - o];

        V0 = (P * V0')';
        V1 = (P * V1')';

        % Make sure all z-coordinates are zero.
        assert(V0(1,3) < 1e-10)
        assert(V0(2,3) < 1e-10)
        assert(V0(3,3) < 1e-10)

        if debug
            figure;
            subplot(211);
            %self.debug_plot(V0, V1, P);
            title('Before rotation');
        end

        % Rotate the second triangle so that it sits in XY plane as well.
        % We have to find the vertex that is not on the common edge.
        other_vert_row = -1;
        for i = 1:3
            if F(fid1, i) ~= EV(eid, 1) && F(fid1, i) ~= EV(eid, 2)
                other_vert_row = i;
            end
        end
        alpha = -atan2(V1(other_vert_row, 3), V1(other_vert_row, 2));
        R = [1, 0, 0; ...
             0, cos(alpha), -sin(alpha); ...
             0, sin(alpha), cos(alpha)];
        V1 = (R*V1')';

        if debug
            subplot(212);
            %self.debug_plot(V0, V1, P);
            title('Rotated');
        end

        % Make sure V1 now sits in the XY plane
        assert(V1(1,3) < 1e-10);
        assert(V1(2,3) < 1e-10);
        assert(V1(3,3) < 1e-10);

        % Calculate the angle between the two local frames.
        % Reminder: for each frame, e0 is the first edge 
        % of the triangle
        f0_e0 = normalize_rows(V0(2, :) - V0(1, :));
        f1_e0 = normalize_rows(V1(2, :) - V1(1, :));
        angle_between_frames = atan2(f1_e0(2), f1_e0(1)) - atan2(f0_e0(2), f0_e0(1));

        % Make sure that rotation is correct.
        R = [cos(angle_between_frames), -sin(angle_between_frames); ...
             sin(angle_between_frames), cos(angle_between_frames)];
        f0_e0_rot = (R*f0_e0(1:2)')';
        assert(norm(f0_e0_rot - f1_e0(1:2)) < 1e-10)

        % It's all good.
        %self.k(eid) = angle_between_frames;
        frame_diffs(eid) = angle_between_frames;
    end
end