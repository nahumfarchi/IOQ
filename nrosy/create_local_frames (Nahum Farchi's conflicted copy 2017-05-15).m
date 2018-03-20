function [local_frames] = create_local_frames(mesh)
    %CREATE_LOCAL_FRAMES Create local frames for each face in the given mesh.
    % local_frames - 2*nF x 3, first nF rows are for e1, last nF rows are for
    % e2.

    e1 = mesh.V(mesh.F(:,2),:) - mesh.V(mesh.F(:,1),:);
    e1 = normalize_rows(e1);
    e2 = normalize_rows(cross(mesh.FNormals, e1, 2));
    local_frames = [e1; e2];
end

