classdef Nodetc
    %NODE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        p_fid
        p_eid
        sign
    end
    
    methods
        function self = Nodetc(p_fid, p_eid, sign)
            self.p_fid = p_fid;
            self.p_eid = p_eid;
            self.sign = sign;
        end
    end
    
end