function [ attr ] = ncget_att_ic( fname, varName, group, att  )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    attname = matlab.lang.makeValidName(att,'Prefix','ZZZ');
    a = ncget_atts_ic(fname, varName, group);
    attr = a.(attname);

end

