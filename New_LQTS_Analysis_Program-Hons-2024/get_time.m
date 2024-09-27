function [ time ] = get_time( val,sample_rate )
%GET_TIME Summary of this function goes here
%   Detailed explanation goes here
  time = val * 1000/sample_rate;
  
end

