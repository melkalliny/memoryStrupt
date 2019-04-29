function [start_idx end_idx count] = consecutive_ones(vector)
%
%[Function Description]
%Finds the number of consecutive ones in a binary signal. Returns the
%starting and ending points of the consecutive set and also the count. Since it
%does not use any for loop. It is pretty fast
%
%[Input]
%vector = A binary signal
%[Output]
%
%Author - Shreyes

temp = diff([0 vector' 0]);
start_idx = find(temp == 1);
end_idx = find(temp == -1);
count = end_idx - start_idx;
end_idx = end_idx -1;