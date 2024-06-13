function [val,isTerminal,direction] = flipEvent(t,u)
% have we crossed zero?
val = u(2);
% terminate when current height - ground height == 0
isTerminal = 0;
% approach from above, below, or either?
direction = 1;