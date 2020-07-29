% What this folder is: a collection of all .m files that we will use to
% extract spike trains of 1s (spike) and 0s (no spike), and combine those
% spike trains across sessions and animals. And finally, after we have a
% cell of spike trains across all animals, we can convert those into spike
% rates per bin of time.
%
% A +folder creates a matlab package. See this https://www.mathworks.com/help/matlab/matlab_oop/scoping-classes-with-packages.html
%
% To run codes in this package, you would do:
% spikes.getSpikeTrain()
% spikes.combineSesssions()
% spikes.combineAnimals()
%
% In a sense, this is a way of organizing functions who are being used for
% some common purpose.
