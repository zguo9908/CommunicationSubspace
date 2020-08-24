function [running, sleeping] = getRunningSessions(animal)

load(animal + "task01.mat");
running = [];
sleeping = [];
task = task{1};

for i = 1:numel(task)
    if task{i}.type == "sleep"
        sleeping = [sleeping, i];
    else
        running = [running, i];
    end
end

end

