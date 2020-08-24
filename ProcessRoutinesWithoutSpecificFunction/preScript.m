animalname = "ER1"; %% change here
animalpath = ("C:\Users\BrainMaker\commsubspace\SingleDayExpt\"+animalname+"_direct");
cd (animalpath)

% animaldef
animalinfo = animaldef(animalname);

% getSpectralBehaviorData
animal = {animalname};
data = raw.getSpectralBehaviorData(animal, "epochtype", "run");
save(animalname+'spectralBehavior.mat','-struct','data', '-v7.3');
% save(animal+"spectralBehavior","data");

% getGlobalRipple
globalripple = generateGlobalRipple(animalname);

save(animalname+"globalripple01", "globalripple");

% get avgeeg