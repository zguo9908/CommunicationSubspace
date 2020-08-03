animalname = "KL8"; %% change here

% animaldef
animalinfo = animaldef(animalname);

% getSpectralBehaviorData
animal = {animalname};
data = raw.getSpectralBehaviorData(animal);

% getGlobalRipple
globalripple = generateGlobalRipple(animalname);
animalpath = ("C:\Users\BrainMaker\commsubspace\SingleDayExpt\"+animalname+"_direct");
cd (animalpath)
save(animal+"globalripple01", "globalripple");

% get avgeeg