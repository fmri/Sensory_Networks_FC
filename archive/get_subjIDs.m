% Steal subjIDs from Vaibhavt's directory 

base = '/projectnb/somerslab/vaibhavt/Projects/Trifloc/Data/SubjectsBehavioralFiles/';

filenames = dir(base);
filenames = {filenames.name};
subjIDs = filenames(3:end);

save('subjIDs.mat', 'subjIDs');

% Create directory for each subjID

for ii = 1:length(subjIDs)

    mkdir(['./data/', subjIDs{ii}])
    mkdir(['./data/', subjIDs{ii}, '/anat'])
    mkdir(['./data/', subjIDs{ii}, '/func'])

end