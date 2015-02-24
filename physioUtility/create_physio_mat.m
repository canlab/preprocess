files = filenames('*.txt');
physio = read_physio_data(files)

% convert pulse to hr


save('physio', 'physio');