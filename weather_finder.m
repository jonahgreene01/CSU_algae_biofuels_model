function [file_ids,weather_script,coordinates_out] = weather_finder(County,State)
    load background_data.mat county_dat
    idx_state = find(strcmp(county_dat(:,4),State));
    idx_county = find(strcmp(county_dat(:,3), County));
    idx_files = intersect(idx_county,idx_state); 
    
    file_ids = county_dat(idx_files,1); 
    weather_script = zeros(length(file_ids)*8760,6);
    coordinates_out = zeros(length(file_ids), 2); 

    for i = 1:length(file_ids)
        file = "Merger/" + file_ids(i) + ".csv";
        weather = readtable(file,'Range','G3:K8762','ReadVariableNames',false); 
        weather_script((8760*(i-1)+1):(i*8760),:)  = [repmat(i,8760,1), weather{:,:}]; 
        coord =  readtable(file,'Range','F1:G1','ReadVariableNames',false);
        coordinates_out(i,:) = coord{:,:};
    end
end

