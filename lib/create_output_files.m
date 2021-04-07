function out_file = create_output_files(destinationfolder,inputsize)
    
    out_file_original = '\EE6413_FinalProject';
    end_file = '.txt';
    fileOutID = zeros(inputsize,1);
    for i = 1:inputsize
        tempfiledir = strcat(destinationfolder,'\',out_file_original,'_file',string(i),end_file);
        fileOutID(i) = fopen(tempfiledir,'w');
    end
    out_file = fileOutID(:);

end



