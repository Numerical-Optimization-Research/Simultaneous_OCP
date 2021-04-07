function status = batt_gen_out(batt_num,destination_folder,input_files,input_name)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    status = -1;
    param = 4;
    
    choice = get_run_choice();
    [file_id,file_name,file_path] = create_file_id(choice,destination_folder,input_name);
    
    switch choice
        case 1
            run_from_file(file_id,file_name,file_path);
            status = fclose(file_id);
        case 2
            run_new_sim(file_id,param);
            produce_plot(file_path,file_name);
            status = fclose(file_id);
        otherwise
            fprintf('You have chosen to quit\n');
            return;
    end
    
    
end


function [file_out_id,file_name,tempfiledir] = create_file_id(choice,destinationfolder,input)

    out_file_original = '\EE6413_FinalProject';
    end_file = '.txt';
    file_name = strcat(out_file_original,input,end_file);
    tempfiledir = strcat(destinationfolder,'\',file_name);
    if isfile(tempfiledir) && (choice == 1)
        file_out_id = fopen(tempfiledir,'r'); 
    else
        file_out_id = fopen(tempfiledir,'w'); 
    end

end


function choice = get_run_choice()
    
    fprintf('How would you like to run.\n');
    fprintf('1)\tUse a previous Battery Model?\n');
    fprintf('2)\tRun a new Model (takes some time)\n');
    fprintf('0)\tQuit\n');
    choice = input('Your choice [0,1,2]:   ');

    fprintf('**************************************\n');
end


function run_from_file(file_id,file_name,file_path)

    if isfile(file_path)
        line = fgetl(file_id);
        if line == ""
            fprintf('The file was empty running a new model\n');
            fprintf('**************************************\n');
            run_new_sim(file_id);
        end
        produce_plot(file_path,file_name);
    else
        fprintf('There is no model stored, creating a new one.\n');
        fprintf('**************************************\n');
        run_new_sim(file_id);
        produce_plot(file_path,file_name);
    end

end


function run_new_sim(file_id,param)

    out = setup(param);
    dim = size(out);
    for i = 1:dim(2)
        
        fprintf(file_id,'%4.4f,',out(1,i));
        fprintf(file_id,'%4.4f,',out(2,i));
        for j = 1:dim(1) - 3
                fprintf(file_id,'%4.4f,',out(j+2,i));
        end
        fprintf(file_id,'%4.4f\n',out(j+3,i));
        
    end

end


function produce_plot(file_path,file_name)

    out = readmatrix(file_path);
    t = out(:,1);
    f = out(:,2);
    x1 = out(:,3);
    x2 = out(:,4);
    x3 = out(:,5);
    max_cost = max(f);
    max_show = ceil(max_cost/2)*2;
    
    figure('name','Optimum Cost');
    subplot(2,1,1)
    plot(t,f,'LineWidth',1.5);
    grid on;
    ylim([0 max_show]);
    xlabel('Time index');
    ylabel('Function Cost');

    subplot(2,1,2)
    plot(t,x1,t,x2,t,x3,'LineWidth',1.5);
    grid on;
    ylim([0 1]);
    xlabel('Time index');
    ylabel('X values');
    legend('DOD','T','SOC');
    file_final = erase(file_path,'.txt');
    name = strcat(file_final,'_plot');
    saveas(gcf,name,'png')
    
end


function out = setup(param)
    socMax = 0.98;
    socMin = 0.2;
    N = 10;
    time_max = 61;
    x(1) = socMax - socMin;     % DoD
    x(2) = 1;     % Temp
    x(3) = socMin + 0.45;    % SoC
    x(3) = enforceBCSOC(x(3));
    x_temp = x;
    x_opt = zeros(3,time_max);
    x_opt(1,1) = x(1);
    x_opt(2,1) = x(2);
    x_opt(3,1) = x(3);
    cost_opt = zeros(1,time_max);
    
    for i = 1:time_max

        if i == 1
            x_in = x;
        else
            x_in = x_opt(:,i-1);
        end

        [cost_brute,range] = brute_force(time_max,N,x_in,param);
        out = find_path(cost_brute,N);

        temp_compare = cost_brute(N/2,N/2,N/2);
        temp1 = cost_brute(out(:,N/2 - 1));
        temp2 = cost_brute(out(:,N/2 + 1));

        if temp1(1) < temp_compare(1)
            x_temp(:) = range(out(1,N/2 - 1),:);
            x_opt(:,i) = x_temp(:);
            cost_opt(i) = cost_brute(out(1,N/2 - 1));
        elseif temp2(1) < temp_compare(1)
            x_temp(:) = range(out(:,N/2 + 1),:);
            x_opt(:,i) = x_temp(:);
            cost_opt(i) = cost_brute(out(1,N/2 - 1));
        else
            x_opt(:,i) = x_in(:);
            cost_opt(i) = cost_opt(i-1);
        end

    end

    out = [0:time_max-1;cost_opt;x_opt];

end


function [cost,out] = brute_force(time_max,N,x,param)

    delta = N/(param*time_max);    
    range = zeros(N,3);
    range(N/2,1) = x(1);
    range(N/2,2) = x(2);
    range(N/2,3) = x(3);
    
    for i = 1:N
        for j = 1:length(x)
            if i < (N/2)
                range((N/2) - i,j) = range((N/2) - i + 1,j) - delta;
                if j == 3
                    range((N/2) - i,j) = enforceBCSOC(range((N/2) - i,j));
                elseif j == 1
                    range((N/2) - i,j) = enforceBCDOD(range((N/2) - i,j));
                end
            
            elseif i > (N/2)
                range(i,j) = range(i - 1,j) + delta; 
                if j == 3
                    range(i,j) = enforceBCSOC(range(i,j));
                elseif j == 1
                    range(i,j) = enforceBCDOD(range(i,j));
                end
            end
        end
    end
    
    cost = policy(range,N);
    out = range;
end


function p = policy(range,N)

    kt = 1.49e-6;
    Ceol = -800000;
    cost = zeros(N,N,N);
    for i = 1:N
        
        for j = 1:N
            
            for k = 1:N
                x_temp = [range(i,1),range(j,2),range(k,3)];
                cost(i,j,k) = Ceol*( -kt*eval(x_temp(2),2)*eval(x_temp(3),3)*( exp( -evalFd(x_temp) ) ) ) ;
            end
            
        end
        
    end
   
    p = cost;

end


function xmin = enforceBCSOC(xmin)
    socMin = 0.2;
    socMax = 0.98;
    for i = 1:length(xmin)
        if xmin(i) < socMin 
            xmin(i) = socMin;
        elseif xmin(i) > socMax
            xmin(i) = socMax;
        end
    end

end


function xmin = enforceBCDOD(xmin)
    socMin = 0.2;
    socMax = 0.98;
    dodMax = 0.85;
    dodMin = socMin;
    for i = 1:length(xmin)
        if xmin(i) < dodMin 
            xmin(i) = dodMin;
        elseif xmin(i) > dodMax
            xmin(i) = dodMax;
        end
    end

end


function funcVal = eval(x,i)
    a = 5.7905;
    b = -6.8292;
    c = 3.3209;
    d = 5.3696*10^-1;
    e = 6.1638e-2;
    kt = 1.49e-6;
    kT = 0.0693;
    Tref = 0.25;
    kSoC = 1.04;
    socRef = 0.5;
    delta_t = 0.1;
    switch i
        case 1
            funcVal = a*x^4 + b*x^3 + c*x^2 + d*x + e;
        case 2
            funcVal = exp(kT*(x - Tref)*(Tref/x));
        case 3
            funcVal = exp(kSoC*(x - socRef));
        case 4
            funcVal = kt*delta_t;
    end
end


function fd = evalFd(x)
    a = 5.7905;
    b = -6.8292;
    c = 3.3209;
    d = 5.3696*10^-1;
    e = 6.1638e-2;
    kt = 1.49e-6;
    kT = 0.0693;
    Tref = 0.25;
    kSoC = 1.04;
    socRef = 0.5;
    delta_t = 0.1;

    funcVal(1) = a*x(1)^4 + b*x(1)^3 + c*x(1)^2 + d*x(1) + e;
    funcVal(2) = exp(kT*(x(2) - Tref)*(Tref/x(2)));
    funcVal(3) = exp(kSoC*(x(3) - socRef));
    funcVal(4) = kt*delta_t;
    fd = -funcVal(4)*funcVal(3)*funcVal(2) - funcVal(1)*funcVal(3)*funcVal(2);

end


function path = find_path(range,N)
    out = zeros(3,N);
    for i = 1:N
        if i < (N/2)
            temp = range([i N/2],[i N/2],[i N/2]);
            [v,loc] = min(temp(:));
            out(:,i) =  i  + 1 - ind2sub(size(out),loc);
        elseif i > (N/2)
            temp = range([N/2 i],[N/2 i],[N/2 i]);
            [v,loc] = min(temp(:));
            out(:,i) = i - ind2sub(size(out),loc);
        end
    end
    path = out;
end


