classdef battery
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        power_max
        power_current
        direction = 'no';
        currentSOC
        maxSOC 
        minSOC
        inv_effic = 0.950813662333302;
        voc_chrg = 420;  % Volts
        voc_dischrg = 350;  % Volts
        r_chrg = (420/11.9)*10^-3;
        r_dischrg = (350/14.3)*10^-3;
        deltaT
        delta_e
    end
    
    methods
        function obj = battery(soc,deltaT)
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here
            obj.maxSOC = 0.95;
            obj.minSOC = 0.2;
            obj.power_max = 9.3;
            obj.currentSOC = soc;
            obj.checkSOC();
            obj.power_current = 1;
            obj.power_current = obj.currentSOC*obj.power_max;
            obj.deltaT = deltaT;
            obj.delta_e = 0.0;
        end
        
        % Charge power should be negative, and discharge power would be
        % positive.
        function obj = updatePower(obj,power)
            temp = obj.power_current;
            gamma_effic = obj.batt_effic();
            
            delta_e_temp = -gamma_effic^(-1)*obj.inv_effic^(-1)*power*obj.deltaT;
            obj.power_current = temp + delta_e_temp;
            obj.currentSOC = obj.power_current/obj.power_max;
            obj.checkSOC();
            
            if (obj.delta_e) > 0
                obj.direction = 'charge';
            else
                obj.direction = 'discharge';
            end
            
            if obj.currentSOC <= 0.2
                obj.delta_e = 0;
            else
                obj.delta_e = delta_e_temp;
            end
        end
        
        
        function delta_e = get_change(obj)
            delta_e = obj.delta_e;
        end
        
        
        function obj = checkSOC(obj)
            
            if obj.currentSOC > obj.maxSOC
                obj.currentSOC = obj.maxSOC;
            elseif obj.currentSOC < obj.minSOC
                obj.currentSOC = obj.minSOC;
            else
                obj.currentSOC = obj.currentSOC; 
            end
            
        end
        
        
        % calculates gamma_bat seen in equation 8 of fortenbacher
        function efficiency = batt_effic(obj)
            
            if strcmp(obj.direction,'charge')
                
                power_temp = obj.power_current*(10^3);
                temp = (obj.voc_chrg)^2 - 4*obj.r_chrg*power_temp;
                efficiency = 1 - abs((obj.voc_chrg - sqrt(temp))/(2*obj.voc_chrg));
                
            elseif strcmp(obj.direction,'discharge')
                
                power_temp = obj.power_current*(10^3);
                temp = (obj.voc_dischrg)^2 - 4*obj.r_dischrg*power_temp;
                efficiency = 1 - abs((obj.voc_dischrg - sqrt(temp))/(2*obj.voc_dischrg));
                
            else
                
                power_temp = obj.power_current*(10^3);
                temp = (obj.voc_chrg)^2 - 4*obj.r_chrg*power_temp;
                efficiency = 1 - abs((obj.voc_chrg - sqrt(temp))/(2*obj.voc_chrg));
                
            end
            
        end
        
        
        function obj = setSOC(obj,soc)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            obj.currentSOC = soc;
        end
        
        
        function power = returnMaxPower(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            power = obj.power_max;
        end
        
        
        function soc = returnCurrentSOC(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            soc = obj.currentSOC;
        end
        
        
    end
end

