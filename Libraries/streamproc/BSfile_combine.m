function [ output ] = BSfile_combine( BSin1, BSin2, BSout )
    
    unix(['touch ' BSout]);

    fobj = fopen(BSout,'w');
    
    f_in1 = fopen(BSin1,'r');
    
    f1_line = fgets(f_in1);
    
    while ischar( f1_line )
        

        fprintf(fobj, f1_line);
        
        f1_line = fgets(f_in1);
    end
    
    fclose(f_in1);
    
    %% read in file 2
    
    f_in2 = fopen(BSin2,'r');
    
    f2_line = fgets(f_in2);
    
    Iskip = 1;
    
    while ischar( f2_line )
       
        if Iskip<10
            f2_line = fgets(f_in2);
            Iskip = Iskip+1;
            continue
        end
        fprintf(fobj, f2_line);
        f2_line = fgets(f_in2);
    end

    fclose(f_in1);
    
    fclose(fobj);
end