function current = Control(Bfieldnav,pqrdotnav)
    Magtorquer_params
    k = 67200;          %%Calculated by assuming current = 40milliamps
    
    current = k*cross(pqrdotnav,Bfieldnav)/(n*A);
     
    %%If magnetic saturation takes place at 40 milli amps
    if sum(abs(current)) > 0.04             
        current = (current/norm(current))*0.04;
    end
end
