fid = fopen('BZpath.dat', 'w');
steps = 40;
counter = 0;

dq = 1.0/steps;
figure

for j = 1:steps
    counter = counter + 1;
    vec1 = j*dq;
    vec2 = 0.0;
    plot(vec1, vec2, 'b*');
    hold on;
    fprintf(fid, '%2.0f %20.16f %20.16f \n', counter-1, vec1, vec2);
end

for j = 1:steps
    counter = counter + 1;
    vec1 = 1.0;
    vec2 = j*dq;
    plot(vec1, vec2, 'b*');
    hold on;
    fprintf(fid, '%4.0f %20.16f %20.16f \n', counter-1, vec1, vec2);
end

dq = -1.0/steps;
for j = 1:steps
    counter = counter + 1;
    vec1 = 1.0 + j*dq;
    vec2 = 1.0 + j*dq;
    
    if (vec1==0) && (vec2==0)
        vec1=1e-3;
        vec2=1e-3;
    end
    
    plot(vec1, vec2, 'b*');
    hold on;
    fprintf(fid, '%2.0f %20.16f %20.16f \n', counter-1, vec1, vec2);
end

xlabel('k1')
ylabel('k2')
hold off
fclose(fid);


