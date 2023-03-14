% read data from .cfg file
fid=fopen('final.cfg');
line = split(fgetl(fid));
nparticle = str2double(line{5});

skipNRows(fid,2);
H_matrix = zeros(9,1);
for i=1:9
    line = split(fgetl(fid));
    H_matrix(i) = str2double(line{3});
end

skipNRows(fid,5);
pos = zeros(nparticle,3);
type = zeros(nparticle,1);

for i=1:nparticle
    line = split(fgetl(fid));
    s1 = str2double(line{1});
    s2 = str2double(line{2});
    s3 = str2double(line{3});
    pos(i,1) = s1*H_matrix(1)+s2*H_matrix(4)+s3*H_matrix(7);
    pos(i,2) = s1*H_matrix(2)+s2*H_matrix(5)+s3*H_matrix(8);
    pos(i,3) = s1*H_matrix(3)+s2*H_matrix(6)+s3*H_matrix(9);
    type(i) = str2double(line{4});
end

fclose(fid);

% search neighbors, delete overlap particles and assign new type
nneighbors = 36;
nneighbors1 = 20;
nneighbors2 = 16;
radius = 1.43; % Al radius, obtained from ovito particle property
cutoff1 = 2*radius;
cutoff2 = 2*sqrt(2)*radius;

max_type = max(type);
neighbor = {};

for i=1:nparticle
    for j=1:nparticle
        dis = sqrt((pos(i,1)-pos(j,1))^2 + (pos(i,2)-pos(j,2))^2 + (pos(i,3)-pos(j,3))^2);
        if((dis < 1.01*cutoff2) && (i ~= j))
            layer=1;
            if(dis < 1.01*cutoff1)
                layer=0;
            end
            neighbor = [neighbor, [i, j, layer]];
        end
    end
end

function skipNRows(fid, n)
    for i=1:n
        fgetl(fid);
    end
end
