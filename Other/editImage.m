
clear all

a = imread('ines_PV1.png');
% imagesc(a)

b = a;

for y = 1:size(a,1)
    for x = 1:size(a,2)
        cols([1:3]) = a(y,x,:);
        if cols(1) == cols(2) & cols(1) == cols(3) & min(cols) < 255
            b(y,x,:) = 255;
        end
        if cols(1) < 100 & cols(2) < 100 & cols(3) < 100
            b(y,x,:) = 255;
        elseif cols(1) > 10 & cols(2) < 240 & cols(3) < 240
            if cols(1) < 240 & cols(3) > 10 
                b(y,x,:) = 255;
            end
        end
        if cols(1) > 220 & cols(2) > 220 & cols(3) > 220
            b(y,x,:) = 255;
        end
    end
end

figure
imagesc(b)
aesthetics
box off
axis off

%%

clear all

a = imread('ines_PV1.png');
% imagesc(a)

b = a;

for y = 1:size(a,1)
    for x = 1:size(a,2)
        cols([1:3]) = a(y,x,:);
        if cols(1) == cols(2) & cols(1) == cols(3) & min(cols) < 255
            b(y,x,:) = 255;
        end
        if cols(1) < 100 & cols(2) < 100 & cols(3) < 100
            b(y,x,:) = 255;
        elseif cols(1) > 10 & cols(2) < 240 & cols(3) < 240
            if cols(1) < 240 & cols(3) > 10 
                b(y,x,:) = 255;
            end
        end
        if cols(1) > 220 & cols(2) > 220 & cols(3) > 220
            b(y,x,:) = 255;
        end        
        if cols(1) < 150
            b(y,x,:) = 255;
        end
    end
end

figure
imagesc(b)
aesthetics
box off
axis off