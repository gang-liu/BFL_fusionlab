function [dice, dice_global, dice_false_pos, dice_false_neg] = computeLabelDice(labelimg1, labelimg2, labels)
%function [dice, dice_global, false_pos, false_neg] = computeLabelDice(labelimg1, labelimg2, labels)
%assumes first one, i.e. labelimg1, is ground truth

if (nargin < 3)
    labels = unique(labelimg1(:));
end

dice = zeros(length(labels),1);
dice_global = 0;
dice_false_pos = dice;
dice_false_neg = dice;

for i = 1:length(labels)
    
    overlap = sum(sum(sum((labelimg1 == labels(i)).*(labelimg2 == labels(i)))));
    
    size1 = sum(sum(sum(labelimg1 == labels(i))));
    size2 = sum(sum(sum(labelimg2 == labels(i))));
    union = (size1 + size2)/2;
    
    dice(i) = overlap/union;
    dice_global = dice_global + overlap;
    
    if (nargout > 2)
       
        false_pos = sum(sum(sum((labelimg1 ~= labels(i)).*(labelimg2 == labels(i)))));
        false_neg = sum(sum(sum((labelimg1 == labels(i)).*(labelimg2 ~= labels(i)))));
        
        dice_false_pos(i) = false_pos/size2;
        dice_false_neg(i) = false_neg/size1;
    end
end

dice_global = dice_global/numel(labelimg1);