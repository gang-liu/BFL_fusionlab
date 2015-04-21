function dice = computeOverallDice(labelimg1, labelimg2)

%assumes foreground > 0

overlap = sum(sum(sum((labelimg1> 0) & (labelimg2 > 0))));
union = (sum(sum(sum(labelimg1> 0 | labelimg2> 0))));
    
dice = overlap/union;

