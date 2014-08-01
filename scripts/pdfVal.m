% Take in a PDF from KSDENSITY() and interpolate to give an estimate for a
% given value. When using a real distribution (weibull, rayleigh, etc.) you
% have a continuous PDF equation. Since we are using KSDENSITY to get an
% unknown "function"... we don't have such a luxury.

function outval = pdfVal(f,x, inval)


dx = x(2)-x(1);

nearestInd = find(min(abs(x-inval)) == abs(x-inval));

if nearestInd == length(x) || nearestInd == 1
    outval = f(nearestInd);
else
    if abs(x(nearestInd-1)-inval) < abs(x(nearestInd+1)-inval)
        ind = [nearestInd-1 nearestInd];
    else
        ind = [nearestInd nearestInd+1];
    end
    
    slope = f(ind(2)) - f(ind(1));
    
    numDX = (inval-x(ind(1)))/dx;
    
    outval = (numDX * slope) + f(ind(1));
end

end