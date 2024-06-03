function _insert_decimal(val::Int,digits)
    str = lpad(string(val),digits,"0")
    pos = length(str) - digits
    int = rpad(str[1:pos],1,"0")
    dec = str[pos+1:end]
    inserted = int*"."*dec
    return inserted
end
function _errorstring(x,Δx;nsig=2,force_decimal=false)
    @assert Δx > 0
    sgn = x < 0 ? "-" : ""
    x = abs(x)
    # For values larger 10^nsig, increase nsig to the number of integer digits
    digits_x  = ndigits(Int(round(x)))
    digits_Δx = ndigits(Int(round(Δx)))
    if max(digits_x,digits_Δx) > nsig
        nsig = max(digits_x,digits_Δx)
    end
    # round error part to desired number of signficant digits
    # convert to integer if no fractional part exists
    Δx_rounded = round(Δx,sigdigits=nsig) 
    # get number of decimal digits for x  
    floor_log_10 = floor(Int,log10(Δx))
    dec_digits   = (nsig - 1) - floor_log_10
    # round x, to desired number of decimal digits 
    # (standard julia function deals with negative dec_digits) 
    x_rounded = round(x,digits=dec_digits)
    # get decimal and integer part if there is a decimal part
    if dec_digits > 0
        digits_val = Int(round(x_rounded*10.0^(dec_digits)))
        digits_unc = Int(round(Δx_rounded*10.0^(dec_digits)))
        str_val = _insert_decimal(digits_val,dec_digits) 
        str_unc = _insert_decimal(digits_unc,dec_digits)
        str_unc = (nsig > dec_digits || force_decimal) ? str_unc : string(digits_unc)
        return sgn*str_val,str_unc
    else
        return sgn*"$(Int(x_rounded))","$(Int(Δx_rounded))"
    end
end
function _errorstring_asymmetric(x,Δxl,Δxu;nsig=2)
    # if one of the uncertainties is > 1, and the other < 1, make the decimal point explicit
    force_l = (Δxl < 1) && (Δxu > 1)
    force_u = (Δxl > 1) && (Δxu < 1)
    # get mean value and formatted uncertainties
    val1,parl = _errorstring(x, Δxl;nsig,force_decimal=force_l)
    val2,paru = _errorstring(x, Δxu;nsig,force_decimal=force_u)
    # check if we use the same number of significant digits in both cases
    sig1 = length(last(split(val1,".")))
    sig2 = length(last(split(val2,".")))
    # if not then adjust it so that we include more sig digits where needed
    if sig1 != sig2
        if sig1 < sig2
            val1,parl = _errorstring(x, Δxl;force_decimal=force_l,nsig=nsig+sig2-sig1)
            val2,paru = _errorstring(x, Δxu;force_decimal=force_u,nsig)
        else
            val1,parl = _errorstring(x, Δxl;force_decimal=force_l,nsig)
            val2,paru = _errorstring(x, Δxu;force_decimal=force_u,nsig=nsig+sig1-sig2)
        end
        sig1 = length(last(split(val1,".")))
        sig2 = length(last(split(val2,".")))    
    end
    @assert val1 == val2
    return val1, paru, parl 
end
function errorstring_asymmetric(x,Δxl,Δxu)
    val1, paru, parl = _errorstring_asymmetric(x,Δxl,Δxu)
    return "$val1(+$paru)(-$parl)"
end
function errorstring(x,Δx;nsig=2)
    val, par = _errorstring(x,Δx;nsig)
    return "$val($par)"
end
errorstring(x,Δxl,Δxu) = errorstring_asymmetric(x,Δxl,Δxu)