#%%
module BinModule

"""
// ****************************************
// Get Bin Edge Info
//******************************************* 
"""


include("BinEdgeInfo.jl")

export  Getbininfo

function Getbininfo(n)

    if n < 1
        println("Error : Bin number n should be [1,153] ")
    elseif n > 153 
        println("Error : Bin number n should be [1,153] ")
    else
    println("Info Bin number = ",n)
    println("Q2 GeV2 : ")
    println(BinQ2low[n])
    println(BinQ2high[n])
    println("x : " )
    println(Binxlow[n])        
    println(Binxhigh[n])


end



end

end #end of module

#%%