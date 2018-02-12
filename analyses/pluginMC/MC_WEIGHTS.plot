# BEGIN PLOT /MC_WEIGHTS/.*weight.*
Title=Event weight distribution
XLabel=$w$
YLabel=$1/N dN/dw$
# END PLOT


# BEGIN PLOT /MC_XS/xsfraction_neg
Title=Negative weight fraction
XLabel=
YLabel=$\frac{\sum_{w_i<0} |w_i|}{\sum |w_i|}$
LogY=0
ShowZero=0
XMajorTickMarks=20
XMinorTickMarks=0
# END PLOT
