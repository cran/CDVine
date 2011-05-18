BiCopName <- function(family, short=TRUE)
{
if(is.logical(short)==FALSE) stop("'short' has to be a logical variable.")

if(is.numeric(family))	# Zahl zu Name
{
	if(short==TRUE)		# kurzer Name
	{
		if(family==0) fam="I"
		else if(family==1) fam="N"
		else if(family==2) fam="t"
		else if(family==3) fam="C"
		else if(family==4) fam="G"
		else if(family==5) fam="F"
		else if(family==6) fam="J"
		else if(family==7) fam="BB1"
		else if(family==9) fam="BB7"
		else if(family==13) fam="SC"
		else if(family==14) fam="SG"
		else if(family==16) fam="SJ"
		else if(family==23) fam="C90"
		else if(family==24) fam="G90"
		else if(family==26) fam="J90"
		else if(family==33) fam="C270"
		else if(family==34) fam="G270"
		else if(family==36) fam="J270"
		else stop("Family not implemented.")
	}
	else			# langer Name
	{
		if(family==0) fam="Independence"
		else if(family==1) fam="Gaussian"
		else if(family==2) fam="t"
		else if(family==3) fam="Clayton"
		else if(family==4) fam="Gumbel"
		else if(family==5) fam="Frank"
		else if(family==6) fam="Joe"
		else if(family==7) fam="Clayton-Gumbel"
		else if(family==9) fam="Joe-Clayton"
		else if(family==13) fam="Survival Clayton"
		else if(family==14) fam="Survival Gumbel"
		else if(family==16) fam="Survival Joe"
		else if(family==23) fam="Rotated Clayton 90 degrees"
		else if(family==24) fam="Rotated Gumbel 90 degrees"
		else if(family==26) fam="Rotated Joe 90 degrees"
		else if(family==33) fam="Rotated Clayton 270 degrees"
		else if(family==34) fam="Rotated Gumbel 270 degrees"
		else if(family==36) fam="Rotated Joe 270 degrees"
		else stop("Family not implemented.")
	}
}
else	# Name zu Zahl
{
	if(family=="I" || family=="Independence") fam=0
	else if(family=="N" || family=="Gaussian") fam=1
	else if(family=="t") fam=2
	else if(family=="C" || family=="Clayton") fam=3
	else if(family=="G" || family=="Gumbel") fam=4
	else if(family=="F" || family=="Frank") fam=5
	else if(family=="J" || family=="Joe") fam=6
	else if(family=="BB1" || family=="Clayton-Gumbel") fam=7
	else if(family=="BB7" || family=="Joe-Clayton") fam=9
	else if(family=="SC" || family=="Survival Clayton") fam=13
	else if(family=="SG" || family=="Survival Gumbel") fam=14
	else if(family=="SJ" || family=="Survival Joe") fam=16
	else if(family=="C90" || family=="Rotated Clayton 90 degrees") fam=23
	else if(family=="G90" || family=="Rotated Gumbel 90 degrees") fam=24
	else if(family=="J90" || family=="Rotated Joe 90 degrees") fam=26
	else if(family=="C270" || family=="Rotated Clayton 270 degrees") fam=33
	else if(family=="G270" || family=="Rotated Gumbel 270 degrees") fam=34
	else if(family=="J270" || family=="Rotated Joe 270 degrees") fam=36
	else stop("Family not implemented.")
}

return(fam)
}