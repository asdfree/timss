if ( .Platform$OS.type == 'windows' ) memory.limit( 256000 )

options("lodown.cachaca.savecache"=FALSE)

library(lodown)
this_sample_break <- Sys.getenv( "this_sample_break" )
timss_cat <- get_catalog( "timss" , output_dir = file.path( getwd() ) )
record_categories <- ceiling( seq( nrow( timss_cat ) ) / ceiling( nrow( timss_cat ) / 6 ) )
timss_cat <- timss_cat[ record_categories == this_sample_break , ]
timss_cat <- lodown( "timss" , timss_cat )
if( any( timss_cat$year == 2015 ) ){











library(survey)
library(mitools)
library(RSQLite)

# load the ASG (student background) + ASH (home background) merged design
timss_design <- readRDS( file.path( getwd() , "2015/asg_design.rds" ) )

design_weights <- readRDS( file.path( getwd() , "2015/asg_weights.rds" ) )

five_tablenames <- paste0( "asg_2015_" , 1:5 )

timss_design <- lodown:::svyMDBdesign( timss_design )
timss_design <- 
	update( 
		timss_design , 
		
		one = 1 ,
		
		idcntry = factor( idcntry ) ,
		
		sex = factor( itsex , labels = c( "male" , "female" ) ) ,
		
		born_2005_or_later = as.numeric( itbirthy >= 2005 )

	)
timss_MIcombine( with( timss_design , svyby( ~ one , ~ one , unwtd.count ) ) )

timss_MIcombine( with( timss_design , svyby( ~ one , ~ sex , unwtd.count ) ) )
timss_MIcombine( with( timss_design , svytotal( ~ one ) ) )

timss_MIcombine( with( timss_design ,
	svyby( ~ one , ~ sex , svytotal )
) )
timss_MIcombine( with( timss_design , svymean( ~ asmmat ) ) )

timss_MIcombine( with( timss_design ,
	svyby( ~ asmmat , ~ sex , svymean )
) )
timss_MIcombine( with( timss_design , svymean( ~ idcntry ) ) )

timss_MIcombine( with( timss_design ,
	svyby( ~ idcntry , ~ sex , svymean )
) )
timss_MIcombine( with( timss_design , svytotal( ~ asmmat ) ) )

timss_MIcombine( with( timss_design ,
	svyby( ~ asmmat , ~ sex , svytotal )
) )
timss_MIcombine( with( timss_design , svytotal( ~ idcntry ) ) )

timss_MIcombine( with( timss_design ,
	svyby( ~ idcntry , ~ sex , svytotal )
) )
timss_MIcombine( with( timss_design , svyquantile( ~ asmmat , 0.5 , se = TRUE ) ) )

timss_MIcombine( with( timss_design ,
	svyby( 
		~ asmmat , ~ sex , svyquantile , 0.5 ,
		se = TRUE , keep.var = TRUE , ci = TRUE 
) ) )
timss_MIcombine( with( timss_design ,
	svyratio( numerator = ~ asssci , denominator = ~ asmmat )
) )
sub_timss_design <- subset( timss_design , idcntry %in% c( 36 , 40 , 31 , 957 ) )
timss_MIcombine( with( sub_timss_design , svymean( ~ asmmat ) ) )
this_result <-
	timss_MIcombine( with( timss_design ,
		svymean( ~ asmmat )
	) )

coef( this_result )
SE( this_result )
confint( this_result )
cv( this_result )

grouped_result <-
	timss_MIcombine( with( timss_design ,
		svyby( ~ asmmat , ~ sex , svymean )
	) )

coef( grouped_result )
SE( grouped_result )
confint( grouped_result )
cv( grouped_result )
degf( timss_design$designs[[1]] )
timss_MIcombine( with( timss_design , svyvar( ~ asmmat ) ) )
# SRS without replacement
timss_MIcombine( with( timss_design ,
	svymean( ~ asmmat , deff = TRUE )
) )

# SRS with replacement
timss_MIcombine( with( timss_design ,
	svymean( ~ asmmat , deff = "replace" )
) )
MIsvyciprop( ~ born_2001_or_later , timss_design ,
	method = "likelihood" , na.rm = TRUE )
MIsvyttest( asmmat ~ born_2001_or_later , timss_design )
MIsvychisq( ~ born_2001_or_later + idcntry , timss_design )
glm_result <- 
	timss_MIcombine( with( timss_design ,
		svyglm( asmmat ~ born_2001_or_later + idcntry )
	) )
	
summary( glm_result )

}
