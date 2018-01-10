if ( .Platform$OS.type == 'windows' ) memory.limit( 256000 )

library(lodown)
this_sample_break <- Sys.getenv( "this_sample_break" )
timss_cat <- get_catalog( "timss" , output_dir = file.path( getwd() ) )
record_categories <- ceiling( seq( nrow( timss_cat ) ) / ceiling( nrow( timss_cat ) / 6 ) )
timss_cat <- timss_cat[ record_categories == this_sample_break , ]
lodown( "timss" , timss_cat )
if( any( timss_cat$year == 2015 ) ){
library(lodown)
# examine all available TIMSS microdata files
timss_cat <-
	get_catalog( "timss" ,
		output_dir = file.path( getwd() ) )

# 2015 only
timss_cat <- subset( timss_cat , year == 2015 )
# download the microdata to your local computer


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
lodown:::timss_MIcombine( with( timss_design , svyby( ~ one , ~ one , unwtd.count ) ) )

lodown:::timss_MIcombine( with( timss_design , svyby( ~ one , ~ sex , unwtd.count ) ) )
lodown:::timss_MIcombine( with( timss_design , svytotal( ~ one ) ) )

lodown:::timss_MIcombine( with( timss_design ,
	svyby( ~ one , ~ sex , svytotal )
) )
lodown:::timss_MIcombine( with( timss_design , svymean( ~ asmmat ) ) )

lodown:::timss_MIcombine( with( timss_design ,
	svyby( ~ asmmat , ~ sex , svymean )
) )
lodown:::timss_MIcombine( with( timss_design , svymean( ~ idcntry ) ) )

lodown:::timss_MIcombine( with( timss_design ,
	svyby( ~ idcntry , ~ sex , svymean )
) )
lodown:::timss_MIcombine( with( timss_design , svytotal( ~ asmmat ) ) )

lodown:::timss_MIcombine( with( timss_design ,
	svyby( ~ asmmat , ~ sex , svytotal )
) )
lodown:::timss_MIcombine( with( timss_design , svytotal( ~ idcntry ) ) )

lodown:::timss_MIcombine( with( timss_design ,
	svyby( ~ idcntry , ~ sex , svytotal )
) )
lodown:::timss_MIcombine( with( timss_design , svyquantile( ~ asmmat , 0.5 , se = TRUE ) ) )

lodown:::timss_MIcombine( with( timss_design ,
	svyby( 
		~ asmmat , ~ sex , svyquantile , 0.5 ,
		se = TRUE , keep.var = TRUE , ci = TRUE 
) ) )
lodown:::timss_MIcombine( with( timss_design ,
	svyratio( numerator = ~ asssci , denominator = ~ asmmat )
) )
sub_timss_design <- subset( timss_design , idcntry %in% c( 36 , 40 , 31 , 957 ) )
lodown:::timss_MIcombine( with( sub_timss_design , svymean( ~ asmmat ) ) )
this_result <-
	lodown:::timss_MIcombine( with( timss_design ,
		svymean( ~ asmmat )
	) )

coef( this_result )
SE( this_result )
confint( this_result )
cv( this_result )

grouped_result <-
	lodown:::timss_MIcombine( with( timss_design ,
		svyby( ~ asmmat , ~ sex , svymean )
	) )

coef( grouped_result )
SE( grouped_result )
confint( grouped_result )
cv( grouped_result )
degf( timss_design$designs[[1]] )
lodown:::timss_MIcombine( with( timss_design , svyvar( ~ asmmat ) ) )
# SRS without replacement
lodown:::timss_MIcombine( with( timss_design ,
	svymean( ~ asmmat , deff = TRUE )
) )

# SRS with replacement
lodown:::timss_MIcombine( with( timss_design ,
	svymean( ~ asmmat , deff = "replace" )
) )
lodown:::MIsvyciprop( ~ born_2001_or_later , timss_design ,
	method = "likelihood" , na.rm = TRUE )
lodown:::MIsvyttest( asmmat ~ born_2001_or_later , timss_design )
lodown:::MIsvychisq( ~ born_2001_or_later + idcntry , timss_design )
glm_result <- 
	lodown:::timss_MIcombine( with( timss_design ,
		svyglm( asmmat ~ born_2001_or_later + idcntry )
	) )
	
summary( glm_result )

}
