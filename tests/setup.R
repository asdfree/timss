# brando for stella,
# gump's jenny, rock's adrian,
# students toward math test
timss_MIcombine <-
	function (results, variances, call = sys.call(), df.complete = Inf, ...) {
		m <- length(results)
		oldcall <- attr(results, "call")
		if (missing(variances)) {
			variances <- suppressWarnings(lapply(results, vcov))
			results <- lapply(results, coef)
		}
		vbar <- variances[[1]]
		cbar <- results[[1]]
		for (i in 2:m) {
			cbar <- cbar + results[[i]]
			vbar <- vbar + variances[[i]]
		}
		cbar <- cbar/m
		vbar <- vbar/m

		# MODIFICATION
		# evar <- var(do.call("rbind", results))
		evar <- sum( ( unlist( results ) - cbar )^2 / 4 )

		
		r <- (1 + 1/m) * evar/vbar
		df <- (m - 1) * (1 + 1/r)^2
		if (is.matrix(df)) df <- diag(df)
		if (is.finite(df.complete)) {
			dfobs <- ((df.complete + 1)/(df.complete + 3)) * df.complete *
			vbar/(vbar + evar)
			if (is.matrix(dfobs)) dfobs <- diag(dfobs)
			df <- 1/(1/dfobs + 1/df)
		}
		if (is.matrix(r)) r <- diag(r)
		rval <- list(coefficients = cbar, variance = vbar + evar *
		(m + 1)/m, call = c(oldcall, call), nimp = m, df = df,
		missinfo = (r + 2/(df + 3))/(r + 1))
		class(rval) <- "MIresult"
		rval
	}
library(httr)

tf <- tempfile()

this_url <- "https://timss2019.org/international-database/downloads/T19_G4_SPSS%20Data.zip"

GET( this_url , write_disk( tf ) , progress() )

unzipped_files <- unzip( tf , exdir = tempdir() )
library(haven)

# limit unzipped files to those starting with `asg` followed by three letters followed by `m7`
asg_fns <- unzipped_files[ grepl( '^asg[a-z][a-z][a-z]m7' , basename( unzipped_files ) ) ]

# further limit asg files to the first ten countries
countries_thru_canada <- c("alb", "arm", "aus", "aut", "aze", "bhr", "bfl", "bih", "bgr", "can")

fns_thru_canada <- paste0( paste0( '^asg' , countries_thru_canada , 'm7' ) , collapse = "|" )

asg_alb_can_fns <- asg_fns[ grepl( fns_thru_canada , basename( asg_fns ) ) ]

timss_df <- NULL

for( spss_fn in asg_alb_can_fns ){

	this_tbl <- read_spss( spss_fn )
	
	this_tbl <- zap_labels( this_tbl )
	
	this_df <- data.frame( this_tbl )
	
	names( this_df ) <- tolower( names( this_df ) )
	
	timss_df <- rbind( timss_df , this_df )
	
}

# order the data.frame by unique student id
timss_df <- timss_df[ with( timss_df , order( idcntry , idstud ) ) , ]
# timss_fn <- file.path( path.expand( "~" ) , "TIMSS" , "this_file.rds" )
# saveRDS( timss_df , file = timss_fn , compress = FALSE )
# timss_df <- readRDS( timss_fn )
# identify all columns ending with `01` thru `05`
ppv <- grep( "(.*)0[1-5]$" , names( timss_df ) , value = TRUE )

# remove those ending digits
ppv_prefix <- gsub( "0[1-5]$" , "" , ppv )

# identify each of the possibilities with exactly five matches (five implicates)
pv <- names( table( ppv_prefix )[ table( ppv_prefix ) == 5 ] )

# identify each of the `01` thru `05` plausible value columns
pv_columns <-
	grep( 
		paste0( "^" , pv , "0[1-5]$" , collapse = "|" ) , 
		names( timss_df ) , 
		value = TRUE 
	)
pv_wide_df <- timss_df[ c( 'idcntry' , 'idstud' , pv_columns ) ]

timss_df[ pv_columns ] <- NULL
pv_long_df <- 
	reshape( 
		pv_wide_df , 
		varying = lapply( paste0( pv , '0' ) , paste0 , 1:5 ) , 
		direction = 'long' , 
		timevar = 'implicate' , 
		idvar = c( 'idcntry' , 'idstud' ) 
	)

names( pv_long_df ) <- gsub( "01$" , "" , names( pv_long_df ) )
timss_long_df <- merge( timss_df , pv_long_df )

timss_long_df <- timss_long_df[ with( timss_long_df , order( idcntry , idstud ) ) , ]

stopifnot( nrow( timss_long_df ) == nrow( pv_long_df ) )

stopifnot( nrow( timss_long_df ) / 5 == nrow( timss_df ) )
timss_list <- split( timss_long_df , timss_long_df[ , 'implicate' ] )
weights_df <- timss_df[ c( 'jkrep' , 'jkzone' ) ]

for( j in 1:75 ){
	for( i in 0:1 ){
		weights_df[ weights_df[ , 'jkzone' ] != j , paste0( 'rw' , i , j ) ] <- 1
		
		weights_df[ weights_df[ , 'jkzone' ] == j , paste0( 'rw' , i , j ) ] <- 
			2 * ( weights_df[ weights_df[ , 'jkzone' ] == j , 'jkrep' ] == i )
	}
}

weights_df[ c( 'jkrep' , 'jkzone' ) ] <- NULL

library(survey)
library(mitools)

timss_design <- 
	svrepdesign(
		weights = ~totwgt ,
		repweights = weights_df , 
		data = imputationList( timss_list ) ,
		type = "other" ,
		scale = 0.5 ,
		rscales = rep( 1 , 150 ) ,
		combined.weights = FALSE ,
		mse = TRUE
	)
timss_design <- 
	update( 
		timss_design , 
		
		one = 1 ,
		
		countries_thru_canada = 
		
			factor( 
			
				as.numeric( idcntry ) ,
				
				levels = c(8L, 51L, 36L, 40L, 31L, 48L, 956L, 70L, 100L, 124L) ,

				labels =
					c("Albania", "Armenia", "Australia", "Austria", "Azerbaijan", "Bahrain",
					"Belgium (Flemish)", "Bosnia and Herzegovina", "Bulgaria", "Canada")
				
			) ,
		
		sex = factor( asbg01 , levels = 1:2 , labels = c( "female" , "male" ) ) ,
		
		born_in_country = ifelse( asbg07 %in% 1:2 , as.numeric( asbg07 == 1 ) , NA )

	)
timss_MIcombine( with( timss_design , svyby( ~ one , ~ one , unwtd.count ) ) )

timss_MIcombine( with( timss_design , svyby( ~ one , ~ sex , unwtd.count ) ) )
timss_MIcombine( with( timss_design , svytotal( ~ one ) ) )

timss_MIcombine( with( timss_design ,
	svyby( ~ one , ~ sex , svytotal )
) )
timss_MIcombine( with( timss_design , svymean( ~ asmmat , na.rm = TRUE ) ) )

timss_MIcombine( with( timss_design ,
	svyby( ~ asmmat , ~ sex , svymean , na.rm = TRUE )
) )
timss_MIcombine( with( timss_design , svymean( ~ countries_thru_canada ) ) )

timss_MIcombine( with( timss_design ,
	svyby( ~ countries_thru_canada , ~ sex , svymean )
) )
timss_MIcombine( with( timss_design , svytotal( ~ asmmat , na.rm = TRUE ) ) )

timss_MIcombine( with( timss_design ,
	svyby( ~ asmmat , ~ sex , svytotal , na.rm = TRUE )
) )
timss_MIcombine( with( timss_design , svytotal( ~ countries_thru_canada ) ) )

timss_MIcombine( with( timss_design ,
	svyby( ~ countries_thru_canada , ~ sex , svytotal )
) )
timss_MIcombine( with( timss_design ,
	svyquantile(
		~ asmmat ,
		0.5 , se = TRUE , na.rm = TRUE 
) ) )

timss_MIcombine( with( timss_design ,
	svyby(
		~ asmmat , ~ sex , svyquantile ,
		0.5 , se = TRUE ,
		ci = TRUE , na.rm = TRUE
) ) )
timss_MIcombine( with( timss_design ,
	svyratio( numerator = ~ asssci , denominator = ~ asmmat )
) )
sub_timss_design <- subset( timss_design , idcntry %in% c( 36 , 40 , 31 , 956 ) )
timss_MIcombine( with( sub_timss_design , svymean( ~ asmmat , na.rm = TRUE ) ) )
this_result <-
	timss_MIcombine( with( timss_design ,
		svymean( ~ asmmat , na.rm = TRUE )
	) )

coef( this_result )
SE( this_result )
confint( this_result )
cv( this_result )

grouped_result <-
	timss_MIcombine( with( timss_design ,
		svyby( ~ asmmat , ~ sex , svymean , na.rm = TRUE )
	) )

coef( grouped_result )
SE( grouped_result )
confint( grouped_result )
cv( grouped_result )
degf( timss_design$designs[[1]] )
timss_MIcombine( with( timss_design , svyvar( ~ asmmat , na.rm = TRUE ) ) )
# SRS without replacement
timss_MIcombine( with( timss_design ,
	svymean( ~ asmmat , na.rm = TRUE , deff = TRUE )
) )

# SRS with replacement
timss_MIcombine( with( timss_design ,
	svymean( ~ asmmat , na.rm = TRUE , deff = "replace" )
) )
# MIsvyciprop( ~ born_in_country , timss_design ,
# 	method = "likelihood" , na.rm = TRUE )
# MIsvyttest( asmmat ~ born_in_country , timss_design )
# MIsvychisq( ~ born_in_country + countries_thru_canada , timss_design )
glm_result <- 
	timss_MIcombine( with( timss_design ,
		svyglm( asmmat ~ born_in_country + countries_thru_canada )
	) )
	
summary( glm_result )
australia_design <- subset( timss_design , countries_thru_canada %in% "Australia" )

stopifnot( nrow( australia_design ) == 5890 )

result <- timss_MIcombine( with( australia_design , svymean( ~ asmmat ) ) )

stopifnot( round( coef( result ) , 3 ) == 515.880 )

stopifnot( round( SE( result ) , 3 ) == 2.776 )

australia_fn <- unzipped_files[ grepl( 'asgaus' , basename( unzipped_files ) ) ]
australia_tbl <- read_spss( australia_fn )
australia_tbl <- zap_labels( australia_tbl )
australia_df <- data.frame( australia_tbl )
names( australia_df ) <- tolower( names( australia_df ) )

mean_proficiency <-
	mean( c(
		with( australia_df , weighted.mean( asmmat01 , totwgt ) ) ,
		with( australia_df , weighted.mean( asmmat02 , totwgt ) ) ,
		with( australia_df , weighted.mean( asmmat03 , totwgt ) ) ,
		with( australia_df , weighted.mean( asmmat04 , totwgt ) ) ,
		with( australia_df , weighted.mean( asmmat05 , totwgt ) )
	) )

stopifnot( round( mean_proficiency , 3 ) == 515.880 )

for( k in 1:5 ){

	this_variance <- 0
	
	for( j in 1:75 ){
		for( i in 0:1 ){
			this_variance <- 
				this_variance + 
				( 
					weighted.mean( 
						australia_df[ , paste0( 'asmmat0' , k ) ] , 
						ifelse( 
							j == australia_df[ , 'jkzone' ] , 
							australia_df[ , 'totwgt' ] * 2 * ( australia_df[ , 'jkrep' ] == i ) , 
							australia_df[ , 'totwgt' ] 
						)
					) -
					weighted.mean( 
						australia_df[ , paste0( 'asmmat0' , k ) ] , 
						australia_df[ , 'totwgt' ]
					)
				)^2
		}
	}
	
	assign( paste0( 'v' , k ) , this_variance * 0.5 )

}

sampling_variance <- mean( c( v1 , v2 , v3 , v4 , v5 ) )
stopifnot( round( sampling_variance , 3 ) == 7.397 )

t0 <-
	mean( c(
		with( australia_df , weighted.mean( asmmat01 , totwgt ) ) ,
		with( australia_df , weighted.mean( asmmat02 , totwgt ) ) ,
		with( australia_df , weighted.mean( asmmat03 , totwgt ) ) ,
		with( australia_df , weighted.mean( asmmat04 , totwgt ) ) ,
		with( australia_df , weighted.mean( asmmat05 , totwgt ) )
	) )

imputation_variance <-
	( 6 / 5 ) * 
		( 
			( ( with( australia_df , weighted.mean( asmmat01 , totwgt ) ) - t0 )^2 / 4 ) +
			( ( with( australia_df , weighted.mean( asmmat02 , totwgt ) ) - t0 )^2 / 4 ) +
			( ( with( australia_df , weighted.mean( asmmat03 , totwgt ) ) - t0 )^2 / 4 ) +
			( ( with( australia_df , weighted.mean( asmmat04 , totwgt ) ) - t0 )^2 / 4 ) +
			( ( with( australia_df , weighted.mean( asmmat05 , totwgt ) ) - t0 )^2 / 4 ) 
		)

stopifnot( round( imputation_variance , 3 ) == 0.309 )

stopifnot( round( sampling_variance + imputation_variance , 3 ) == 7.706 )

