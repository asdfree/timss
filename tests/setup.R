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

timss_df <- NULL

for( spss_fn in unzipped_files[ grepl( '^asg[a-z][a-z][a-z]m7' , basename( unzipped_files ) ) ] ){

	this_tbl <- read_spss( spss_fn )
	
	this_tbl <- zap_labels( this_tbl )
	
	this_df <- data.frame( this_tbl )
	
	names( this_df ) <- tolower( names( this_df ) )
	
	timss_df <- rbind( timss_df , this_df )
	
}

timss_df <- timss_df[ with( timss_df , order( idcntry , idstud ) ) , ]
# timss_fn <- file.path( path.expand( "~" ) , "TIMSS" , "this_file.rds" )
# saveRDS( timss_df , file = timss_fn , compress = FALSE )
# timss_df <- readRDS( timss_fn )
ppv <- grep( "(.*)0[1-5]$" , names( timss_df ) , value = TRUE )
ppv_prefix <- gsub( "0[1-5]$" , "" , ppv )
pv <- names( table( ppv_prefix )[ table( ppv_prefix ) == 5 ] )
pv_columns <- grep( paste0( "^" , pv , "0[1-5]$" , collapse = "|" ) , names( timss_df ) , value = TRUE )

id_columns <- c( 'idcntry' , 'idstud' )

pv_wide_df <- timss_df[ c( id_columns , pv_columns ) ]

timss_df[ pv_columns ] <- NULL

pv_long_df <- 
	reshape( 
		pv_wide_df , 
		varying = lapply( paste0( pv , '0' ) , paste0 , 1:5 ) , 
		direction = 'long' , 
		timevar = 'implicate' , 
		idvar = id_columns 
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
timss_MIcombine( with( timss_design ,
	svyquantile(
		~ asmmat ,
		0.5 , se = TRUE 
) ) )

timss_MIcombine( with( timss_design ,
	svyby(
		~ asmmat , ~ sex , svyquantile ,
		0.5 , se = TRUE ,
		ci = TRUE 
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
# MIsvyciprop( ~ born_2001_or_later , timss_design ,
# 	method = "likelihood" , na.rm = TRUE )
# MIsvyttest( asmmat ~ born_2001_or_later , timss_design )
# MIsvychisq( ~ born_2001_or_later + idcntry , timss_design )
glm_result <- 
	timss_MIcombine( with( timss_design ,
		svyglm( asmmat ~ born_2001_or_later + idcntry )
	) )
	
summary( glm_result )
timss_design <- update( timss_design , one = 1 )
w <- subset( timss_design , idcntry == 36 )
nrow(w)
lapply( w$designs , function( u ) vcov( svymean( ~ asmmat , u ) ) )
mean( unlist( lapply( w$designs , function( u ) vcov( svymean( ~ asmmat , u ) ) ) ) )
timss_MIcombine( with( w , svymean( ~ asmmat ) ) )

library(haven)

this_df <- read_spss( grep( 'asgaus' , unzipped_files , value = TRUE ) )
this_df <- zap_labels( this_df )
this_df <- data.frame( this_df )
names( this_df ) <- tolower( names( this_df ) )

mean( c(
	with( this_df , weighted.mean( asmmat01 , totwgt ) ) ,
	with( this_df , weighted.mean( asmmat02 , totwgt ) ) ,
	with( this_df , weighted.mean( asmmat03 , totwgt ) ) ,
	with( this_df , weighted.mean( asmmat04 , totwgt ) ) ,
	with( this_df , weighted.mean( asmmat05 , totwgt ) )
) )

for( k in 1:5 ){

	this_variance <- 0
	
	for( j in 1:75 ){
		for( i in 0:1 ){
			this_variance <- 
				this_variance + 
				( 
					weighted.mean( 
						this_df[ , paste0( 'asmmat0' , k ) ] , 
						ifelse( j == this_df[ , 'jkzone' ] , this_df[ , 'totwgt' ] * 2 * ( this_df[ , 'jkrep' ] == i ) , this_df[ , 'totwgt' ] )
					) -
					weighted.mean( 
						this_df[ , paste0( 'asmmat0' , k ) ] , 
						this_df[ , 'totwgt' ]
					)
				)^2
		}
	}
	
	assign( paste0( 'v' , k ) , this_variance * 0.5 )

}

# sampling variance
mean( c( v1 , v2 , v3 , v4 , v5 ) )

t0 <-
	mean( c(
		with( this_df , weighted.mean( asmmat01 , totwgt ) ) ,
		with( this_df , weighted.mean( asmmat02 , totwgt ) ) ,
		with( this_df , weighted.mean( asmmat03 , totwgt ) ) ,
		with( this_df , weighted.mean( asmmat04 , totwgt ) ) ,
		with( this_df , weighted.mean( asmmat05 , totwgt ) )
	) )

# imputation variance
( 6 / 5 ) * 
	( 
		( ( with( this_df , weighted.mean( asmmat01 , totwgt ) ) - t0 )^2 / 4 ) +
		( ( with( this_df , weighted.mean( asmmat02 , totwgt ) ) - t0 )^2 / 4 ) +
		( ( with( this_df , weighted.mean( asmmat03 , totwgt ) ) - t0 )^2 / 4 ) +
		( ( with( this_df , weighted.mean( asmmat04 , totwgt ) ) - t0 )^2 / 4 ) +
		( ( with( this_df , weighted.mean( asmmat05 , totwgt ) ) - t0 )^2 / 4 ) 
	)

