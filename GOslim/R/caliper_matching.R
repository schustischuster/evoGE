
set.seed(11)
a <- data.frame(value=rnorm(500,7), sign=1)
b <- data.frame(value=rnorm(500,8), sign=0)
df <- rbind(a,b)



		# Create background gene set with caliper matching

		calpr <- 1

        matchSample <- function(x) {

            success <- FALSE
            while (!success) {

                # Create background gene set
                match_res <- matchit(sign ~ value, x, method = "nearest", caliper = calpr, 
                	std.caliper = TRUE, distance = "logit", replace = FALSE, m.order = "data", 
                    ratio = 1)
                match_res_m <- match_res$match.matrix

                # Extract standard mean difference from matchIt summary data
                comp <- as.data.frame(summary(match_res, standardize = TRUE)["sum.matched"])
                stmdif <- abs(comp[1,3])
                varR <- abs(comp[1,4])

                calpr <- calpr-0.01

                # check for success
                success <- ((stmdif <= 0.01) && (varR >= 0.99))
            }

            return(match_res_m)
        }

        match_res_df <- matchSample(df)

        matched <- df[match_res_df,]

        boxplot(a[,1],b[,1],matched[,1])







